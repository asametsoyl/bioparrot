import { browser } from 'wxt/browser';
import { CLS } from '@/core/constants';
import { isLlmEngine } from '@/core/engines/meta';
import { getAllExperts, suggestExpertForPage } from '@/core/experts';
import { t } from '@/core/i18n';
import { isRtl, isSameLanguage, normalizeLang } from '@/core/languages';
import { markupToPlain, tidyMarkup } from '@/core/markup';
import { lookupSectionHeading } from '@/core/section-headings';
import { makeBatches } from '@/core/batch';
import { TranslateError, type DisplayMode, type Settings, type TranslationContext } from '@/core/types';
import { sendToBackground } from '@/lib/messaging';
import type { PageInfo } from './page-info';
import { buildFragment, createError, createLoading, createWrapper, plainFragment, wrapperClass } from './render';
import { Segmenter, insideOurs, isOurs, shouldTranslate, type RawUnit, type Serialized } from './segmenter';
import { GLOBAL_EXCLUDE, REFERENCE_SELECTORS, resolveSiteRule, safeSelector, type ResolvedRule } from './site-rules';
import { applyStyleVariables, injectPageStyle } from './styles';

type UnitStatus = 'idle' | 'queued' | 'loading' | 'done' | 'error' | 'skipped';

interface Unit {
  id: number;
  root: Element;
  parent: Element | ShadowRoot;
  nodes: Node[];
  whole: boolean;
  ser: Serialized;
  status: UnitStatus;
  /** Yükleniyor/hata/çeviri öğesi. */
  el: HTMLElement | null;
  translation?: string;
  rich?: boolean;
  inline: boolean;
  /** Fareyle tek tek çevrilen birim. */
  manual: boolean;
  /** "Yalnız çeviri" modunda gizlenen özgün düğümler (geri yüklemek için). */
  hidden?: { texts: { node: Text; data: string }[]; elements: Element[] };
}

const MAIN_EXCLUDE = [
  'nav',
  'aside',
  '[role="navigation"]',
  '[role="banner"]',
  '[role="contentinfo"]',
  '[role="search"]',
  '[role="menu"]',
  '[role="menubar"]',
  '[role="toolbar"]',
  '[role="tablist"]',
  'header:not(article header):not(main header):not([role="main"] header)',
  'footer:not(article footer):not(main footer)',
  '.breadcrumb',
  '.breadcrumbs',
  '.skip-link',
];

export interface TranslatorHost {
  settings(): Settings;
  pageInfo: PageInfo;
  onChange(): void;
  onError(err: TranslateError): void;
}

/** Web sayfasının iki dilli çevirisini yöneten ana sınıf. */
export class PageTranslator {
  active = false;
  private units = new Map<number, Unit>();
  private rootUnits = new WeakMap<Element, Unit[]>();
  private nodeUnit = new WeakMap<Node, Unit>();
  private ownText = new WeakMap<Text, string>();
  private queue: Unit[] = [];
  private flushTimer: ReturnType<typeof setTimeout> | null = null;
  private inflight = 0;
  private io: IntersectionObserver | null = null;
  private mo: MutationObserver | null = null;
  private observedShadows = new WeakSet<ShadowRoot>();
  private segmenter: Segmenter | null = null;
  private rule: ResolvedRule | null = null;
  private includeSelector: string | null = null;
  private nextId = 1;
  private epoch = 0;
  private dirtyContainers = new Set<Element | ShadowRoot>();
  private dirtySubtrees = new Set<Element>();
  private dirtyTimer: ReturnType<typeof setTimeout> | null = null;
  private termsPromise: Promise<{ source: string; target: string }[]> | null = null;
  private styleId = 'underline';
  /** Seçim çevirisiyle paragrafların altına eklenen bloklar. */
  private selectionWraps = new Set<HTMLElement>();
  expertId = 'general';

  constructor(private host: TranslatorHost) {}

  get mode(): DisplayMode {
    return this.host.settings().displayMode;
  }

  get pending(): number {
    return this.queue.length + this.inflight;
  }

  // ---- Kurulum ---------------------------------------------------------------

  private setup(): Segmenter {
    const s = this.host.settings();
    injectPageStyle(document);
    applyStyleVariables(s.style);
    this.styleId = s.style.id;
    this.expertId = this.resolveExpert(s);
    this.rule = resolveSiteRule(location.hostname, s.rules.custom);
    const exclude = [...GLOBAL_EXCLUDE, ...this.rule.exclude];
    if (s.page.scope === 'main') exclude.push(...MAIN_EXCLUDE);
    if (s.academic.skipReferences) exclude.push(...REFERENCE_SELECTORS);
    const opts = { exclude: safeSelector(exclude), skipCode: s.page.skipCode, skipReferences: s.academic.skipReferences };
    this.includeSelector = s.page.scope === 'main' && this.rule.include.length ? safeSelector(this.rule.include) : null;
    if (this.segmenter) this.segmenter.reset(opts);
    else this.segmenter = new Segmenter(opts);
    return this.segmenter;
  }

  resolveExpert(s: Settings = this.host.settings()): string {
    const engine = s.engines[s.engineId];
    if (!isLlmEngine(engine)) return s.expertId;
    if (s.academic.autoExpert && s.expertId === 'general') {
      const suggested = suggestExpertForPage(location.hostname, this.host.pageInfo.isPaper, s.academic.expertId);
      if (suggested && getAllExperts(s).some((e) => e.id === suggested)) return suggested;
    }
    return s.expertId;
  }

  private scanRoots(): (Element | ShadowRoot)[] {
    if (this.includeSelector) {
      const matches = Array.from(document.querySelectorAll(this.includeSelector));
      const outer = matches.filter((el) => !matches.some((other) => other !== el && other.contains(el)));
      if (outer.length) return outer;
    }
    return document.body ? [document.body] : [];
  }

  // ---- Başlat / durdur --------------------------------------------------------

  async start(): Promise<void> {
    if (this.active) return;
    this.active = true;
    const epoch = ++this.epoch;
    const segmenter = this.setup();
    const s = this.host.settings();

    if (s.page.lazy && 'IntersectionObserver' in window) {
      this.io = new IntersectionObserver((entries) => this.onIntersect(entries), { rootMargin: '400px 0px 800px 0px' });
    }
    this.mo = new MutationObserver((records) => this.onMutations(records));
    if (document.body) this.mo.observe(document.body, { childList: true, subtree: true, characterData: true });
    this.host.onChange();

    for (const root of this.scanRoots()) {
      await segmenter.scan(root, (units) => this.register(units), () => epoch !== this.epoch);
    }
    this.host.onChange();
  }

  stop(): void {
    if (!this.active && !this.units.size) return;
    this.active = false;
    this.epoch++;
    this.io?.disconnect();
    this.io = null;
    this.mo?.disconnect();
    this.mo = null;
    this.observedShadows = new WeakSet();
    if (this.flushTimer) clearTimeout(this.flushTimer);
    if (this.dirtyTimer) clearTimeout(this.dirtyTimer);
    this.flushTimer = this.dirtyTimer = null;
    this.queue = [];
    this.dirtyContainers.clear();
    this.dirtySubtrees.clear();
    for (const unit of this.units.values()) this.removeUnitDom(unit);
    for (const wrap of this.selectionWraps) wrap.remove();
    this.selectionWraps.clear();
    this.units.clear();
    this.rootUnits = new WeakMap();
    this.nodeUnit = new WeakMap();
    this.termsPromise = null;
    this.host.onChange();
  }

  toggle(): void {
    if (this.active) this.stop();
    else void this.start();
  }

  /** Ayarlar değişince görünümü günceller (stil, mod) veya gerekirse yeniden çevirir. */
  onSettingsChanged(prev: Settings, next: Settings): void {
    applyStyleVariables(next.style);
    if (prev.style.id !== next.style.id) {
      this.styleId = next.style.id;
      for (const u of this.units.values()) {
        if (u.status === 'done' && u.el) u.el.className = wrapperClass(u.inline, this.styleId);
      }
    }
    if (prev.displayMode !== next.displayMode) {
      for (const u of this.units.values()) this.applyMode(u);
    }
    const needsRetranslate =
      prev.targetLanguage !== next.targetLanguage ||
      prev.sourceLanguage !== next.sourceLanguage ||
      prev.engineId !== next.engineId ||
      prev.expertId !== next.expertId ||
      prev.page.scope !== next.page.scope ||
      prev.page.skipCode !== next.page.skipCode ||
      prev.academic.skipReferences !== next.academic.skipReferences;
    if (needsRetranslate && this.active) {
      this.stop();
      void this.start();
    }
  }

  // ---- Kayıt ve kuyruk --------------------------------------------------------

  private register(raws: RawUnit[], manual = false): Unit[] {
    const segmenter = this.segmenter!;
    const created: Unit[] = [];
    for (const raw of raws) {
      const existing = this.nodeUnit.get(raw.nodes[0]);
      // "Yalnız çeviri" modunda özgün metni boşaltılmış birimler yeniden taranmaz.
      if (existing?.hidden) continue;
      const ser = segmenter.serialize(raw.nodes);
      if (existing && existing.ser.markup === ser.markup) continue;
      if (existing) this.removeUnit(existing);
      if (!shouldTranslate(ser.prose, ser.plain)) continue;
      const unit: Unit = {
        id: this.nextId++,
        root: raw.root,
        parent: raw.parent,
        nodes: raw.nodes,
        whole: raw.whole,
        ser,
        status: 'idle',
        el: null,
        inline: this.decideInline(raw, ser),
        manual,
      };
      this.units.set(unit.id, unit);
      for (const n of raw.nodes) this.nodeUnit.set(n, unit);
      const list = this.rootUnits.get(raw.root);
      if (list) list.push(unit);
      else this.rootUnits.set(raw.root, [unit]);
      if (raw.parent instanceof ShadowRoot) this.watchShadow(raw.parent);
      created.push(unit);

      if (manual || !this.io || getComputedStyle(raw.root).display === 'contents') this.enqueue(unit);
      else this.io.observe(raw.root);
    }
    return created;
  }

  private watchShadow(root: ShadowRoot): void {
    if (this.observedShadows.has(root)) return;
    this.observedShadows.add(root);
    injectPageStyle(root);
    this.mo?.observe(root, { childList: true, subtree: true, characterData: true });
  }

  private decideInline(raw: RawUnit, ser: Serialized): boolean {
    if (!this.host.settings().page.inlineShort) return false;
    if (/^H[1-6]$/.test(raw.root.tagName)) return false;
    const text = ser.plain.trim();
    if (text.length > 28 || text.includes('\n')) return false;
    return text.split(/\s+/).length <= 4;
  }

  private onIntersect(entries: IntersectionObserverEntry[]): void {
    for (const entry of entries) {
      if (!entry.isIntersecting) continue;
      this.io?.unobserve(entry.target);
      for (const unit of this.rootUnits.get(entry.target) ?? []) {
        if (unit.status === 'idle') this.enqueue(unit);
      }
    }
  }

  private enqueue(unit: Unit): void {
    unit.status = 'queued';
    this.showLoading(unit);
    this.queue.push(unit);
    if (!this.flushTimer) this.flushTimer = setTimeout(() => void this.flush(), 60);
    this.host.onChange();
  }

  private async flush(): Promise<void> {
    this.flushTimer = null;
    const units = this.queue.splice(0);
    if (!units.length) return;
    const s = this.host.settings();
    const epoch = this.epoch;
    const target = s.targetLanguage;

    // Zaten hedef dildeki (veya çevrilmesi istenmeyen dillerdeki) paragrafları atla.
    await Promise.all(
      units.map(async (u) => {
        if (u.ser.prose.length < 16) return;
        const { lang, reliable } = await detectLanguage(u.ser.prose);
        if (!reliable || !lang) return;
        if (isSameLanguage(lang, target) || s.rules.neverLanguages.some((l) => isSameLanguage(l, lang))) {
          u.status = 'skipped';
          this.removeIndicator(u);
        }
      }),
    );
    if (epoch !== this.epoch) return;

    const todo = units.filter((u) => u.status === 'queued' && u.nodes.some((n) => n.isConnected));
    for (const u of units) if (u.status === 'queued' && !todo.includes(u)) this.removeUnit(u);
    if (!todo.length) {
      this.host.onChange();
      return;
    }
    const engine = s.engines[s.engineId];
    // Standart makale bölüm başlıkları klasik motorlara gönderilmeden doğru karşılıkla gösterilir.
    if (!isLlmEngine(engine)) {
      for (const u of todo) {
        const heading = u.ser.elements.length === 0 ? lookupSectionHeading(u.ser.plain, target) : undefined;
        if (heading) this.render(u, heading, false);
      }
    }
    const pendingUnits = todo.filter((u) => u.status === 'queued');
    if (!pendingUnits.length) {
      this.host.onChange();
      return;
    }
    const context = await this.getContext(s);
    if (epoch !== this.epoch) return;
    const batches = makeBatches(
      pendingUnits.map((_, i) => i),
      pendingUnits.map((u) => u.ser.markup),
      engine?.batchSize ?? 20,
      engine?.batchChars ?? 3000,
    );
    await Promise.all(batches.map((b) => this.sendBatch(b.map((i) => pendingUnits[i]), context, epoch)));
  }

  private async sendBatch(units: Unit[], context: TranslationContext | undefined, epoch: number): Promise<void> {
    const s = this.host.settings();
    units.forEach((u) => (u.status = 'loading'));
    this.inflight += units.length;
    this.host.onChange();
    try {
      const res = await sendToBackground({
        type: 'translate',
        req: {
          texts: units.map((u) => u.ser.markup),
          from: s.sourceLanguage,
          to: s.targetLanguage,
          engineId: s.engineId,
          expertId: this.expertId,
          rich: true,
          context,
        },
      });
      if (epoch !== this.epoch) return;
      units.forEach((u, i) => {
        const translation = res.translations[i];
        const detected = res.detected?.[i];
        if (!translation || (detected && isSameLanguage(detected, s.targetLanguage))) {
          u.status = 'skipped';
          this.removeIndicator(u);
          return;
        }
        const plain = res.rich ? markupToPlain(translation) : translation;
        if (normalizeText(plain) === normalizeText(u.ser.plain)) {
          u.status = 'skipped';
          this.removeIndicator(u);
          return;
        }
        this.render(u, translation, res.rich);
      });
    } catch (err) {
      if (epoch !== this.epoch) return;
      const e = err instanceof TranslateError ? err : new TranslateError(String(err));
      units.forEach((u) => this.renderError(u, e));
      this.host.onError(e);
    } finally {
      this.inflight -= units.length;
      this.host.onChange();
    }
  }

  // ---- Bağlam ve terim tutarlılığı -------------------------------------------

  private async getContext(s: Settings): Promise<TranslationContext | undefined> {
    const engine = s.engines[s.engineId];
    const llm = isLlmEngine(engine);
    if (!s.page.contextAware || (!llm && engine?.type !== 'deepl')) return undefined;
    const info = this.host.pageInfo;
    const context: TranslationContext = { title: info.title, summary: info.summary, url: location.href };
    const academic = getAllExperts(s).find((e) => e.id === this.expertId)?.group === 'academic';
    if (llm && academic && s.academic.termConsistency && (info.summary || info.title)) {
      this.termsPromise ??= sendToBackground({
        type: 'extractTerms',
        text: [info.title, info.summary || this.sampleText()].filter(Boolean).join('\n\n'),
        to: s.targetLanguage,
        engineId: s.engineId,
        expertId: this.expertId,
      }).catch(() => []);
      const terms = await Promise.race([this.termsPromise, new Promise<[]>((r) => setTimeout(() => r([]), 12_000))]);
      if (terms.length) context.terms = terms;
    }
    return context;
  }

  private sampleText(): string {
    const roots = this.scanRoots();
    return roots
      .map((r) => (r instanceof Element ? (r as HTMLElement).innerText ?? '' : ''))
      .join('\n')
      .slice(0, 3000);
  }

  // ---- DOM işlemleri ----------------------------------------------------------

  private insertionPoint(unit: Unit): { parent: Node; before: Node | null } {
    if (unit.whole) return { parent: unit.root, before: null };
    const last = unit.nodes[unit.nodes.length - 1];
    return { parent: last.parentNode ?? unit.parent, before: last.nextSibling };
  }

  private place(unit: Unit, el: HTMLElement): void {
    if (unit.el && unit.el.isConnected) {
      unit.el.replaceWith(el);
    } else {
      const { parent, before } = this.insertionPoint(unit);
      parent.insertBefore(el, before && before.parentNode === parent ? before : null);
    }
    unit.el = el;
  }

  private showLoading(unit: Unit): void {
    this.place(unit, createLoading(document, unit.id));
  }

  private removeIndicator(unit: Unit): void {
    unit.el?.remove();
    unit.el = null;
  }

  private render(unit: Unit, translation: string, rich: boolean): void {
    const s = this.host.settings();
    const wrap = createWrapper(document, { inline: unit.inline, styleId: this.styleId, unitId: unit.id });
    const text = document.createElement('font');
    text.className = CLS.text;
    const fragment = rich ? buildFragment(tidyMarkup(translation), unit.ser.elements) : plainFragment(translation);
    capitalizeFirst(fragment, unit.ser.plain, s.targetLanguage);
    text.appendChild(fragment);
    if (unit.inline) {
      const sp = document.createElement('font');
      sp.className = 'bp-sp';
      sp.textContent = ' ';
      wrap.appendChild(sp);
    }
    wrap.appendChild(text);
    wrap.lang = s.targetLanguage;
    if (isRtl(s.targetLanguage)) wrap.dir = 'rtl';
    else if (getComputedStyle(unit.root).direction === 'rtl') wrap.dir = 'ltr';
    this.place(unit, wrap);
    unit.translation = translation;
    unit.rich = rich;
    unit.status = 'done';
    this.applyMode(unit);
  }

  private renderError(unit: Unit, err: TranslateError): void {
    const el = createError(document, unit.id, err.message, t('cnt.failed'));
    el.title = `${err.message}\n${t('cnt.clickRetry')}`;
    el.addEventListener('click', (ev) => {
      ev.preventDefault();
      ev.stopPropagation();
      unit.status = 'idle';
      this.enqueue(unit);
    });
    this.place(unit, el);
    unit.status = 'error';
  }

  private applyMode(unit: Unit): void {
    const translationOnly = this.mode === 'translation' && unit.status === 'done';
    if (unit.el) {
      unit.el.title = translationOnly ? `${t('cnt.original')}: ${unit.ser.plain}` : '';
      unit.el.classList.toggle('bp-flush', translationOnly);
    }
    if (translationOnly) this.hideOriginal(unit);
    else this.showOriginal(unit);
  }

  /** Özgün metni gizler: metin düğümleri boşaltılır, satır içi öğeler gizlenir (sayfa düzeni bozulmaz). */
  private hideOriginal(unit: Unit): void {
    if (unit.hidden) return;
    const hidden: NonNullable<Unit['hidden']> = { texts: [], elements: [] };
    for (const node of unit.nodes) {
      if (node.nodeType === Node.TEXT_NODE) {
        const text = node as Text;
        if (!text.data) continue;
        hidden.texts.push({ node: text, data: text.data });
        this.ownText.set(text, '');
        text.data = '';
      } else if (node.nodeType === Node.ELEMENT_NODE && !isOurs(node)) {
        (node as Element).classList.add(CLS.origHidden);
        hidden.elements.push(node as Element);
      }
    }
    unit.hidden = hidden;
  }

  private showOriginal(unit: Unit): void {
    const hidden = unit.hidden;
    if (!hidden) return;
    for (const { node, data } of hidden.texts) {
      if (node.data === '') {
        this.ownText.set(node, data);
        node.data = data;
      }
    }
    for (const el of hidden.elements) el.classList.remove(CLS.origHidden);
    unit.hidden = undefined;
  }

  private removeUnitDom(unit: Unit): void {
    unit.el?.remove();
    unit.el = null;
    this.showOriginal(unit);
  }

  private removeUnit(unit: Unit): void {
    this.removeUnitDom(unit);
    this.units.delete(unit.id);
    const list = this.rootUnits.get(unit.root);
    if (list) {
      const i = list.indexOf(unit);
      if (i >= 0) list.splice(i, 1);
    }
    for (const n of unit.nodes) if (this.nodeUnit.get(n) === unit) this.nodeUnit.delete(n);
    unit.status = 'skipped';
  }

  // ---- Dinamik içerik ---------------------------------------------------------

  private onMutations(records: MutationRecord[]): void {
    if (!this.active || !this.segmenter) return;
    const segmenter = this.segmenter;
    for (const r of records) {
      if (r.type === 'characterData') {
        const node = r.target as Text;
        if (this.ownText.get(node) === node.data) {
          this.ownText.delete(node);
          continue;
        }
        if (insideOurs(node)) continue;
        const unit = this.nodeUnit.get(node);
        const container = unit ? (unit.parent instanceof ShadowRoot ? unit.parent : unit.root) : segmenter.containerOf(node);
        if (container) this.dirtyContainers.add(container);
        continue;
      }
      if (r.type !== 'childList') continue;
      const target = r.target;
      if ((target.nodeType === Node.ELEMENT_NODE && isOurs(target)) || insideOurs(target)) continue;
      let relevant = false;
      r.addedNodes.forEach((n) => {
        if (isOurs(n)) return;
        relevant = true;
        if (n.nodeType === Node.ELEMENT_NODE) this.dirtySubtrees.add(n as Element);
      });
      r.removedNodes.forEach((n) => {
        if (!isOurs(n)) relevant = true;
      });
      if (!relevant) continue;
      if (target instanceof ShadowRoot) this.dirtyContainers.add(target);
      else if (target.nodeType === Node.ELEMENT_NODE) {
        segmenter.forget(target as Element);
        const container = segmenter.containerOf(target);
        if (container) this.dirtyContainers.add(container);
      }
    }
    if ((this.dirtyContainers.size || this.dirtySubtrees.size) && !this.dirtyTimer) {
      this.dirtyTimer = setTimeout(() => void this.processDirty(), 300);
    }
  }

  private async processDirty(): Promise<void> {
    this.dirtyTimer = null;
    if (!this.active || !this.segmenter) return;
    const segmenter = this.segmenter;
    const epoch = this.epoch;

    // DOM'dan çıkan birimleri temizle
    for (const unit of Array.from(this.units.values())) {
      if (!unit.nodes.some((n) => n.isConnected)) this.removeUnit(unit);
    }

    const containers = Array.from(this.dirtyContainers);
    const subtrees = Array.from(this.dirtySubtrees).filter((el) => el.isConnected);
    this.dirtyContainers.clear();
    this.dirtySubtrees.clear();

    for (const c of containers) {
      if (!(c instanceof ShadowRoot ? c.host.isConnected : c.isConnected)) continue;
      if (!this.inScope(c instanceof ShadowRoot ? c.host : c)) continue;
      this.rescanContainer(c);
    }
    const roots = subtrees.filter((el) => !subtrees.some((other) => other !== el && other.contains(el)));
    for (const el of roots) {
      if (epoch !== this.epoch) return;
      if (!this.inScope(el)) {
        // Kapsam dışı yeni öğe, kapsamdaki öğeleri içeriyor olabilir (örn. yeni yüklenen tweetler).
        if (this.includeSelector) {
          const inner = Array.from(el.querySelectorAll(this.includeSelector));
          for (const i of inner) await segmenter.scan(i, (u) => this.register(u), () => epoch !== this.epoch);
        }
        continue;
      }
      const container = segmenter.containerOf(el);
      if (container && container !== el) this.rescanContainer(container);
      await segmenter.scan(el, (u) => this.register(u), () => epoch !== this.epoch);
    }
  }

  private inScope(el: Element): boolean {
    if (!this.includeSelector) return true;
    try {
      return !!el.closest(this.includeSelector);
    } catch {
      return true;
    }
  }

  private rescanContainer(container: Element | ShadowRoot): void {
    const segmenter = this.segmenter!;
    const root = container instanceof ShadowRoot ? container.host : container;
    const existing = (this.rootUnits.get(root) ?? []).filter((u) => u.parent === container || container === root);
    const { units: fresh } = segmenter.runsOf(container);
    const keep = new Set<Unit>();
    const toRegister: RawUnit[] = [];
    for (const raw of fresh) {
      const owner = this.nodeUnit.get(raw.nodes[0]);
      if (owner?.hidden) {
        keep.add(owner);
        continue;
      }
      const ser = segmenter.serialize(raw.nodes);
      const match = existing.find((u) => !keep.has(u) && u.ser.markup === ser.markup);
      if (match) {
        keep.add(match);
        for (const n of match.nodes) if (this.nodeUnit.get(n) === match) this.nodeUnit.delete(n);
        match.nodes = raw.nodes;
        for (const n of raw.nodes) this.nodeUnit.set(n, match);
        continue;
      }
      toRegister.push(raw);
    }
    for (const u of existing) {
      if (keep.has(u) || u.manual) continue;
      if (u.hidden && u.nodes.some((n) => n.isConnected)) continue;
      this.removeUnit(u);
    }
    if (toRegister.length) this.register(toRegister);
  }

  // ---- Fareyle çeviri ---------------------------------------------------------

  /** Verilen öğeyi içeren paragrafı çevirir; zaten çevrilmişse kaldırır. */
  toggleAt(target: Element): void {
    if (!this.segmenter) this.setup();
    const segmenter = this.segmenter!;
    const ours = target.closest(`.${CLS.wrap}`);
    if (ours) {
      const id = Number(ours.getAttribute('data-bp-unit'));
      const unit = this.units.get(id);
      if (unit) this.removeUnit(unit);
      this.host.onChange();
      return;
    }
    const container = segmenter.containerOf(target);
    if (!container) return;
    const { units } = segmenter.runsOf(container);
    const raw = units.find((u) => u.nodes.some((n) => n === target || n.contains(target))) ?? units[0];
    if (!raw) return;
    const existing = this.nodeUnit.get(raw.nodes[0]);
    // Çeviri sarmalayıcısı sayfa tarafından kaldırılmışsa birim "çevrilmemiş" sayılır.
    const visible = !!existing?.el?.isConnected;
    if (existing && visible && (existing.status === 'done' || existing.status === 'loading' || existing.status === 'queued')) {
      this.removeUnit(existing);
      this.host.onChange();
      return;
    }
    if (existing) this.removeUnit(existing);
    if (!this.active && !this.mo) {
      injectPageStyle(document);
      applyStyleVariables(this.host.settings().style);
    }
    this.register([raw], true);
  }

  /**
   * Seçilen metnin çevirisini paragrafın altına ekler; aynı yerde tekrar çağrılırsa kaldırır.
   * Seçim paragrafın neredeyse tamamıysa paragraf biçimlendirmesi (bağlantılar, atıflar) korunarak çevrilir.
   */
  async translateSelection(range: Range, text: string): Promise<void> {
    if (!this.segmenter) this.setup();
    const segmenter = this.segmenter!;
    const s = this.host.settings();
    injectPageStyle(document);
    applyStyleVariables(s.style);
    this.styleId = s.style.id;
    // Üç tıklama gibi seçimler bir sonraki öğenin başında bitebilir; seçimdeki gerçek metin sınırları kullanılır.
    const bounds = rangeTextBounds(range);
    const endNode: Node = bounds.last ?? range.endContainer;
    const startNode: Node = bounds.first ?? range.startContainer;
    const anchor = endNode.nodeType === Node.ELEMENT_NODE ? (endNode as Element) : endNode.parentElement;
    if (!anchor || insideOurs(anchor)) return;
    const containerNode = segmenter.containerOf(anchor);
    if (!containerNode) return;
    const container = containerNode instanceof ShadowRoot ? containerNode.host : containerNode;

    // Aynı paragrafta önceki seçim çevirisi varsa kaldır (aç/kapat).
    for (const wrap of this.selectionWraps) {
      if (wrap.isConnected && container.contains(wrap)) {
        wrap.remove();
        this.selectionWraps.delete(wrap);
        if (wrap.dataset.bpSource === text) return;
      }
    }

    const { units } = segmenter.runsOf(containerNode);
    const contains = (n: Node, target: Node) => n === target || n.contains(target);
    const raw = units.find((u) => u.nodes.some((n) => contains(n, endNode)));
    if (raw) {
      const ser = segmenter.serialize(raw.nodes);
      const startsInside = raw.nodes.some((n) => contains(n, startNode));
      const selected = text.replace(/\s+/g, '').length;
      if (startsInside && selected >= ser.plain.replace(/\s+/g, '').length * 0.8) {
        this.toggleAt(anchor);
        return;
      }
    }

    const wrap = createWrapper(document, { inline: false, styleId: this.styleId, unitId: 0 });
    wrap.classList.add('bp-sel-wrap');
    wrap.dataset.bpSource = text;
    wrap.appendChild(createLoading(document, 0));
    if (raw && !raw.whole) {
      const last = raw.nodes[raw.nodes.length - 1];
      last.parentNode?.insertBefore(wrap, last.nextSibling);
    } else {
      container.appendChild(wrap);
    }
    this.selectionWraps.add(wrap);
    try {
      const res = await sendToBackground({
        type: 'translate',
        req: {
          texts: [text],
          from: s.sourceLanguage,
          to: s.targetLanguage,
          engineId: s.selection.engineId || s.engineId,
          expertId: this.expertId,
          rich: false,
          context: s.page.contextAware ? { title: this.host.pageInfo.title, summary: this.host.pageInfo.summary } : undefined,
        },
      });
      const out = document.createElement('font');
      out.className = CLS.text;
      out.appendChild(plainFragment(res.translations[0] ?? ''));
      wrap.replaceChildren(out);
      wrap.lang = s.targetLanguage;
      if (isRtl(s.targetLanguage)) wrap.dir = 'rtl';
    } catch (err) {
      const e = err instanceof TranslateError ? err : new TranslateError(String(err));
      wrap.replaceChildren(createError(document, 0, e.message, t('cnt.failed')));
      this.host.onError(e);
    }
  }

  /** Fareyle çevirinin hedefleyeceği paragrafın kapsayıcısı (vurgulamak için). */
  hoverTarget(target: Element): Element | null {
    if (!this.segmenter) this.setup();
    const c = this.segmenter!.containerOf(target);
    return c instanceof ShadowRoot ? c.host : c;
  }
}

function normalizeText(s: string): string {
  return s.replace(/[\s ]+/g, ' ').trim().toLowerCase();
}

const detectCache = new Map<string, { lang: string; reliable: boolean }>();

async function detectLanguage(text: string): Promise<{ lang: string; reliable: boolean }> {
  const sample = text.slice(0, 600);
  const cached = detectCache.get(sample);
  if (cached) return cached;
  let result = { lang: '', reliable: false };
  try {
    const r = await browser.i18n.detectLanguage(sample);
    const top = r?.languages?.[0];
    result = { lang: normalizeLang(top?.language), reliable: !!r?.isReliable && (top?.percentage ?? 0) >= 70 };
  } catch {
    /* algılama yoksa motor karar verir */
  }
  if (detectCache.size > 2000) detectCache.clear();
  detectCache.set(sample, result);
  return result;
}

/** Kaynak büyük harfle başlıyorsa çevirinin ilk harfini büyütür (motorlar başlıkları bazen küçük harfle döndürür). */
function capitalizeFirst(fragment: DocumentFragment, source: string, lang: string): void {
  const first = source.trimStart().match(/\p{L}/u)?.[0];
  if (!first || first === first.toLocaleLowerCase() || first !== first.toLocaleUpperCase()) return;
  const walker = document.createTreeWalker(fragment, NodeFilter.SHOW_TEXT);
  for (let node = walker.nextNode() as Text | null; node; node = walker.nextNode() as Text | null) {
    const i = node.data.search(/\p{L}/u);
    if (i < 0) continue;
    const ch = node.data[i];
    if (ch === ch.toLocaleLowerCase(lang) && ch !== ch.toLocaleUpperCase(lang)) {
      node.data = node.data.slice(0, i) + ch.toLocaleUpperCase(lang) + node.data.slice(i + 1);
    }
    return;
  }
}

/** Seçimin kapsadığı ilk ve son (boş olmayan) metin düğümleri. */
function rangeTextBounds(range: Range): { first: Text | null; last: Text | null } {
  const common = range.commonAncestorContainer;
  const root = common.nodeType === Node.TEXT_NODE ? (common.parentNode ?? common) : common;
  const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
  let first: Text | null = null;
  let last: Text | null = null;
  for (let n = walker.nextNode() as Text | null; n; n = walker.nextNode() as Text | null) {
    if (!/\S/.test(n.data) || !range.intersectsNode(n) || insideOurs(n)) continue;
    if (n === range.endContainer && range.endOffset === 0) continue;
    if (n === range.startContainer && range.startOffset >= n.data.length) continue;
    first ??= n;
    last = n;
  }
  return { first, last };
}

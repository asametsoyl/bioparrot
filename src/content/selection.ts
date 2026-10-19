import { markSvg } from '@/core/brand';
import { t } from '@/core/i18n';
import { languageName } from '@/core/languages';
import { getLocale } from '@/core/i18n';
import type { Settings } from '@/core/types';
import { sendToBackground } from '@/lib/messaging';
import { updateSettings } from '@/lib/settings';
import type { PageInfo } from './page-info';
import { ICONS, createUiHost, eventTarget, isEditable } from './ui-host';

const CSS = `
.bp-sel-icon { position: fixed; z-index: 2147483647; width: 22px; height: 22px; display: grid; place-items: center; cursor: pointer; animation: bp-pop .12s ease-out; }
.bp-sel-icon::before { content: ''; width: 12px; height: 12px; border-radius: 50%; background: radial-gradient(circle at 35% 35%, #4ade80, #16a34a 70%); border: 2px solid #fff; box-shadow: 0 1px 4px rgba(0,0,0,.35); transition: transform .12s ease; }
.bp-sel-icon:hover::before { transform: scale(1.35); }
@keyframes bp-pop { from { transform: scale(.85); opacity: 0; } }
.bp-card { position: fixed; z-index: 2147483647; width: min(440px, calc(100vw - 24px)); max-height: min(70vh, 560px); display: flex; flex-direction: column; background: var(--bp-bg); color: var(--bp-fg); border-radius: 14px; box-shadow: var(--bp-shadow); border: 1px solid var(--bp-line); overflow: hidden; animation: bp-pop .12s ease-out; }
.bp-card-head { display: flex; align-items: center; gap: 8px; padding: 9px 10px 9px 12px; border-bottom: 1px solid var(--bp-line); cursor: move; }
.bp-card-head .title { font-weight: 600; font-size: 13px; flex: 1; display: flex; align-items: center; gap: 6px; }
.bp-card-head .langs { font-size: 12px; color: var(--bp-muted); font-weight: 500; }
.bp-icon-btn { width: 28px; height: 28px; border-radius: 7px; display: grid; place-items: center; color: var(--bp-muted); }
.bp-icon-btn:hover { background: var(--bp-soft); color: var(--bp-fg); }
.bp-card-body { padding: 12px 14px; overflow: auto; }
.bp-source { color: var(--bp-muted); font-size: 13px; margin-bottom: 8px; max-height: 5.5em; overflow: hidden; display: -webkit-box; -webkit-line-clamp: 3; -webkit-box-orient: vertical; }
.bp-result { font-size: 15px; line-height: 1.6; white-space: pre-wrap; word-break: break-word; }
.bp-result.err { color: #dc2626; font-size: 13px; }
.bp-explain { margin-top: 12px; padding-top: 10px; border-top: 1px dashed var(--bp-line); font-size: 13.5px; line-height: 1.6; }
.bp-explain p { margin: 0 0 6px; }
.bp-explain ul { margin: 0 0 6px; padding-left: 18px; }
.bp-card-foot { display: flex; align-items: center; gap: 4px; padding: 6px 8px; border-top: 1px solid var(--bp-line); background: var(--bp-soft); }
.bp-chip { display: inline-flex; align-items: center; gap: 5px; height: 28px; padding: 0 10px; border-radius: 8px; font-size: 12.5px; color: var(--bp-fg); }
.bp-chip:hover { background: var(--bp-bg); color: var(--bp-green); }
.bp-chip[disabled] { opacity: .5; pointer-events: none; }
.bp-spacer { flex: 1; }
.bp-engine { font-size: 11.5px; color: var(--bp-muted); padding-right: 6px; }
`;

/** Seçili metin için çeviri simgesi ve çeviri kartı. */
export class SelectionPopup {
  private host: HTMLElement | null = null;
  private root: ShadowRoot | null = null;
  private icon: HTMLElement | null = null;
  private card: HTMLElement | null = null;
  private lastRect: DOMRect | null = null;
  private lastText = '';

  private lastRange: Range | null = null;

  constructor(
    private settings: () => Settings,
    private pageInfo: PageInfo,
    /** "Paragrafın altında" gösterim modu için çeviriyi sayfaya ekleyen geri çağırma. */
    private onInline?: (range: Range, text: string) => void,
  ) {}

  private wantsInline(): boolean {
    return !!this.onInline && this.settings().selection.display === 'inline';
  }

  attach(): void {
    document.addEventListener('mouseup', this.onMouseUp, true);
    document.addEventListener('mousedown', this.onMouseDown, true);
    document.addEventListener('keydown', this.onKeyDown, true);
    window.addEventListener('scroll', this.onScroll, { passive: true, capture: true });
  }

  detach(): void {
    document.removeEventListener('mouseup', this.onMouseUp, true);
    document.removeEventListener('mousedown', this.onMouseDown, true);
    document.removeEventListener('keydown', this.onKeyDown, true);
    window.removeEventListener('scroll', this.onScroll, true);
    this.host?.remove();
    this.host = this.root = this.icon = this.card = null;
  }

  private ensureHost(): ShadowRoot {
    if (!this.root) {
      const { host, root } = createUiHost('selection', CSS);
      this.host = host;
      this.root = root;
      const wrap = document.createElement('div');
      wrap.className = 'bp-root';
      root.appendChild(wrap);
    }
    return this.root;
  }

  private get container(): HTMLElement {
    return this.ensureHost().querySelector('.bp-root') as HTMLElement;
  }

  private isOwnEvent(e: Event): boolean {
    return !!this.host && e.composedPath().includes(this.host);
  }

  private onMouseDown = (e: MouseEvent) => {
    if (this.isOwnEvent(e)) return;
    this.hideIcon();
    if (this.card && !e.composedPath().includes(this.card)) this.closeCard();
  };

  private onMouseUp = (e: MouseEvent) => {
    if (this.isOwnEvent(e) || e.button !== 0) return;
    const s = this.settings();
    if (s.selection.mode === 'off') return;
    const target = eventTarget(e);
    if (!s.selection.inInputs && isEditable(target)) return;
    setTimeout(() => {
      const sel = window.getSelection();
      const text = sel?.toString().trim() ?? '';
      if (!sel || !text || text.length > 5000 || sel.rangeCount === 0) return;
      const range = sel.getRangeAt(0);
      const rect = range.getBoundingClientRect();
      this.lastRect = rect.width || rect.height ? rect : new DOMRect(e.clientX, e.clientY, 1, 1);
      this.lastText = text;
      this.lastRange = range.cloneRange();
      if (s.selection.mode === 'direct') {
        if (this.wantsInline()) this.onInline!(this.lastRange, text);
        else this.openCard(text, this.lastRect);
      } else this.showIcon(this.lastRect, e.clientX, e.clientY);
    }, 10);
  };

  private onKeyDown = (e: KeyboardEvent) => {
    if (e.key === 'Escape') {
      this.hideIcon();
      this.closeCard();
    }
  };

  private onScroll = () => {
    this.hideIcon();
  };

  private showIcon(rect: DOMRect, mx: number, my: number): void {
    this.hideIcon();
    const icon = document.createElement('div');
    icon.className = 'bp-sel-icon';
    icon.title = t('sel.translate');
    const x = Math.min(Math.max((mx || rect.right) + 4, 8), window.innerWidth - 30);
    const y = Math.min(Math.max((my || rect.bottom) + 6, 8), window.innerHeight - 30);
    icon.style.left = `${x}px`;
    icon.style.top = `${y}px`;
    icon.addEventListener('mousedown', (e) => e.preventDefault());
    icon.addEventListener('click', () => {
      this.hideIcon();
      if (this.wantsInline() && this.lastRange) {
        this.onInline!(this.lastRange, this.lastText);
        window.getSelection()?.removeAllRanges();
      } else this.openCard(this.lastText, this.lastRect ?? rect);
    });
    this.container.appendChild(icon);
    this.icon = icon;
  }

  private hideIcon(): void {
    this.icon?.remove();
    this.icon = null;
  }

  /** Mevcut seçimi (veya verilen metni) kart içinde çevirir. */
  openForSelection(text?: string): void {
    const sel = window.getSelection();
    const value = (text ?? sel?.toString() ?? '').trim();
    if (!value) return;
    let rect: DOMRect | null = null;
    if (sel && sel.rangeCount) {
      const range = sel.getRangeAt(0);
      if (this.wantsInline() && !range.collapsed && range.toString().trim()) {
        this.onInline!(range.cloneRange(), range.toString().trim());
        return;
      }
      rect = range.getBoundingClientRect();
    }
    this.openCard(value, rect && (rect.width || rect.height) ? rect : new DOMRect(window.innerWidth / 2 - 220, 80, 1, 1));
  }

  private closeCard(): void {
    this.card?.remove();
    this.card = null;
    window.speechSynthesis?.cancel();
  }

  private openCard(text: string, rect: DOMRect): void {
    this.closeCard();
    const s = this.settings();
    const ui = getLocale();
    const card = document.createElement('div');
    card.className = 'bp-card';
    card.innerHTML = `
      <div class="bp-card-head">
        <div class="title"><span class="mark"></span><span>BioParrot</span><span class="langs"></span></div>
        <button class="bp-icon-btn" data-act="settings">${ICONS.settings}</button>
        <button class="bp-icon-btn" data-act="close">${ICONS.close}</button>
      </div>
      <div class="bp-card-body">
        <div class="bp-source"></div>
        <div class="bp-result"><span class="bp-spin"></span></div>
        <div class="bp-explain" hidden></div>
      </div>
      <div class="bp-card-foot">
        <button class="bp-chip" data-act="copy">${ICONS.copy}<span></span></button>
        <button class="bp-chip" data-act="speak">${ICONS.speak}<span></span></button>
        <button class="bp-chip" data-act="explain">${ICONS.spark}<span></span></button>
        <button class="bp-chip" data-act="glossary" hidden>${ICONS.book}<span></span></button>
        <span class="bp-spacer"></span><span class="bp-engine"></span>
      </div>`;
    card.querySelector('.mark')!.innerHTML = markSvg(18);
    card.querySelector('.langs')!.textContent = `→ ${languageName(s.targetLanguage, ui)}`;
    card.querySelector('.bp-source')!.textContent = text;
    const label = (act: string, key: Parameters<typeof t>[0]) => {
      const btn = card.querySelector<HTMLElement>(`[data-act="${act}"]`)!;
      btn.querySelector('span')!.textContent = t(key);
      btn.title = t(key);
    };
    label('copy', 'common.copy');
    label('speak', 'common.speak');
    // Kopyala/seslendir yalnızca simge olarak gösterilir (alt çubuk sığsın diye).
    card.querySelectorAll<HTMLElement>('[data-act="copy"] span, [data-act="speak"] span').forEach((el) => (el.hidden = true));
    label('explain', 'sel.explain');
    label('glossary', 'sel.addGlossary');
    card.querySelector<HTMLElement>('[data-act="settings"]')!.title = t('common.settings');
    card.querySelector<HTMLElement>('[data-act="close"]')!.title = t('common.close');
    card.querySelector('.bp-engine')!.textContent = s.engines[s.selection.engineId || s.engineId]?.name ?? '';

    // Konum: seçimin altına; altta yer yoksa ve üstte daha çok yer varsa üstüne.
    const width = Math.min(440, window.innerWidth - 24);
    const left = Math.min(Math.max(rect.left, 12), window.innerWidth - width - 12);
    const spaceBelow = window.innerHeight - rect.bottom - 22;
    const spaceAbove = rect.top - 22;
    card.style.left = `${left}px`;
    if (spaceBelow >= 280 || spaceBelow >= spaceAbove) {
      card.style.top = `${Math.max(rect.bottom + 10, 12)}px`;
      card.style.maxHeight = `${Math.max(spaceBelow, 180)}px`;
    } else {
      card.style.bottom = `${window.innerHeight - rect.top + 10}px`;
      card.style.maxHeight = `${Math.max(spaceAbove, 180)}px`;
    }

    this.container.appendChild(card);
    this.card = card;
    this.makeDraggable(card);

    const result = card.querySelector<HTMLElement>('.bp-result')!;
    let translation = '';
    const glossaryBtn = card.querySelector<HTMLElement>('[data-act="glossary"]')!;
    const short = text.split(/\s+/).length <= 6 && text.length <= 80;

    card.addEventListener('click', async (e) => {
      const btn = (e.target as HTMLElement).closest<HTMLElement>('[data-act]');
      if (!btn) return;
      const act = btn.dataset.act;
      if (act === 'close') this.closeCard();
      else if (act === 'settings') void sendToBackground({ type: 'openOptions' }).catch(() => {});
      else if (act === 'copy') {
        await navigator.clipboard.writeText(translation || text).catch(() => {});
        btn.title = t('common.copied');
        btn.style.color = 'var(--bp-green)';
        setTimeout(() => {
          btn.title = t('common.copy');
          btn.style.color = '';
        }, 1200);
      } else if (act === 'speak') {
        const synth = window.speechSynthesis;
        if (!synth) return;
        synth.cancel();
        const u = new SpeechSynthesisUtterance(translation || text);
        u.lang = translation ? s.targetLanguage : '';
        synth.speak(u);
      } else if (act === 'explain') {
        btn.setAttribute('disabled', '');
        void this.explain(card, text);
      } else if (act === 'glossary' && translation) {
        await updateSettings((d) => {
          d.glossary.push({ id: crypto.randomUUID(), source: text, target: translation.trim() });
        });
        btn.querySelector('span')!.textContent = t('sel.addedGlossary');
        btn.setAttribute('disabled', '');
      }
    });

    sendToBackground({
      type: 'translate',
      req: {
        texts: [text],
        from: s.sourceLanguage,
        to: s.targetLanguage,
        engineId: s.selection.engineId || s.engineId,
        rich: false,
        context: s.page.contextAware ? { title: this.pageInfo.title, summary: this.pageInfo.summary } : undefined,
      },
    })
      .then((res) => {
        translation = res.translations[0] ?? '';
        result.textContent = translation;
        if (short && translation) glossaryBtn.hidden = false;
      })
      .catch((err) => {
        result.classList.add('err');
        result.textContent = t('err.generic', { msg: String(err?.message ?? err) });
      });
  }

  private async explain(card: HTMLElement, text: string): Promise<void> {
    const s = this.settings();
    const box = card.querySelector<HTMLElement>('.bp-explain')!;
    box.hidden = false;
    box.innerHTML = '<span class="bp-spin"></span>';
    try {
      const out = await sendToBackground({
        type: 'explain',
        text,
        to: s.targetLanguage,
        context: { title: this.pageInfo.title },
      });
      box.replaceChildren(renderMiniMarkdown(out));
    } catch (err) {
      box.textContent = String((err as Error)?.message ?? err);
    }
  }

  private makeDraggable(card: HTMLElement): void {
    const head = card.querySelector<HTMLElement>('.bp-card-head')!;
    let sx = 0;
    let sy = 0;
    let ox = 0;
    let oy = 0;
    const move = (e: PointerEvent) => {
      card.style.left = `${ox + e.clientX - sx}px`;
      card.style.top = `${oy + e.clientY - sy}px`;
      card.style.bottom = 'auto';
    };
    head.addEventListener('pointerdown', (e) => {
      if ((e.target as HTMLElement).closest('button')) return;
      const r = card.getBoundingClientRect();
      sx = e.clientX;
      sy = e.clientY;
      ox = r.left;
      oy = r.top;
      head.setPointerCapture(e.pointerId);
      head.addEventListener('pointermove', move);
      head.addEventListener(
        'pointerup',
        () => {
          head.removeEventListener('pointermove', move);
        },
        { once: true },
      );
    });
  }
}

/** Model çıktısını güvenli biçimde (innerHTML kullanmadan) basit paragraf/liste/kalın olarak gösterir. */
export function renderMiniMarkdown(text: string): DocumentFragment {
  const frag = document.createDocumentFragment();
  const inline = (parent: HTMLElement, line: string) => {
    const parts = line.split(/(\*\*[^*]+\*\*)/g);
    for (const part of parts) {
      if (/^\*\*[^*]+\*\*$/.test(part)) {
        const b = document.createElement('strong');
        b.textContent = part.slice(2, -2);
        parent.appendChild(b);
      } else if (part) {
        parent.appendChild(document.createTextNode(part.replace(/^#+\s*/, '')));
      }
    }
  };
  let list: HTMLUListElement | null = null;
  for (const raw of text.split('\n')) {
    const line = raw.trim();
    if (!line) {
      list = null;
      continue;
    }
    const bullet = /^([-*•]|\d+[.)])\s+(.*)$/.exec(line);
    if (bullet) {
      if (!list) {
        list = document.createElement('ul');
        frag.appendChild(list);
      }
      const li = document.createElement('li');
      inline(li, bullet[2]);
      list.appendChild(li);
    } else {
      list = null;
      const p = document.createElement('p');
      inline(p, line);
      frag.appendChild(p);
    }
  }
  return frag;
}

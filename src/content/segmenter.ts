import { CLS } from '@/core/constants';
import { escapeHtml } from '@/core/markup';
import { REFERENCE_HEADING } from './site-rules';

/** Hiçbir zaman çevrilmeyen ve içine girilmeyen etiketler. */
const SKIP_TAGS = new Set([
  'SCRIPT', 'STYLE', 'NOSCRIPT', 'TEMPLATE', 'HEAD', 'META', 'LINK', 'TITLE', 'BASE',
  'TEXTAREA', 'INPUT', 'SELECT', 'OPTION', 'OPTGROUP', 'DATALIST', 'BUTTON',
  'CANVAS', 'VIDEO', 'AUDIO', 'IFRAME', 'OBJECT', 'EMBED', 'PICTURE', 'SOURCE', 'TRACK', 'MAP', 'AREA',
  'SVG', 'MATH', 'HR', 'DIALOG',
]);

/** Paragraf içinde olduğu gibi korunan (çevrilmeyen) satır içi öğeler. */
const ATOMIC_TAGS = new Set([
  'CODE', 'KBD', 'SAMP', 'VAR', 'TT', 'IMG', 'SVG', 'MATH', 'MJX-CONTAINER', 'SUP', 'SUB',
  'INPUT', 'SELECT', 'TEXTAREA', 'BUTTON', 'CANVAS', 'VIDEO', 'AUDIO', 'PICTURE', 'OBJECT', 'EMBED', 'IFRAME', 'WBR',
]);
const ATOMIC_SELECTOR = '.katex, .MathJax, .MathJax_SVG, .MathJax_CHTML, mjx-container, .mwe-math-element, .ltx_Math, .ltx_math_unparsed, [translate="no"], .notranslate';

/** Satır içi ama kaynak/çeviride biçim etiketi olarak korunanlar. */
const PAIRED: Record<string, 'a' | 'b' | 'i'> = { A: 'a', B: 'b', STRONG: 'b', I: 'i', EM: 'i', CITE: 'i' };

/** Bir öğenin içinde blok öğe olup olmadığını hızlıca anlamak için. */
const BLOCK_DESCENDANT_SELECTOR =
  'div,p,ul,ol,li,dl,dt,dd,table,tr,td,th,section,article,header,footer,nav,aside,main,h1,h2,h3,h4,h5,h6,blockquote,pre,figure,figcaption,form,fieldset,address,details,summary,hr';

export type NodeClass = 'ours' | 'skip' | 'hidden' | 'atomic' | 'inline' | 'block';

export interface RawUnit {
  /** Gözlemlenen ve stil uygulanan kapsayıcı öğe. */
  root: Element;
  /** Düğümlerin gerçek ebeveyni (gölge kökte ShadowRoot olabilir). */
  parent: Element | ShadowRoot;
  nodes: Node[];
  /** Koşu kapsayıcının tüm anlamlı içeriğini kapsıyor. */
  whole: boolean;
}

export interface Serialized {
  /** Motorlara gönderilen zengin metin (HTML kaçışlı, etiketli). */
  markup: string;
  /** Okunabilir düz metin (kod/formül metni dahil). */
  plain: string;
  /** Kod/formül hariç düz metin (dil algılama ve filtreler için). */
  prose: string;
  /** Etiket indeksi → kaynak öğe. */
  elements: Element[];
}

export interface SegmenterOptions {
  exclude: string | null;
  skipCode: boolean;
  skipReferences: boolean;
}

export function isOurs(node: Node): boolean {
  if (node.nodeType !== Node.ELEMENT_NODE) return false;
  const el = node as Element;
  return el.classList.contains(CLS.wrap) || el.hasAttribute('data-bp-ui');
}

export function insideOurs(node: Node): boolean {
  const el = node.nodeType === Node.ELEMENT_NODE ? (node as Element) : node.parentElement;
  return !!el?.closest(`.${CLS.wrap}, [data-bp-ui]`);
}

function isTrivial(n: Node): boolean {
  if (n.nodeType === Node.TEXT_NODE) return !/\S/.test((n as Text).data);
  return n.nodeType === Node.ELEMENT_NODE && (n as Element).tagName === 'BR';
}

export class Segmenter {
  private cache = new WeakMap<Element, NodeClass>();

  constructor(public opts: SegmenterOptions) {}

  reset(opts?: SegmenterOptions): void {
    if (opts) this.opts = opts;
    this.cache = new WeakMap();
  }

  forget(el: Element): void {
    this.cache.delete(el);
  }

  classify(el: Element): NodeClass {
    const cached = this.cache.get(el);
    if (cached) return cached;
    const c = this.computeClass(el);
    this.cache.set(el, c);
    return c;
  }

  private computeClass(el: Element): NodeClass {
    if (isOurs(el)) return 'ours';
    const tag = el.tagName.toUpperCase();
    if (tag === 'BR') return 'atomic';
    let display = '';
    const getDisplay = () => (display ||= getComputedStyle(el).display);

    const excluded = !!this.opts.exclude && safeMatches(el, this.opts.exclude);
    const atomicLike = ATOMIC_TAGS.has(tag) || safeMatches(el, ATOMIC_SELECTOR);
    if (excluded || atomicLike || (this.opts.skipCode && tag === 'PRE')) {
      if (getDisplay() === 'none') return 'hidden';
      // Satır içinde duranlar paragrafın parçası olarak korunur, bloklar tamamen atlanır.
      return isInlineDisplay(getDisplay()) && tag !== 'PRE' ? 'atomic' : 'skip';
    }
    if (SKIP_TAGS.has(tag)) return 'skip';
    if (el.getAttribute('contenteditable') === 'true' || el.getAttribute('contenteditable') === '') return 'skip';

    const d = getDisplay();
    if (d === 'none') return 'hidden';
    if (!isInlineDisplay(d)) return 'block';
    // İçinde blok öğe barındıran satır içi öğeler (örn. kart şeklindeki bağlantılar) kapsayıcı gibi taranır.
    if (el.firstElementChild && el.querySelector(BLOCK_DESCENDANT_SELECTOR)) {
      const inner = el.querySelectorAll(BLOCK_DESCENDANT_SELECTOR);
      for (const child of Array.from(inner).slice(0, 20)) {
        const cd = getComputedStyle(child).display;
        if (!isInlineDisplay(cd) && cd !== 'none') return 'block';
      }
    }
    return 'inline';
  }

  /** Bir kapsayıcının doğrudan çocuklarındaki satır içi koşuları ve taranacak blok çocukları bulur. */
  runsOf(container: Element | ShadowRoot): { units: RawUnit[]; blocks: (Element | ShadowRoot)[] } {
    const root: Element = container instanceof ShadowRoot ? container.host : container;
    const units: RawUnit[] = [];
    const blocks: (Element | ShadowRoot)[] = [];
    let run: Node[] = [];
    let sawBlock = false;
    let refMode = false;
    let refLevel = 7;

    const flush = () => {
      let s = 0;
      let e = run.length;
      while (s < e && isTrivial(run[s])) s++;
      while (e > s && isTrivial(run[e - 1])) e--;
      if (e > s) units.push({ root, parent: container, nodes: run.slice(s, e), whole: false });
      run = [];
    };

    for (const child of Array.from(container.childNodes)) {
      if (child.nodeType === Node.TEXT_NODE) {
        if (refMode) continue;
        if (run.length || /\S/.test((child as Text).data)) run.push(child);
        continue;
      }
      if (child.nodeType !== Node.ELEMENT_NODE) continue;
      const el = child as Element;
      const c = this.classify(el);
      if (c === 'ours') continue;
      if (c === 'skip') {
        flush();
        sawBlock = true;
        continue;
      }
      if (c === 'hidden') {
        // Görünür hale gelince çevrilsin diye içi yine de taranır (kesişim gözlemcisi tetikler).
        if (!refMode) blocks.push(el);
        if (el.shadowRoot) blocks.push(el.shadowRoot);
        continue;
      }
      if (c === 'atomic' || c === 'inline') {
        if (refMode) continue;
        run.push(el);
        if (c === 'inline' && el.shadowRoot) blocks.push(el.shadowRoot);
        continue;
      }
      // Blok çocuk
      flush();
      sawBlock = true;
      if (this.opts.skipReferences) {
        const level = headingLevel(el);
        if (level) {
          const text = (el.textContent ?? '').trim();
          if (REFERENCE_HEADING.test(text)) {
            refMode = true;
            refLevel = level;
            blocks.push(el);
            continue;
          }
          if (refMode && level <= refLevel) refMode = false;
        }
        if (refMode) continue;
      }
      blocks.push(el);
      if (el.shadowRoot) blocks.push(el.shadowRoot);
    }
    flush();
    if (units.length === 1 && !sawBlock && container === root) units[0].whole = true;
    return { units, blocks };
  }

  /** Bir alt ağacı parça parça tarar; ana iş parçacığını kilitlememek için arada nefes alır. */
  async scan(start: Element | ShadowRoot, onUnits: (units: RawUnit[]) => void, isCancelled: () => boolean): Promise<void> {
    const stack: (Element | ShadowRoot)[] = [start];
    let deadline = performance.now() + 12;
    while (stack.length) {
      if (isCancelled()) return;
      const c = stack.pop()!;
      if (c instanceof Element && this.opts.skipReferences && isReferenceContainer(c)) continue;
      const { units, blocks } = this.runsOf(c);
      if (units.length) onUnits(units);
      for (let i = blocks.length - 1; i >= 0; i--) stack.push(blocks[i]);
      if (performance.now() > deadline) {
        await yieldToMain();
        deadline = performance.now() + 12;
      }
    }
  }

  /** Bir öğeyi içeren en yakın blok kapsayıcı. */
  containerOf(node: Node): Element | ShadowRoot | null {
    let el: Element | null = node.nodeType === Node.ELEMENT_NODE ? (node as Element) : node.parentElement;
    if (el && node.nodeType === Node.ELEMENT_NODE) {
      const c = this.classify(el);
      if (c === 'block' || c === 'hidden') return el;
    }
    while (el) {
      const parent: Node | null = el.parentNode;
      if (!parent) return null;
      if (parent instanceof ShadowRoot) return parent;
      el = parent as Element;
      if (el === document.body || el === document.documentElement) return el;
      const c = this.classify(el);
      if (c === 'block' || c === 'hidden') return el;
      if (c === 'skip' || c === 'ours') return null;
    }
    return null;
  }

  /** Koşudaki düğümleri motorlara gönderilecek zengin metne çevirir. */
  serialize(nodes: Node[]): Serialized {
    const elements: Element[] = [];
    let markup = '';
    let plain = '';
    let prose = '';

    const walk = (node: Node) => {
      if (node.nodeType === Node.TEXT_NODE) {
        const text = (node as Text).data.replace(/[\s ]+/g, ' ');
        markup += escapeHtml(text);
        plain += text;
        prose += text;
        return;
      }
      if (node.nodeType !== Node.ELEMENT_NODE) return;
      const el = node as Element;
      const c = this.classify(el);
      if (c === 'ours' || c === 'hidden' || c === 'skip') return;
      const tag = el.tagName.toUpperCase();
      if (tag === 'BR') {
        markup += '<br>';
        plain += '\n';
        prose += '\n';
        return;
      }
      if (c === 'atomic') {
        const index = elements.push(el) - 1;
        const raw = (el.textContent || el.getAttribute('alt') || el.getAttribute('aria-label') || '').replace(/\s+/g, ' ').trim();
        if (!raw || isCitationLike(el, raw)) {
          // Atıf numaraları ve görseller için içeriksiz yer tutucu: motorlar onu doğru yere koyar ve metni çoğaltmaz.
          markup += `<img i="${index}">`;
        } else {
          const text = raw.length > 80 ? raw.slice(0, 80) + '…' : raw;
          markup += `<code class="notranslate" translate="no" i="${index}">${escapeHtml(text)}</code>`;
        }
        plain += raw;
        prose += ' ';
        return;
      }
      const paired = PAIRED[tag];
      if (paired && (paired !== 'a' || el.hasAttribute('href'))) {
        const index = elements.push(el) - 1;
        markup += `<${paired} i="${index}">`;
        el.childNodes.forEach(walk);
        markup += `</${paired}>`;
        return;
      }
      el.childNodes.forEach(walk);
    };
    nodes.forEach(walk);

    markup = markup.replace(/^(?:\s|<br>)+|(?:\s|<br>)+$/g, '').replace(/ {2,}/g, ' ');
    plain = plain.replace(/[ \t]{2,}/g, ' ').trim();
    prose = prose.replace(/\s{2,}/g, ' ').trim();
    return { markup, plain, prose, elements };
  }
}

/** Atıf/dipnot işareti gibi duran kısa öğeler (¹, [12], 3–5, †). */
function isCitationLike(el: Element, text: string): boolean {
  const tag = el.tagName.toUpperCase();
  if (tag !== 'SUP' && tag !== 'A' && !el.classList.contains('reference')) return false;
  return /^[\s\[(]*[\d\s,;\-–—]+[\])\s]*$/.test(text) || /^[*†‡§¶#]+$/.test(text);
}

/** Satır içi akışta kalan display değerleri (MathML'in "math" değeri dahil). */
function isInlineDisplay(d: string): boolean {
  return d.startsWith('inline') || d === 'math' || d === 'contents' || d.startsWith('ruby');
}

function safeMatches(el: Element, selector: string): boolean {
  try {
    return el.matches(selector);
  } catch {
    return false;
  }
}

function headingLevel(el: Element): number {
  const m = /^H([1-6])$/.exec(el.tagName);
  if (m) return Number(m[1]);
  // <div class="mw-heading"><h2>…</h2></div> gibi sarmalanmış başlıklar
  if (el.childElementCount <= 3) {
    const h = el.querySelector(':scope > h1, :scope > h2, :scope > h3, :scope > h4, :scope > h5, :scope > h6');
    if (h) return Number(h.tagName[1]);
  }
  return 0;
}

const REFERENCE_CONTAINER_SELECTOR = [
  '#references', '#bibliography', '#reference-list', '#ref-list', '.references', '.bibliography', '.ref-list', '.reflist',
  '.ltx_bibliography', '.c-article-references', '.article-section__references', '.citation-list', '.html-references',
  '#Bib1-content', 'ol.references', 'ul.references', '[data-title="References"]', '.mw-references-wrap', '.csl-bib-body',
].join(',');

function isReferenceContainer(el: Element): boolean {
  if (safeMatches(el, REFERENCE_CONTAINER_SELECTOR)) return true;
  if (el.tagName === 'SECTION' || el.tagName === 'DIV') {
    const first = el.firstElementChild;
    if (first && headingLevel(first) && REFERENCE_HEADING.test((first.textContent ?? '').trim())) return true;
  }
  return false;
}

function yieldToMain(): Promise<void> {
  const sched = (globalThis as any).scheduler;
  if (sched?.yield) return sched.yield();
  return new Promise((resolve) => setTimeout(resolve, 0));
}

/** Çevrilmeye değer metin mi (harf içerir, adres/e-posta/kod tanımlayıcısı değil). */
export function shouldTranslate(prose: string, plain: string): boolean {
  const t = prose.trim();
  if (!t) return false;
  if (!/\p{L}/u.test(t)) return false;
  const cjk = /[぀-ヿ㐀-鿿가-힯]/.test(t);
  if (t.length < 2 && !cjk) return false;
  if (/^(https?:\/\/|www\.)\S+$/i.test(t)) return false;
  if (/^[\w.+-]+@[\w-]+\.[\w.-]+$/.test(t)) return false;
  // Tek kelimelik kod tanımlayıcıları (camelCase, snake_case, fonksiyon())
  if (!/\s/.test(t) && (/^[a-z]+[A-Z][\w]*$/.test(t) || /^\w+_\w+$/.test(t) || /^[\w.]+\(\)$/.test(t))) return false;
  // Yalnızca kısaltma/sembol gibi görünen çok kısa büyük harfli metinler (DNA, PCR) çevrilmez
  if (!/\s/.test(plain.trim()) && /^[A-Z0-9][A-Z0-9\-/+.]{0,7}$/.test(t)) return false;
  return true;
}

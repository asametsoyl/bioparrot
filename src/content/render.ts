import { CLS } from '@/core/constants';
import { tokenizeMarkup } from '@/core/markup';

/** Kopyalanan öğelerden kimlik ve olay özniteliklerini temizler. */
function sanitizeClone<T extends Element>(el: T): T {
  const all = [el, ...Array.from(el.querySelectorAll('*'))];
  for (const node of all) {
    node.removeAttribute('id');
    for (const attr of Array.from(node.attributes)) {
      if (attr.name.startsWith('on')) node.removeAttribute(attr.name);
    }
  }
  return el;
}

interface Frame {
  node: Node;
  tag: string;
  /** İçeriği yok sayılır (atomik öğenin modeldeki metni). */
  sink: boolean;
}

/**
 * Zengin çeviri metninden DOM parçası üretir. Etiketler kaynak öğelerin kopyalarıyla
 * değiştirilir (bağlantılar tıklanabilir kalır), metin her zaman textNode olarak eklenir.
 */
export function buildFragment(markup: string, elements: Element[], doc: Document = document): DocumentFragment {
  const fragment = doc.createDocumentFragment();
  const stack: Frame[] = [{ node: fragment, tag: '#root', sink: false }];
  const top = () => stack[stack.length - 1];
  const sinking = () => stack.some((f) => f.sink);

  for (const token of tokenizeMarkup(markup)) {
    if (token.type === 'text') {
      if (!sinking() && token.text) top().node.appendChild(doc.createTextNode(token.text));
      continue;
    }
    if (token.type === 'br') {
      if (!sinking()) top().node.appendChild(doc.createElement('br'));
      continue;
    }
    if (token.type === 'open') {
      const index = token.index !== null && /^\d+$/.test(token.index) ? Number(token.index) : -1;
      const source = index >= 0 ? elements[index] : undefined;
      if (sinking()) {
        if (!token.selfClosing) stack.push({ node: top().node, tag: token.tag, sink: true });
        continue;
      }
      if (source && (token.tag === 'code' || token.tag === 'img')) {
        top().node.appendChild(sanitizeClone(source.cloneNode(true) as Element));
        if (!token.selfClosing) stack.push({ node: top().node, tag: token.tag, sink: true });
        continue;
      }
      if (token.tag === 'img') continue;
      if (source && (token.tag === 'a' || token.tag === 'b' || token.tag === 'i' || token.tag === 'span')) {
        const clone = sanitizeClone(source.cloneNode(false) as Element);
        top().node.appendChild(clone);
        if (!token.selfClosing) stack.push({ node: clone, tag: token.tag, sink: false });
        continue;
      }
      // Bilinmeyen etiket (örn. sözlük sarmalayıcısı): içeriği mevcut düğüme aktarılır.
      if (!token.selfClosing) stack.push({ node: top().node, tag: token.tag, sink: false });
      continue;
    }
    // close
    for (let i = stack.length - 1; i > 0; i--) {
      if (stack[i].tag === token.tag) {
        stack.length = i;
        break;
      }
    }
  }
  return fragment;
}

/** Düz metin çeviriyi (satır sonlarıyla) DOM parçasına çevirir. */
export function plainFragment(text: string, doc: Document = document): DocumentFragment {
  const fragment = doc.createDocumentFragment();
  text.split('\n').forEach((line, i) => {
    if (i > 0) fragment.appendChild(doc.createElement('br'));
    if (line) fragment.appendChild(doc.createTextNode(line));
  });
  return fragment;
}

export function createWrapper(doc: Document, opts: { inline: boolean; styleId: string; unitId: number }): HTMLElement {
  const wrap = doc.createElement('font');
  wrap.className = wrapperClass(opts.inline, opts.styleId);
  wrap.setAttribute('data-bp-unit', String(opts.unitId));
  wrap.setAttribute('translate', 'no');
  return wrap;
}

export function wrapperClass(inline: boolean, styleId: string): string {
  return `${CLS.wrap} ${inline ? CLS.inline : CLS.block} notranslate bp-style-${styleId}`;
}

export function createLoading(doc: Document, unitId: number): HTMLElement {
  const wrap = doc.createElement('font');
  wrap.className = `${CLS.wrap} ${CLS.inline} notranslate`;
  wrap.setAttribute('data-bp-unit', String(unitId));
  const spin = doc.createElement('font');
  spin.className = CLS.loading;
  wrap.appendChild(spin);
  return wrap;
}

export function createError(doc: Document, unitId: number, message: string, label: string): HTMLElement {
  const wrap = doc.createElement('font');
  wrap.className = `${CLS.wrap} ${CLS.inline} notranslate`;
  wrap.setAttribute('data-bp-unit', String(unitId));
  const err = doc.createElement('font');
  err.className = CLS.error;
  err.textContent = `⚠ ${label}`;
  err.title = message;
  err.setAttribute('role', 'button');
  err.tabIndex = 0;
  wrap.appendChild(err);
  return wrap;
}

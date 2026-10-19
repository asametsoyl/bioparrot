import type { GlossaryEntry } from './types';

/**
 * Zengin metin biçimi (motorlar arasında ortak):
 *   <a i="0">…</a>, <b i="1">…</b>, <i i="2">…</i>, <span i="3">…</span>  → eşli satır içi öğeler
 *   <code class="notranslate" translate="no" i="4">…</code>              → hiç çevrilmeyen atomik öğeler
 *   <br>                                                                  → satır sonu
 * Metin HTML kaçışlıdır (&amp; &lt; &gt;). Google, Microsoft, DeepL (HTML modu) ve
 * yapay zeka modelleri bu biçimi etiketleri koruyarak çevirir.
 */

const ENTITY_MAP: Record<string, string> = {
  amp: '&',
  lt: '<',
  gt: '>',
  quot: '"',
  apos: "'",
  nbsp: ' ',
  ensp: ' ',
  emsp: ' ',
  thinsp: ' ',
  ndash: '–',
  mdash: '—',
  hellip: '…',
  laquo: '«',
  raquo: '»',
  lsquo: '‘',
  rsquo: '’',
  ldquo: '“',
  rdquo: '”',
};

export function escapeHtml(s: string): string {
  return s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
}

export function escapeAttr(s: string): string {
  return escapeHtml(s).replace(/"/g, '&quot;');
}

export function unescapeHtml(s: string): string {
  if (!s.includes('&')) return s;
  return s.replace(/&(#x[0-9a-f]+|#\d+|[a-z]+);/gi, (m, code: string) => {
    if (code[0] === '#') {
      const n = code[1] === 'x' || code[1] === 'X' ? parseInt(code.slice(2), 16) : parseInt(code.slice(1), 10);
      return Number.isFinite(n) && n > 0 && n < 0x110000 ? String.fromCodePoint(n) : m;
    }
    return ENTITY_MAP[code.toLowerCase()] ?? m;
  });
}

/** Düz metni zengin biçime çevirir (kaçış + satır sonları). */
export function plainToMarkup(s: string): string {
  return escapeHtml(s).replace(/\r?\n/g, '<br>');
}

export type MarkupToken =
  | { type: 'text'; text: string }
  | { type: 'open'; tag: string; index: string | null; selfClosing: boolean }
  | { type: 'close'; tag: string }
  | { type: 'br' };

const TAG_RE = /<\s*(\/)?\s*([a-zA-Z][a-zA-Z0-9-]*)((?:\s+[^>]*?)?)\s*(\/)?\s*>/g;
const VOID_TAGS = new Set(['img', 'hr', 'wbr', 'input']);
const INDEX_RE = /\bi\s*=\s*["']?([\w-]+)["']?/;

export function tokenizeMarkup(s: string): MarkupToken[] {
  const tokens: MarkupToken[] = [];
  let last = 0;
  TAG_RE.lastIndex = 0;
  let m: RegExpExecArray | null;
  while ((m = TAG_RE.exec(s))) {
    if (m.index > last) tokens.push({ type: 'text', text: unescapeHtml(s.slice(last, m.index)) });
    last = TAG_RE.lastIndex;
    const tag = m[2].toLowerCase();
    if (tag === 'br') {
      tokens.push({ type: 'br' });
      continue;
    }
    if (m[1]) {
      tokens.push({ type: 'close', tag });
      continue;
    }
    const idx = INDEX_RE.exec(m[3] ?? '');
    tokens.push({ type: 'open', tag, index: idx ? idx[1] : null, selfClosing: !!m[4] || VOID_TAGS.has(tag) });
  }
  if (last < s.length) tokens.push({ type: 'text', text: unescapeHtml(s.slice(last)) });
  return tokens;
}

/** Etiketleri kaldırıp düz metne çevirir (<br> → satır sonu). */
export function markupToPlain(s: string): string {
  if (!s.includes('<') && !s.includes('&')) return s;
  return tokenizeMarkup(s)
    .map((t) => (t.type === 'text' ? t.text : t.type === 'br' ? '\n' : ''))
    .join('');
}

export function hasMarkupTags(s: string): boolean {
  return /<\s*\/?\s*[a-zA-Z]/.test(s);
}

/** Metnin etiket dışındaki bölümlerine fonksiyon uygular. */
function mapTextParts(s: string, fn: (text: string) => string): string {
  return s
    .split(/(<[^>]*>)/)
    .map((part) => (part.startsWith('<') && part.endsWith('>') ? part : fn(part)))
    .join('');
}

function escapeRegExp(s: string): string {
  return s.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
}

function termRegex(entry: GlossaryEntry): RegExp {
  return new RegExp(`(?<![\\p{L}\\p{N}_])${escapeRegExp(entry.source.trim())}(?![\\p{L}\\p{N}_])`, entry.caseSensitive ? 'gu' : 'giu');
}

/** Metinlerde geçen sözlük girdilerini döner. */
export function relevantGlossary(texts: string[], glossary: GlossaryEntry[] | undefined, max = 60): GlossaryEntry[] {
  if (!glossary?.length) return [];
  const joined = markupToPlain(texts.join('\n'));
  const out: GlossaryEntry[] = [];
  for (const entry of glossary) {
    if (!entry.source?.trim() || !entry.target?.trim()) continue;
    if (termRegex(entry).test(joined)) out.push(entry);
    if (out.length >= max) break;
  }
  return out;
}

/**
 * Klasik çeviri motorları için: sözlük terimlerini hedef karşılığıyla değiştirip
 * "çevirme" etiketiyle sarar; motor bu kısmı olduğu gibi bırakır.
 */
export function applyGlossaryMarkup(s: string, glossary: GlossaryEntry[]): string {
  if (!glossary.length) return s;
  let n = 0;
  return mapTextParts(s, (text) => {
    let out = text;
    for (const entry of glossary) {
      const re = termRegex({ ...entry, source: escapeHtml(entry.source.trim()) });
      out = out.replace(re, () => `<span class="notranslate" translate="no" i="g${n++}">${escapeHtml(entry.target.trim())}</span>`);
    }
    return out;
  });
}

/** applyGlossaryMarkup'ın eklediği sarmalayıcıları kaldırır (içeriği bırakır). */
export function stripGlossaryMarkup(s: string): string {
  if (!s.includes('i="g')) return s;
  return s.replace(/<span\b[^>]*\bi="g\d+"[^>]*>([\s\S]*?)<\/span>/g, '$1');
}

/**
 * Klasik motorların HTML modunda etiketlerin çevresine eklediği gereksiz boşlukları temizler:
 * "( <i>p</i> < 0,05 )" → "(<i>p</i> < 0,05)", "ölçüldü<img i="0"> ." → "ölçüldü<img i="0">.".
 */
export function tidyMarkup(s: string): string {
  return s
    .replace(/([(\[{“‘«])\s+(?=<[a-z])/g, '$1')
    .replace(/(<\/[a-z]+>|<img\b[^>]*>)\s+(?=[)\]}”’».,;:!?])/g, '$1')
    .replace(/ {2,}/g, ' ');
}

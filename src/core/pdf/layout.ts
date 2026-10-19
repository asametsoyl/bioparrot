/**
 * PDF sayfa düzeni analizi: pdf.js metin öğelerinden satırları ve paragraf bloklarını çıkarır,
 * iki sütunlu makaleleri, formülleri, başlık/altbilgiyi, şekil içi metinleri ve kaynakçayı ayırt eder.
 * pdf.js'e bağımlı değildir; giriş olarak ham metin öğelerini ve görünüm dönüşümünü alır.
 */
import { escapeHtml } from '../markup';

export interface RawTextItem {
  str: string;
  /** PDF metin matrisi [a, b, c, d, e, f]. */
  transform: number[];
  width: number;
  height: number;
  fontName: string;
  hasEOL?: boolean;
}

export interface RawPage {
  index: number;
  /** Ölçek 1'deki görünüm genişliği/yüksekliği (PDF noktası). */
  width: number;
  height: number;
  /** pdf.js viewport.transform (ölçek 1). */
  viewportTransform: number[];
  items: RawTextItem[];
  /** fontName → gerçek yazı tipi adı (örn. "ABCDEF+CMMI10"). */
  fonts: Record<string, string>;
}

export interface Glyph {
  text: string;
  x: number;
  /** Taban çizgisi (üstten). */
  baseline: number;
  width: number;
  size: number;
  math: boolean;
  mono: boolean;
  bold: boolean;
  italic: boolean;
}

export interface Line {
  glyphs: Glyph[];
  x0: number;
  x1: number;
  top: number;
  bottom: number;
  baseline: number;
  size: number;
  text: string;
  mathRatio: number;
  monoRatio: number;
  boldRatio: number;
}

export type BlockKind =
  | 'text'
  | 'heading'
  | 'caption'
  | 'math'
  | 'code'
  | 'header'
  | 'footer'
  | 'page-number'
  | 'figure-text'
  | 'meta'
  | 'reference';

export interface Block {
  id: string;
  page: number;
  kind: BlockKind;
  x0: number;
  y0: number;
  x1: number;
  y1: number;
  /** Baskın yazı boyutu. */
  size: number;
  lines: Line[];
  /** Çeviri motoruna gidecek zengin metin (formüller/üst simgeler yer tutucu). */
  markup: string;
  plain: string;
  /** Yer tutucuların özgün metinleri (i="n" → atoms[n]). */
  atoms: { text: string; sup: boolean }[];
  bold: boolean;
  translatable: boolean;
  /** Okuma sırası. */
  order: number;
}

export interface PageLayout {
  index: number;
  width: number;
  height: number;
  blocks: Block[];
}

export interface LayoutOptions {
  skipReferences: boolean;
}

const MATH_FONT = /CMMI|CMSY|CMEX|CMBSY|MSAM|MSBM|LMMath|Math|Symbol|STIX|esint|rsfs|EUFM|EUSM|txsy|txmi|pxsy|pxmi|MTMI|MTSY|MTEX|Euclid|Cambria.?Math|XITS|Asana|wasy|stmary/i;
const MONO_FONT = /Mono|Courier|CMTT|Consol|Menlo|Inconsolata|SFMono|LucidaTypewriter|TeX-?Gyre-?Cursor|Typewriter/i;
const BOLD_FONT = /Bold|Black|Heavy|Semibold|Demi|CMBX|CMB\d|-B\b|Medi/i;
const ITALIC_FONT = /Italic|Oblique|CMTI|CMMI|-It\b|Slant/i;
const CAPTION = /^\s*(fig(ure|\.)?|table|tab\.|scheme|algorithm|chart|plate|supplementary\s+(fig(ure)?|table)|extended\s+data\s+(fig(ure)?|table))\s*[.:]?\s*[\dSIVX]/i;
const REF_HEADING = /^\s*(\d+(\.\d+)*\.?\s*)?(references?|bibliography|literature cited|works cited|cited literature|kaynak[çc]a|kaynaklar|referanslar)\s*:?\s*$/i;
const RESUME_HEADING = /^\s*([A-Z]\.?|\d+\.?)?\s*(appendix|appendices|supplementary|supporting information|ek\b|ekler)/i;

function mul(m1: number[], m2: number[]): number[] {
  return [
    m1[0] * m2[0] + m1[2] * m2[1],
    m1[1] * m2[0] + m1[3] * m2[1],
    m1[0] * m2[2] + m1[2] * m2[3],
    m1[1] * m2[2] + m1[3] * m2[3],
    m1[0] * m2[4] + m1[2] * m2[5] + m1[4],
    m1[1] * m2[4] + m1[3] * m2[5] + m1[5],
  ];
}

/** Ham öğeleri görünüm koordinatlarındaki (üst-sol) karakter gruplarına çevirir. */
export function toGlyphs(page: RawPage): Glyph[] {
  const out: Glyph[] = [];
  for (const item of page.items) {
    if (!item.str) continue;
    const tx = mul(page.viewportTransform, item.transform);
    // Döndürülmüş metinler (kenar filigranları, dikey etiketler) atlanır.
    if (Math.abs(tx[1]) > 0.05 * Math.abs(tx[0]) + 0.01 || Math.abs(tx[2]) > 0.05 * Math.abs(tx[3]) + 0.01) continue;
    const size = Math.hypot(tx[2], tx[3]);
    if (!size || size < 1) continue;
    const realName = page.fonts[item.fontName] ?? item.fontName;
    const scaleX = Math.hypot(page.viewportTransform[0], page.viewportTransform[1]) || 1;
    out.push({
      text: item.str,
      x: tx[4],
      baseline: tx[5],
      width: Math.abs(item.width * scaleX),
      size,
      math: MATH_FONT.test(realName),
      mono: MONO_FONT.test(realName),
      bold: BOLD_FONT.test(realName),
      italic: ITALIC_FONT.test(realName),
    });
  }
  return out;
}

/**
 * Sayfadaki dikey sütun boşluğunu (gutter) bulur: sayfanın orta bölgesinde, satırların büyük
 * çoğunluğunun hiç dokunmadığı en geniş dikey şerit. İki sütunlu makalelerde satırları bölmek için.
 */
export function findGutter(glyphs: Glyph[], pageWidth: number): { x0: number; x1: number } | null {
  const visible = glyphs.filter((g) => g.text.trim() && g.width > 0);
  if (visible.length < 40) return null;
  const W = Math.ceil(pageWidth);
  const cover = new Float32Array(W + 1);
  let total = 0;
  for (const g of visible) {
    const a = Math.max(0, Math.floor(g.x));
    const b = Math.min(W, Math.ceil(g.x + g.width));
    for (let x = a; x < b; x++) cover[x] += 1;
    total += 1;
  }
  const lo = Math.floor(pageWidth * 0.3);
  const hi = Math.ceil(pageWidth * 0.7);
  let peak = 0;
  for (let x = 0; x <= W; x++) peak = Math.max(peak, cover[x]);
  // Başlık ve yazar satırları ortadan geçebilir; bunları tolere etmek için eşik görece yüksek tutulur.
  const limit = peak * 0.15;
  let best: { x0: number; x1: number } | null = null;
  let start = -1;
  for (let x = lo; x <= hi; x++) {
    if (cover[x] <= limit) {
      if (start < 0) start = x;
    } else if (start >= 0) {
      if (!best || x - start > best.x1 - best.x0) best = { x0: start, x1: x };
      start = -1;
    }
  }
  if (start >= 0 && (!best || hi - start > best.x1 - best.x0)) best = { x0: start, x1: hi };
  if (!best || best.x1 - best.x0 < 6) return null;
  // Her iki yanda da yeterince metin olmalı (gerçekten iki sütun).
  let left = 0;
  let right = 0;
  for (let x = 0; x < best.x0; x++) left += cover[x];
  for (let x = best.x1; x <= W; x++) right += cover[x];
  if (left < total * 0.15 || right < total * 0.15) return null;
  return best;
}

/** Karakter gruplarını satırlara ayırır (sütun boşluklarında böler, üst/alt simgeleri satıra katar). */
export function buildLines(glyphs: Glyph[], gutter: { x0: number; x1: number } | null = null): Line[] {
  const sorted = glyphs.filter((g) => g.text.trim() || g.text === ' ').sort((a, b) => a.baseline - b.baseline || a.x - b.x);
  // 1) Taban çizgisine göre satır grupları
  const rows: Glyph[][] = [];
  for (const g of sorted) {
    if (!g.text.trim()) continue;
    let placed = false;
    for (let i = rows.length - 1; i >= Math.max(0, rows.length - 4); i--) {
      const row = rows[i];
      const ref = row.reduce((m, x) => (x.size > m.size ? x : m), row[0]);
      const tol = 0.42 * Math.min(ref.size, g.size);
      const supTol = 0.62 * ref.size;
      const dy = Math.abs(g.baseline - ref.baseline);
      const smaller = g.size < 0.86 * ref.size || ref.size < 0.86 * g.size;
      if (dy <= tol || (smaller && dy <= supTol)) {
        row.push(g);
        placed = true;
        break;
      }
    }
    if (!placed) rows.push([g]);
  }
  // 2) Her grubu yatay boşluklara göre satır parçalarına böl
  const lines: Line[] = [];
  for (const row of rows) {
    row.sort((a, b) => a.x - b.x);
    let current: Glyph[] = [];
    let lastX1 = -Infinity;
    const body = median(row.map((g) => g.size));
    for (const g of row) {
      const gap = g.x - lastX1;
      const crossesGutter = !!gutter && lastX1 <= gutter.x0 + 1 && g.x >= gutter.x1 - 1 && gap > 0.6 * body;
      if (current.length && (gap > Math.max(1.6 * body, 9) || crossesGutter)) {
        lines.push(makeLine(current));
        current = [];
      }
      current.push(g);
      lastX1 = Math.max(lastX1, g.x + g.width);
    }
    if (current.length) lines.push(makeLine(current));
  }
  return lines.sort((a, b) => a.top - b.top || a.x0 - b.x0);
}

function median(values: number[]): number {
  if (!values.length) return 0;
  const s = [...values].sort((a, b) => a - b);
  return s[Math.floor(s.length / 2)];
}

function makeLine(glyphs: Glyph[]): Line {
  const sizeByChars = new Map<number, number>();
  let chars = 0;
  let mathChars = 0;
  let monoChars = 0;
  let boldChars = 0;
  for (const g of glyphs) {
    const n = g.text.replace(/\s/g, '').length;
    chars += n;
    if (g.math) mathChars += n;
    if (g.mono) monoChars += n;
    if (g.bold) boldChars += n;
    const key = Math.round(g.size * 2) / 2;
    sizeByChars.set(key, (sizeByChars.get(key) ?? 0) + n);
  }
  const size = [...sizeByChars.entries()].sort((a, b) => b[1] - a[1])[0]?.[0] ?? glyphs[0].size;
  const main = glyphs.filter((g) => Math.abs(g.size - size) < 0.6);
  const baseline = median((main.length ? main : glyphs).map((g) => g.baseline));
  const x0 = Math.min(...glyphs.map((g) => g.x));
  const x1 = Math.max(...glyphs.map((g) => g.x + g.width));
  return {
    glyphs,
    x0,
    x1,
    top: baseline - size * 0.85,
    bottom: baseline + size * 0.25,
    baseline,
    size,
    text: glyphs.map((g) => g.text).join(''),
    mathRatio: chars ? mathChars / chars : 0,
    monoRatio: chars ? monoChars / chars : 0,
    boldRatio: chars ? boldChars / chars : 0,
  };
}

interface BlockDraft {
  lines: Line[];
  x0: number;
  x1: number;
  y0: number;
  y1: number;
  size: number;
  gap: number;
}

function overlapRatio(a0: number, a1: number, b0: number, b1: number): number {
  const ov = Math.min(a1, b1) - Math.max(a0, b0);
  const minW = Math.min(a1 - a0, b1 - b0);
  return minW > 0 ? ov / minW : 0;
}

/** Satırları paragraf bloklarına gruplar (sütunlar ayrı tutulur). */
export function buildBlocks(lines: Line[], pageWidth: number): BlockDraft[] {
  const blocks: BlockDraft[] = [];
  for (const line of lines) {
    const isDisplayMath = line.mathRatio > 0.55;
    let best: BlockDraft | null = null;
    let bestDist = Infinity;
    for (const b of blocks) {
      const last = b.lines[b.lines.length - 1];
      const dist = line.top - b.y1;
      const maxGap = Math.max(b.gap * 1.25, 0.9 * line.size);
      if (dist > maxGap || dist < -0.6 * line.size) continue;
      if (overlapRatio(line.x0, line.x1, b.x0, b.x1) < 0.35) continue;
      if (Math.abs(line.size - b.size) > 0.18 * b.size) continue;
      const lastMath = last.mathRatio > 0.55;
      if (isDisplayMath !== lastMath) continue;
      // Daktilo yazısı (e-posta, kod) ile normal metin ve birkaç sütuna yayılan satırlar ayrı tutulur.
      if (Math.abs(line.monoRatio - last.monoRatio) > 0.5) continue;
      const blockWidth = b.x1 - b.x0;
      if (line.x1 - line.x0 > blockWidth * 1.6 && blockWidth < pageWidth * 0.4) continue;
      // Paragraf başı girintisi veya önceki satırın erken bitmesi yeni paragraf demektir.
      const prevEndsSentence = /[.!?:;。]["'”’)]?\s*$/.test(last.text);
      const indented = line.x0 - b.x0 > 0.9 * line.size && line.x0 - b.x0 < 6 * line.size;
      const prevShort = last.x1 < b.x1 - 2.5 * line.size;
      if (prevEndsSentence && (indented || prevShort) && b.lines.length >= 1 && !isDisplayMath) continue;
      // Kalın tek satırlık başlık ile altındaki gövde metni ayrılır.
      if (last.boldRatio > 0.8 && line.boldRatio < 0.3 && b.lines.length === 1) continue;
      if (dist < bestDist) {
        best = b;
        bestDist = dist;
      }
    }
    if (best) {
      const last = best.lines[best.lines.length - 1];
      const step = line.top - last.top;
      if (best.lines.length === 1 && step > 0) best.gap = step - (last.bottom - last.top) + 0.4 * line.size;
      best.lines.push(line);
      best.x0 = Math.min(best.x0, line.x0);
      best.x1 = Math.max(best.x1, line.x1);
      best.y1 = Math.max(best.y1, line.bottom);
    } else {
      blocks.push({ lines: [line], x0: line.x0, x1: line.x1, y0: line.top, y1: line.bottom, size: line.size, gap: 0.6 * line.size });
    }
  }
  return blocks;
}

/** Blok satırlarını tek metne birleştirir; formül ve üst simgeleri yer tutucuya çevirir. */
function blockText(lines: Line[], bodySize: number): { markup: string; plain: string; atoms: Block['atoms'] } {
  const atoms: Block['atoms'] = [];
  let markup = '';
  let plain = '';
  lines.forEach((line, li) => {
    let lastX1: number | null = null;
    let lastMath = false;
    let mathRun: string[] = [];
    const flushMath = () => {
      if (!mathRun.length) return;
      const text = mathRun.join('').trim();
      mathRun = [];
      if (!text) return;
      const i = atoms.push({ text, sup: false }) - 1;
      markup += `<code class="notranslate" translate="no" i="${i}">${escapeHtml(text)}</code>`;
      plain += text;
    };
    for (const g of line.glyphs) {
      // Formül ile metin arasındaki geçişlerde boşluk daha dar olabilir.
      const threshold = lastX1 !== null && lastMath !== g.math ? 0.06 : 0.14;
      const endsWithSpace = mathRun.length ? /\s$/.test(mathRun[mathRun.length - 1]) : /\s$/.test(plain);
      const space = lastX1 !== null && g.x - lastX1 > threshold * g.size && !endsWithSpace && !/^\s/.test(g.text);
      const sup = g.size < 0.82 * line.size && g.baseline < line.baseline - 0.18 * line.size && /^[\d,\-–*†‡§\s]+$/.test(g.text.trim());
      if (g.math && !sup) {
        if (space) {
          if (mathRun.length) mathRun.push(' ');
          else {
            markup += ' ';
            plain += ' ';
          }
        }
        mathRun.push(g.text);
      } else {
        flushMath();
        if (space) {
          markup += ' ';
          plain += ' ';
        }
        if (sup && g.text.trim()) {
          const i = atoms.push({ text: g.text.trim(), sup: true }) - 1;
          markup += `<img i="${i}">`;
          plain += g.text.trim();
        } else {
          markup += escapeHtml(g.text);
          plain += g.text;
        }
      }
      lastX1 = g.x + g.width;
      lastMath = g.math;
    }
    flushMath();
    if (li < lines.length - 1) {
      // Satır sonu tirelemesini kaldır ("transcrip-\ntion" → "transcription").
      if (/[a-zçğıöşü]-$/i.test(plain) && /^\s*[a-zçğıöşü]/.test(lines[li + 1].text)) {
        markup = markup.replace(/-$/, '');
        plain = plain.replace(/-$/, '');
      } else if (!/\s$/.test(plain)) {
        markup += ' ';
        plain += ' ';
      }
    }
  });
  void bodySize;
  return { markup: markup.replace(/\s+/g, ' ').trim(), plain: plain.replace(/\s+/g, ' ').trim(), atoms };
}

function classify(draft: BlockDraft, text: string, page: { width: number; height: number }, bodySize: number): BlockKind {
  const lines = draft.lines;
  const mathRatio = lines.reduce((s, l) => s + l.mathRatio, 0) / lines.length;
  const monoRatio = lines.reduce((s, l) => s + l.monoRatio, 0) / lines.length;
  const boldRatio = lines.reduce((s, l) => s + l.boldRatio, 0) / lines.length;
  const trimmed = text.trim();
  const inTopBand = draft.y1 < page.height * 0.075;
  const inBottomBand = draft.y0 > page.height * 0.925;
  if (/^(page\s*)?\d{1,4}(\s*(of|\/)\s*\d{1,4})?$/i.test(trimmed) && (inTopBand || inBottomBand)) return 'page-number';
  const words = trimmed.match(/\p{L}{3,}/gu)?.length ?? 0;
  if (mathRatio > 0.55 || (!/\p{L}{3,}/u.test(trimmed) && /[=+\-−×·∑∫√≤≥]/.test(trimmed))) return 'math';
  // "Attention(Q, K, V) = softmax(QKᵀ/√d)V" gibi düz yazı tipli işlev adları içeren denklemler
  if (lines.length <= 3 && mathRatio > 0.1 && words <= 4 && /[=≤≥≈∝∈→⇒≡]/.test(trimmed)) return 'math';
  if (monoRatio > 0.8) return 'code';
  if (/[\w.+-]+@[\w-]+\.[\w.]+/.test(trimmed) && trimmed.length < 260) return 'meta';
  if (CAPTION.test(trimmed)) return 'caption';
  if ((inTopBand || inBottomBand) && trimmed.length < 160 && draft.size <= bodySize * 1.05) return inTopBand ? 'header' : 'footer';
  const short = trimmed.length < 60;
  if (draft.size < bodySize * 0.8 && short) return 'figure-text';
  if (lines.length <= 3 && short && (draft.size > bodySize * 1.12 || boldRatio > 0.7)) return 'heading';
  if (!/\p{L}{2,}/u.test(trimmed)) return 'figure-text';
  return 'text';
}

/** Blokları okuma sırasına koyar: tam genişlikli bloklar bant ayırıcıdır, bant içinde önce sol sütun. */
function orderBlocks(blocks: Block[], pageWidth: number): void {
  const mid = pageWidth / 2;
  const full = (b: Block) => b.x0 < mid - pageWidth * 0.06 && b.x1 > mid + pageWidth * 0.06;
  const sorted = [...blocks].sort((a, b) => a.y0 - b.y0 || a.x0 - b.x0);
  const ordered: Block[] = [];
  let band: Block[] = [];
  const flush = () => {
    const left = band.filter((b) => b.x0 + b.x1 < pageWidth).sort((a, b) => a.y0 - b.y0);
    const right = band.filter((b) => b.x0 + b.x1 >= pageWidth).sort((a, b) => a.y0 - b.y0);
    ordered.push(...left, ...right);
    band = [];
  };
  for (const b of sorted) {
    if (full(b)) {
      flush();
      ordered.push(b);
    } else band.push(b);
  }
  flush();
  ordered.forEach((b, i) => (b.order = i));
  blocks.sort((a, b) => a.order - b.order);
}

/** Bir sayfanın bloklarını çıkarır (belge düzeyindeki işaretlemeler analyzeDocument'ta yapılır). */
export function analyzePage(page: RawPage, bodySize?: number): PageLayout {
  const glyphs = toGlyphs(page);
  const lines = buildLines(glyphs, findGutter(glyphs, page.width));
  const size = bodySize ?? estimateBodySize(lines);
  const drafts = buildBlocks(lines, page.width);
  const blocks: Block[] = drafts.map((d, i) => {
    const { markup, plain, atoms } = blockText(d.lines, size);
    const kind = classify(d, plain, page, size);
    return {
      id: `${page.index}-${i}`,
      page: page.index,
      kind,
      x0: d.x0,
      y0: d.y0,
      x1: d.x1,
      y1: d.y1,
      size: d.size,
      lines: d.lines,
      markup,
      plain,
      atoms,
      bold: d.lines.every((l) => l.boldRatio > 0.7),
      translatable: kind === 'text' || kind === 'heading' || kind === 'caption',
      order: i,
    };
  });
  orderBlocks(blocks, page.width);
  return { index: page.index, width: page.width, height: page.height, blocks };
}

export function estimateBodySize(lines: Line[]): number {
  const counts = new Map<number, number>();
  for (const l of lines) {
    const key = Math.round(l.size * 2) / 2;
    counts.set(key, (counts.get(key) ?? 0) + l.text.length);
  }
  return [...counts.entries()].sort((a, b) => b[1] - a[1])[0]?.[0] ?? 10;
}

/**
 * Belge düzeyinde işaretleme: tekrarlanan üst/alt bilgileri ve kaynakça bölümünü bulur.
 * `pages` okuma sırasına göre dizilmiş olmalıdır.
 */
export function finalizeDocument(pages: PageLayout[], opts: LayoutOptions): void {
  // Tekrarlanan üst/alt bilgiler (rakamlar normalleştirilerek)
  const norm = (s: string) => s.toLowerCase().replace(/\d+/g, '#').replace(/\s+/g, ' ').trim();
  const freq = new Map<string, number>();
  for (const p of pages) {
    const seen = new Set<string>();
    for (const b of p.blocks) {
      if (b.y1 < p.height * 0.1 || b.y0 > p.height * 0.9) {
        const k = norm(b.plain);
        if (k.length > 2 && !seen.has(k)) {
          seen.add(k);
          freq.set(k, (freq.get(k) ?? 0) + 1);
        }
      }
    }
  }
  const minRepeat = Math.max(2, Math.ceil(pages.length * 0.4));
  for (const p of pages) {
    for (const b of p.blocks) {
      if ((b.y1 < p.height * 0.1 || b.y0 > p.height * 0.9) && (freq.get(norm(b.plain)) ?? 0) >= minRepeat && pages.length > 1) {
        b.kind = b.y1 < p.height * 0.1 ? 'header' : 'footer';
        b.translatable = false;
      }
    }
  }
  // İlk sayfada başlık ile "Abstract" arasındaki yazar/kurum blokları çevrilmez.
  const first = pages[0];
  if (first) {
    const abstractHead = first.blocks.find((b) => /^\s*(abstract|summary|özet)\b/i.test(b.plain));
    const title = [...first.blocks].filter((b) => b.y0 < first.height * 0.45 && (b.kind === 'heading' || b.kind === 'text')).sort((a, b) => b.size - a.size)[0];
    if (abstractHead && title && abstractHead.y0 > title.y1) {
      for (const b of first.blocks) {
        if (b === title || b === abstractHead) continue;
        if (b.y0 >= title.y1 - 1 && b.y1 <= abstractHead.y0 + 1 && b.plain.length < 300 && b.kind !== 'caption') {
          b.kind = 'meta';
          b.translatable = false;
        }
      }
    }
  }
  if (!opts.skipReferences) return;
  // Kaynakça: başlığından sonraki metin blokları (ek bölümüne kadar)
  let inRefs = false;
  for (const p of pages) {
    for (const b of p.blocks) {
      const t = b.plain.trim();
      if ((b.kind === 'heading' || b.kind === 'text') && t.length < 40 && REF_HEADING.test(t)) {
        inRefs = true;
        continue;
      }
      if (inRefs && b.kind === 'heading' && RESUME_HEADING.test(t)) inRefs = false;
      if (inRefs && (b.kind === 'text' || b.kind === 'heading')) {
        b.kind = 'reference';
        b.translatable = false;
      }
    }
  }
}

/** İlk sayfalardan başlık ve özeti tahmin eder (bağlam ve terim çıkarımı için). */
export function guessTitleAndAbstract(pages: PageLayout[]): { title: string; abstract: string } {
  const first = pages[0];
  if (!first) return { title: '', abstract: '' };
  const candidates = first.blocks.filter((b) => b.kind === 'heading' || b.kind === 'text');
  const title = [...candidates].filter((b) => b.y0 < first.height * 0.4).sort((a, b) => b.size - a.size)[0]?.plain ?? '';
  let abstract = '';
  const all = pages.slice(0, 2).flatMap((p) => p.blocks);
  const idx = all.findIndex((b) => /^\s*(abstract|summary|özet)\b/i.test(b.plain));
  if (idx >= 0) {
    const head = all[idx];
    abstract = head.plain.replace(/^\s*(abstract|summary|özet)\s*[:.—-]?\s*/i, '');
    if (abstract.length < 100 && all[idx + 1]) abstract += ' ' + all[idx + 1].plain;
  } else {
    abstract = candidates.filter((b) => b.kind === 'text' && b.plain.length > 300)[0]?.plain ?? '';
  }
  return { title: title.slice(0, 300), abstract: abstract.slice(0, 2500) };
}

import fontkit from '@pdf-lib/fontkit';
import { PDFDocument, rgb, type PDFFont, type PDFPage } from 'pdf-lib';
import type { Block, PageLayout } from './layout';
import { runsToPlain, translationRuns } from './text';

export type ExportMode = 'dual' | 'translated' | 'interleaved';

export interface ExportFonts {
  regular: Uint8Array;
  bold?: Uint8Array;
  math?: Uint8Array;
}

export interface ExportInput {
  original: Uint8Array;
  layouts: PageLayout[];
  /** blok kimliği → çeviri */
  translations: Map<string, string>;
  rich: boolean;
  mode: ExportMode;
  fonts: ExportFonts;
  /** Sayfa başlığına eklenecek bilgi notu (örn. "BioParrot · Google → Türkçe"). */
  credit?: string;
  onProgress?: (done: number, total: number) => void;
}

interface FontSet {
  regular: PDFFont;
  bold: PDFFont;
  math?: PDFFont;
  covers: Map<PDFFont, Set<number>>;
}

class TextMeasurer {
  private cache = new Map<string, number>();
  constructor(private fonts: FontSet) {}

  fontFor(ch: string, bold: boolean): PDFFont {
    const cp = ch.codePointAt(0)!;
    const primary = bold ? this.fonts.bold : this.fonts.regular;
    if (this.fonts.covers.get(primary)?.has(cp)) return primary;
    if (this.fonts.covers.get(this.fonts.regular)?.has(cp)) return this.fonts.regular;
    if (this.fonts.math && this.fonts.covers.get(this.fonts.math)?.has(cp)) return this.fonts.math;
    return primary;
  }

  /** Yazı tipinde olmayan karakterleri benzerleriyle değiştirir. */
  sanitize(text: string, bold: boolean): string {
    let out = '';
    for (const ch of text) {
      const font = this.fontFor(ch, bold);
      out += this.fonts.covers.get(font)?.has(ch.codePointAt(0)!) ? ch : /\s/.test(ch) ? ' ' : '?';
    }
    return out;
  }

  width(text: string, size: number, bold: boolean): number {
    let w = 0;
    for (const ch of text) {
      const key = `${bold ? 'b' : 'r'}${ch}`;
      let unit = this.cache.get(key);
      if (unit === undefined) {
        const font = this.fontFor(ch, bold);
        try {
          unit = font.widthOfTextAtSize(ch, 1);
        } catch {
          unit = 0.5;
        }
        this.cache.set(key, unit);
      }
      w += unit * size;
    }
    return w;
  }
}

/** Metni verilen genişliğe sarar. Boşluksuz yazılarda (Çince vb.) karakter bazında böler. */
function wrap(text: string, maxWidth: number, size: number, bold: boolean, m: TextMeasurer): string[] {
  const lines: string[] = [];
  for (const paragraph of text.split('\n')) {
    const tokens = /\s/.test(paragraph) ? paragraph.split(/(\s+)/) : Array.from(paragraph);
    let line = '';
    for (const token of tokens) {
      const candidate = line + token;
      if (m.width(candidate, size, bold) <= maxWidth || !line.trim()) {
        if (m.width(candidate, size, bold) > maxWidth && !line.trim()) {
          // Tek başına sığmayan uzun kelimeyi karakterlerden böl.
          let part = '';
          for (const ch of token) {
            if (m.width(part + ch, size, bold) > maxWidth && part) {
              lines.push(part);
              part = '';
            }
            part += ch;
          }
          line = part;
        } else {
          line = candidate;
        }
      } else {
        lines.push(line.trimEnd());
        line = token.trimStart();
      }
    }
    if (line.trim() || !lines.length) lines.push(line.trimEnd());
  }
  return lines;
}

interface Fitted {
  size: number;
  lines: string[];
  lineHeight: number;
}

/** Çeviriyi bloğun kutusuna sığacak en büyük yazı boyutunu bulur. */
function fit(text: string, block: Block, m: TextMeasurer, bold: boolean): Fitted {
  const width = Math.max(block.x1 - block.x0, 20);
  const height = Math.max(block.y1 - block.y0, block.size * 1.1);
  const start = Math.min(block.size, 30);
  const min = Math.max(start * 0.5, 4);
  let size = start;
  for (let i = 0; i < 14; i++) {
    const lineHeight = size * 1.16;
    const lines = wrap(text, width, size, bold, m);
    if (lines.length * lineHeight <= height + size * 0.35 || size <= min) return { size, lines, lineHeight };
    size = Math.max(min, size * 0.93);
  }
  return { size, lines: wrap(text, width, size, bold, m), lineHeight: size * 1.16 };
}

function drawBlock(page: PDFPage, block: Block, text: string, ox: number, oy: number, pageHeight: number, fonts: FontSet, m: TextMeasurer): void {
  const bold = block.kind === 'heading' || block.bold;
  const clean = m.sanitize(text, bold);
  const fitted = fit(clean, block, m, bold);
  const pad = 1.2;
  const boxTop = block.y0 - pad;
  const usedHeight = Math.max(block.y1 - block.y0 + pad * 2, fitted.lines.length * fitted.lineHeight + pad * 2);
  // Özgün metnin üzerini beyazla kapat.
  page.drawRectangle({
    x: ox + block.x0 - pad,
    y: oy + pageHeight - boxTop - usedHeight,
    width: block.x1 - block.x0 + pad * 2,
    height: usedHeight,
    color: rgb(1, 1, 1),
  });
  let baseline = block.y0 + fitted.size * 0.9;
  for (const line of fitted.lines) {
    let x = ox + block.x0;
    // Aynı yazı tipini kullanan karakterleri tek seferde çiz.
    let run = '';
    let runFont: PDFFont | null = null;
    const flush = () => {
      if (!run || !runFont) return;
      page.drawText(run, { x, y: oy + pageHeight - baseline, size: fitted.size, font: runFont, color: rgb(0.08, 0.08, 0.1) });
      x += runFont.widthOfTextAtSize(run, fitted.size);
      run = '';
    };
    for (const ch of line) {
      const f = m.fontFor(ch, bold);
      if (f !== runFont) {
        flush();
        runFont = f;
      }
      run += ch;
    }
    flush();
    baseline += fitted.lineHeight;
  }
}

/** Çevirileri özgün PDF'in üzerine yerleştirerek yeni bir PDF üretir. */
export async function buildTranslatedPdf(input: ExportInput): Promise<Uint8Array> {
  const src = await PDFDocument.load(input.original, { ignoreEncryption: true, updateMetadata: false });
  const out = await PDFDocument.create();
  out.registerFontkit(fontkit);
  // Not: fontkit'in alt kümeleme (subset) özelliği bazı yazı tiplerinde karakterleri kaybettiği için tam gömülür.
  const regular = await out.embedFont(input.fonts.regular, { subset: false });
  const covers = new Map<PDFFont, Set<number>>([[regular, new Set(regular.getCharacterSet())]]);
  const texts = [...input.translations.values()].join(' ');
  const needsBold = input.layouts.some((l) => l.blocks.some((b) => (b.kind === 'heading' || b.bold) && input.translations.has(b.id)));
  const bold = needsBold && input.fonts.bold ? await out.embedFont(input.fonts.bold, { subset: false }) : regular;
  if (bold !== regular) covers.set(bold, new Set(bold.getCharacterSet()));
  let math: PDFFont | undefined;
  const regularSet = covers.get(regular)!;
  const atomsText = input.layouts.flatMap((l) => l.blocks.flatMap((b) => b.atoms.map((a) => a.text))).join('');
  const missing = [...new Set(Array.from(texts + atomsText))].filter((ch) => !/\s/.test(ch) && !regularSet.has(ch.codePointAt(0)!));
  if (missing.length && input.fonts.math) {
    const candidate = fontkit.create(input.fonts.math as unknown as Buffer) as unknown as { hasGlyphForCodePoint(cp: number): boolean };
    if (missing.some((ch) => candidate.hasGlyphForCodePoint(ch.codePointAt(0)!))) {
      math = await out.embedFont(input.fonts.math, { subset: false });
      covers.set(math, new Set(math.getCharacterSet()));
    }
  }
  const fonts: FontSet = { regular, bold, math, covers };
  const measurer = new TextMeasurer(fonts);

  const total = src.getPageCount();
  const srcPages = src.getPages();
  for (let i = 0; i < total; i++) {
    const layout = input.layouts[i];
    const srcPage = srcPages[i];
    const crop = srcPage.getCropBox();
    const rotation = srcPage.getRotation().angle % 360;
    const blocks = (layout?.blocks ?? []).filter((b) => b.translatable && input.translations.has(b.id));
    const canOverlay = rotation === 0 && !!layout;
    const textOf = (b: Block) => runsToPlain(translationRuns(input.translations.get(b.id)!, input.rich, b));

    if (input.mode === 'dual') {
      const [embedded] = await out.embedPages([srcPage], [{ left: crop.x, bottom: crop.y, right: crop.x + crop.width, top: crop.y + crop.height }]);
      const page = out.addPage([crop.width * 2, crop.height]);
      page.drawPage(embedded, { x: 0, y: 0, width: crop.width, height: crop.height });
      page.drawPage(embedded, { x: crop.width, y: 0, width: crop.width, height: crop.height });
      if (canOverlay) for (const b of blocks) drawBlock(page, b, textOf(b), crop.width, 0, crop.height, fonts, measurer);
      page.drawLine({ start: { x: crop.width, y: 0 }, end: { x: crop.width, y: crop.height }, thickness: 0.5, color: rgb(0.8, 0.8, 0.82) });
    } else {
      if (input.mode === 'interleaved') {
        const [orig] = await out.copyPages(src, [i]);
        out.addPage(orig);
      }
      const [copy] = await out.copyPages(src, [i]);
      const page = out.addPage(copy);
      if (canOverlay) for (const b of blocks) drawBlock(page, b, textOf(b), crop.x, crop.y, crop.height, fonts, measurer);
    }
    input.onProgress?.(i + 1, total);
  }
  if (input.credit) out.setProducer(input.credit);
  out.setCreator('BioParrot');
  return out.save({ useObjectStreams: true });
}

import { describe, expect, it } from 'vitest';
import { analyzePage, buildLines, findGutter, finalizeDocument, type Glyph, type RawPage } from '@/core/pdf/layout';

/** Sentetik sayfa: her satır tek bir metin öğesi, soldan sağa (görünüm = PDF koordinatı, y aşağı doğru). */
function page(lines: { x: number; y: number; text: string; size?: number; font?: string }[], width = 600, height = 800): RawPage {
  const fonts: Record<string, string> = { R: 'ABCDEF+NimbusRomNo9L-Regu', B: 'ABCDEF+NimbusRomNo9L-Medi', M: 'ABCDEF+CMMI10' };
  return {
    index: 0,
    width,
    height,
    viewportTransform: [1, 0, 0, -1, 0, height],
    fonts,
    items: lines.map((l) => {
      const size = l.size ?? 10;
      return { str: l.text, transform: [size, 0, 0, size, l.x, height - l.y], width: l.text.length * size * 0.5, height: size, fontName: l.font ?? 'R' };
    }),
  };
}

describe('PDF layout', () => {
  it('detects the gutter of a two-column page and keeps columns apart', () => {
    const items = [];
    for (let i = 0; i < 30; i++) {
      items.push({ x: 50, y: 100 + i * 12, text: 'left column text that fills the line abc' });
      items.push({ x: 320, y: 100 + i * 12, text: 'right column text fills the line abcd' });
    }
    const raw = page(items);
    const layout = analyzePage(raw);
    const lefts = layout.blocks.filter((b) => b.x1 < 300);
    const rights = layout.blocks.filter((b) => b.x0 > 300);
    expect(lefts.length).toBeGreaterThan(0);
    expect(rights.length).toBeGreaterThan(0);
    expect(layout.blocks.every((b) => b.x1 < 300 || b.x0 > 300)).toBe(true);
    // Okuma sırası: önce sol sütun
    expect(layout.blocks[0].x0).toBe(50);
  });

  it('splits rows at the gutter even with a narrow gap', () => {
    const glyphs: Glyph[] = [];
    for (let i = 0; i < 50; i++) {
      glyphs.push({ text: 'x'.repeat(40), x: 50, baseline: 100 + i * 12, width: 230, size: 10, math: false, mono: false, bold: false, italic: false });
      glyphs.push({ text: 'y'.repeat(40), x: 292, baseline: 100 + i * 12, width: 230, size: 10, math: false, mono: false, bold: false, italic: false });
    }
    const gutter = findGutter(glyphs, 600);
    expect(gutter).not.toBeNull();
    const lines = buildLines(glyphs, gutter);
    expect(lines.length).toBe(100);
  });

  it('separates headings, removes hyphenation and keeps inline math as atoms', () => {
    const raw = page([
      { x: 50, y: 100, text: '1 Introduction', size: 12, font: 'B' },
      { x: 50, y: 130, text: 'Gene expression was measured in transcrip-' },
      { x: 50, y: 142, text: 'tional assays with' },
      { x: 150, y: 142, text: 'N', font: 'M' },
      { x: 160, y: 142, text: 'replicates per condition and sample.' },
    ]);
    const layout = analyzePage(raw, 10);
    expect(layout.blocks[0].kind).toBe('heading');
    const body = layout.blocks[1];
    expect(body.plain).toContain('transcriptional');
    expect(body.atoms.map((a) => a.text)).toEqual(['N']);
    expect(body.markup).toContain('<code class="notranslate" translate="no" i="0">N</code>');
  });

  it('marks the reference section as untranslatable', () => {
    const raw = page([
      { x: 50, y: 100, text: 'Discussion', size: 12, font: 'B' },
      { x: 50, y: 120, text: 'Our results indicate a strong association between the two variables.' },
      { x: 50, y: 160, text: 'References', size: 12, font: 'B' },
      { x: 50, y: 180, text: '[1] Smith J. The cell. Nature 2020; 580: 1-10.' },
    ]);
    const layout = analyzePage(raw, 10);
    finalizeDocument([layout], { skipReferences: true });
    const ref = layout.blocks.find((b) => b.plain.startsWith('[1]'))!;
    expect(ref.kind).toBe('reference');
    expect(ref.translatable).toBe(false);
    expect(layout.blocks.find((b) => b.plain.startsWith('Our results'))!.translatable).toBe(true);
  });
});

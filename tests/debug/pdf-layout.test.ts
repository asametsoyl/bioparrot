// @vitest-environment node
// Hata ayıklama: BP_PDF=/yol/makale.pdf npx vitest run tests/debug
import fs from 'node:fs';
import { test } from 'vitest';
import * as pdfjs from 'pdfjs-dist/legacy/build/pdf.mjs';
import { analyzePage, estimateBodySize, finalizeDocument, type RawPage } from '@/core/pdf/layout';

test.runIf(!!process.env.BP_PDF)('pdf layout dump', async () => {
  const data = new Uint8Array(fs.readFileSync(process.env.BP_PDF!));
  const doc = await pdfjs.getDocument({ data, fontExtraProperties: true }).promise;
  const raws: RawPage[] = [];
  const pages = Number(process.env.BP_PAGES ?? 3);
  for (let i = 0; i < Math.min(doc.numPages, pages); i++) {
    const page = await doc.getPage(i + 1);
    const vp = page.getViewport({ scale: 1 });
    const tc = await page.getTextContent();
    await page.getOperatorList();
    const items = tc.items.filter((x: any) => typeof x.str === 'string') as any[];
    const fonts: Record<string, string> = {};
    for (const it of items) {
      try { const f = page.commonObjs.get(it.fontName) as any; fonts[it.fontName] = (f?.name ?? it.fontName) + (f?.bold || f?.black ? '-Bold' : '') + (f?.italic ? '-Italic' : ''); } catch { fonts[it.fontName] = it.fontName; }
    }
    raws.push({ index: i, width: vp.width, height: vp.height, viewportTransform: vp.transform, items: items.map((x) => ({ str: x.str, transform: x.transform, width: x.width, height: x.height, fontName: x.fontName })), fonts });
  }
  const prelim = raws.map((r) => analyzePage(r));
  const body = estimateBodySize(prelim.flatMap((p) => p.blocks.flatMap((b) => b.lines)));
  const layouts = raws.map((r) => analyzePage(r, body));
  finalizeDocument(layouts, { skipReferences: true });
  const want = (process.env.BP_SHOW ?? '').split(',').filter(Boolean).map(Number);
  for (const l of layouts) {
    if (want.length && !want.includes(l.index + 1)) continue;
    console.log(`\n=== Page ${l.index + 1} (${l.width.toFixed(0)}x${l.height.toFixed(0)}) body=${body}`);
    for (const b of l.blocks) {
      console.log(`[${b.order}] ${b.kind.padEnd(11)} x=${b.x0.toFixed(0)}-${b.x1.toFixed(0)} y=${b.y0.toFixed(0)}-${b.y1.toFixed(0)} s=${b.size} bold=${b.bold ? 1 : 0} lines=${b.lines.length} | ${b.plain.slice(0, 90)}`);
      if (process.env.BP_GLYPHS) for (const ln of b.lines) if (ln.text.includes(process.env.BP_GLYPHS)) console.log(ln.glyphs.map((g) => `«${g.text}» x=${g.x.toFixed(1)} w=${g.width.toFixed(1)} end=${(g.x + g.width).toFixed(1)} s=${g.size.toFixed(1)} m=${g.math ? 1 : 0}`).join('\n'));
      if (process.env.BP_LINES) for (const ln of b.lines) console.log(`      line x=${ln.x0.toFixed(0)}-${ln.x1.toFixed(0)} top=${ln.top.toFixed(1)} s=${ln.size} bold=${ln.boldRatio.toFixed(2)} mono=${ln.monoRatio.toFixed(2)} math=${ln.mathRatio.toFixed(2)} | ${ln.text.slice(0, 70)}`);
    }
  }
  console.log('\nfonts:', [...new Set(raws.flatMap((r) => Object.values(r.fonts)))].join(', '));
}, 60000);

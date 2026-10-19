import { memo, useEffect, useLayoutEffect, useRef, useState } from 'react';
import type { Block, PageLayout } from '@/core/pdf/layout';
import { pdfjs, type PDFDocumentProxy } from '@/core/pdf/pdfjs';
import { translationRuns } from '@/core/pdf/text';
import type { PageStatus } from '@/core/pdf/job';

export type ViewMode = 'dual' | 'translated' | 'original';

interface Props {
  doc: PDFDocumentProxy;
  index: number;
  width: number;
  height: number;
  scale: number;
  view: ViewMode;
  layout?: PageLayout;
  translations: Map<string, string>;
  rich: boolean;
  status: PageStatus;
  version: number;
}

function BlockOverlay({ block, text, rich, scale }: { block: Block; text: string; rich: boolean; scale: number }) {
  const ref = useRef<HTMLDivElement>(null);
  const runs = translationRuns(text, rich, block);
  const heading = block.kind === 'heading' || block.bold;

  useLayoutEffect(() => {
    const el = ref.current;
    if (!el) return;
    const base = Math.min(block.size, 30) * scale;
    const min = base * 0.5;
    let fs = base;
    el.style.height = `${(block.y1 - block.y0) * scale + 2}px`;
    el.style.fontSize = `${fs}px`;
    for (let i = 0; i < 14 && el.scrollHeight > el.clientHeight + 1 && fs > min; i++) {
      fs = Math.max(min, fs * 0.92);
      el.style.fontSize = `${fs}px`;
    }
    if (el.scrollHeight > el.clientHeight + 1) el.style.height = 'auto';
  }, [text, scale, block]);

  return (
    <div
      ref={ref}
      className={`pv-block ${heading ? 'heading' : ''} ${block.kind}`}
      style={{
        left: block.x0 * scale - 1,
        top: block.y0 * scale - 1,
        width: (block.x1 - block.x0) * scale + 2,
      }}
    >
      {runs.map((r, i) =>
        r.text === '\n' ? (
          <br key={i} />
        ) : r.atom ? (
          r.sup ? (
            <sup key={i} className="atom">
              {r.text}
            </sup>
          ) : (
            <span key={i} className="atom">
              {r.text}
            </span>
          )
        ) : r.bold ? (
          <b key={i}>{r.text}</b>
        ) : r.italic ? (
          <i key={i}>{r.text}</i>
        ) : (
          <span key={i}>{r.text}</span>
        ),
      )}
    </div>
  );
}

/** Tek bir PDF sayfası: solda özgün sayfa (seçilebilir metinle), sağda çevirisi yerleştirilmiş sayfa. */
export const PageView = memo(function PageView({ doc, index, width, height, scale, view, layout, translations, rich, status, version }: Props) {
  const rootRef = useRef<HTMLDivElement>(null);
  const origCanvas = useRef<HTMLCanvasElement>(null);
  const transCanvas = useRef<HTMLCanvasElement>(null);
  const textLayerRef = useRef<HTMLDivElement>(null);
  const [visible, setVisible] = useState(false);
  const renderedScale = useRef(0);

  useEffect(() => {
    const el = rootRef.current;
    if (!el) return;
    const io = new IntersectionObserver((entries) => entries.forEach((e) => e.isIntersecting && setVisible(true)), { rootMargin: '1200px 0px' });
    io.observe(el);
    return () => io.disconnect();
  }, []);

  // Görünüm değişince yeni tuvaller boş gelir; yeniden çizilmeleri gerekir.
  useEffect(() => {
    renderedScale.current = 0;
  }, [view]);

  useEffect(() => {
    if (!visible || Math.abs(renderedScale.current - scale) < 0.01) return;
    let cancelled = false;
    let task: ReturnType<Awaited<ReturnType<PDFDocumentProxy['getPage']>>['render']> | null = null;
    void (async () => {
      const page = await doc.getPage(index + 1);
      if (cancelled) return;
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      const viewport = page.getViewport({ scale: scale * dpr });
      const canvas = origCanvas.current ?? transCanvas.current;
      if (!canvas) return;
      canvas.width = Math.floor(viewport.width);
      canvas.height = Math.floor(viewport.height);
      task = page.render({ canvas, viewport });
      try {
        await task.promise;
      } catch {
        return;
      }
      if (cancelled) return;
      // Çeviri tarafı aynı görüntüyü kullanır.
      const other = canvas === origCanvas.current ? transCanvas.current : null;
      if (other) {
        other.width = canvas.width;
        other.height = canvas.height;
        other.getContext('2d')!.drawImage(canvas, 0, 0);
      }
      renderedScale.current = scale;
      // Özgün metni seçilebilir yapan metin katmanı
      const container = textLayerRef.current;
      if (container) {
        container.replaceChildren();
        const textLayer = new pdfjs.TextLayer({
          textContentSource: page.streamTextContent(),
          container,
          viewport: page.getViewport({ scale }),
        });
        await textLayer.render().catch(() => {});
      }
    })();
    return () => {
      cancelled = true;
      task?.cancel();
    };
  }, [visible, scale, doc, index, view]);

  const w = width * scale;
  const h = height * scale;
  const blocks = layout?.blocks.filter((b) => b.translatable) ?? [];

  const original = (
    <div className="pv-page" style={{ width: w, height: h, ['--total-scale-factor' as string]: scale }}>
      <canvas ref={origCanvas} style={{ width: w, height: h }} />
      <div className="textLayer" ref={textLayerRef} />
    </div>
  );

  const translated = (
    <div className={`pv-page pv-translated ${status}`} style={{ width: w, height: h }}>
      <canvas ref={transCanvas} style={{ width: w, height: h }} />
      <div className="pv-overlay" data-version={version}>
        {blocks.map((b) => {
          const text = translations.get(b.id);
          if (text) return <BlockOverlay key={b.id} block={b} text={text} rich={rich} scale={scale} />;
          if (status === 'translating')
            return (
              <div
                key={b.id}
                className="pv-pending"
                style={{ left: b.x0 * scale, top: b.y0 * scale, width: (b.x1 - b.x0) * scale, height: (b.y1 - b.y0) * scale }}
              />
            );
          return null;
        })}
      </div>
    </div>
  );

  return (
    <div ref={rootRef} className={`pv-row view-${view}`} data-page={index + 1}>
      {view !== 'translated' && original}
      {view !== 'original' && translated}
      <div className="pv-num">{index + 1}</div>
    </div>
  );
});

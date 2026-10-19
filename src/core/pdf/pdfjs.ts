import * as pdfjs from 'pdfjs-dist';
import type { PDFDocumentProxy, PDFPageProxy } from 'pdfjs-dist';
import type { TextItem } from 'pdfjs-dist/types/src/display/api';
import type { RawPage } from './layout';

export type { PDFDocumentProxy, PDFPageProxy };
export { pdfjs };

let workerConfigured = false;

/** pdf.js işçi betiğinin adresini ayarlar (platforma göre değişir). */
export function configurePdfWorker(workerSrc: string): void {
  if (workerConfigured && pdfjs.GlobalWorkerOptions.workerSrc === workerSrc) return;
  pdfjs.GlobalWorkerOptions.workerSrc = workerSrc;
  workerConfigured = true;
}

export async function openPdf(data: Uint8Array): Promise<PDFDocumentProxy> {
  // pdf.js veriyi işçiye aktarırken kopyayı tüketir; özgün baytlar dışa aktarım için saklanır.
  // fontExtraProperties: yazı tipi adları (formül/kalın yazı tiplerini tanımak için) ana iş parçacığına da gelsin.
  return pdfjs.getDocument({ data: data.slice(), enableXfa: false, fontExtraProperties: true }).promise;
}

function isTextItem(item: unknown): item is TextItem {
  return !!item && typeof (item as TextItem).str === 'string' && Array.isArray((item as TextItem).transform);
}

/** Bir sayfanın metin öğelerini ve yazı tipi adlarını düzen analizine hazır biçimde okur. */
export async function readRawPage(doc: PDFDocumentProxy, index: number): Promise<RawPage> {
  const page = await doc.getPage(index + 1);
  const viewport = page.getViewport({ scale: 1 });
  const content = await page.getTextContent();
  const items = content.items.filter(isTextItem);
  // Yazı tipi nesneleri işlem listesi oluşturulunca çözülür; formül yazı tiplerini tanımak için gerçek adlar gerekir.
  try {
    await page.getOperatorList();
  } catch {
    /* bazı bozuk sayfalarda başarısız olabilir */
  }
  const fonts: Record<string, string> = {};
  for (const name of new Set(items.map((i) => i.fontName))) {
    try {
      const font = page.commonObjs.has(name) ? (page.commonObjs.get(name) as { name?: string; bold?: boolean; black?: boolean; italic?: boolean } | undefined) : undefined;
      // Bayraklar ada eklenir; düzen analizi ad üzerinden kalın/italik tespiti yapar.
      const flags = `${font?.bold || font?.black ? '-Bold' : ''}${font?.italic ? '-Italic' : ''}`;
      fonts[name] = (font?.name ?? content.styles[name]?.fontFamily ?? name) + flags;
    } catch {
      fonts[name] = content.styles[name]?.fontFamily ?? name;
    }
  }
  return {
    index,
    width: viewport.width,
    height: viewport.height,
    viewportTransform: viewport.transform,
    items: items.map((i) => ({ str: i.str, transform: i.transform, width: i.width, height: i.height, fontName: i.fontName, hasEOL: i.hasEOL })),
    fonts,
  };
}

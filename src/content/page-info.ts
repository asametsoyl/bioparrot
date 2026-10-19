/** Sayfanın başlık, özet ve makale meta verileri (yapay zekaya bağlam ve akademik mod için). */
export interface PageInfo {
  title: string;
  summary: string;
  /** Sayfa bir akademik makale mi (Highwire/Dublin Core meta verileri). */
  isPaper: boolean;
  pdfUrl?: string;
  doi?: string;
  htmlLang?: string;
  isPdf: boolean;
}

function meta(...names: string[]): string {
  for (const name of names) {
    const el = document.querySelector<HTMLMetaElement>(`meta[name="${name}" i], meta[property="${name}" i]`);
    const v = el?.content?.trim();
    if (v) return v;
  }
  return '';
}

const ABSTRACT_SELECTORS = [
  'blockquote.abstract',
  '#abstract',
  '#Abs1-content',
  '#abstracts',
  '#eng-abstract',
  '.abstract-content',
  'section.abstract',
  'div.abstract',
  '.ltx_abstract',
  '[class*="abstract" i]',
];

export function getPageInfo(): PageInfo {
  const citationTitle = meta('citation_title', 'dc.title', 'DC.Title', 'prism.title');
  const isPaper =
    !!citationTitle ||
    !!meta('citation_doi', 'citation_arxiv_id', 'citation_journal_title', 'citation_pdf_url', 'bepress_citation_title');
  const title = citationTitle || meta('og:title', 'twitter:title') || document.title || '';
  let summary = meta('citation_abstract', 'dc.description', 'DC.Description', 'og:description', 'description', 'twitter:description');
  if (isPaper && summary.length < 300) {
    for (const sel of ABSTRACT_SELECTORS) {
      const el = document.querySelector(sel);
      const text = el?.textContent?.replace(/\s+/g, ' ').trim();
      if (text && text.length > summary.length && text.length > 120) {
        summary = text;
        break;
      }
    }
  }
  return {
    title: title.replace(/\s+/g, ' ').trim(),
    summary: summary.replace(/\s+/g, ' ').trim().slice(0, 3000),
    isPaper,
    pdfUrl: meta('citation_pdf_url') || undefined,
    doi: meta('citation_doi', 'dc.identifier') || undefined,
    htmlLang: document.documentElement.lang || undefined,
    isPdf: document.contentType === 'application/pdf',
  };
}

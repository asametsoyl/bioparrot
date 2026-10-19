import { hostMatches } from '@/core/experts';
import type { SiteRule } from '@/core/types';
import { siteMatches } from '@/lib/settings';

export interface ResolvedRule {
  /** "Ana içerik" kapsamında yalnızca bu öğelerin içi çevrilir. */
  include: string[];
  /** Her kapsamda atlanan öğeler. */
  exclude: string[];
}

interface BuiltinRule {
  hosts: string[];
  include?: string[];
  exclude?: string[];
}

/** Yaygın akademik ve popüler siteler için yerleşik kurallar. */
const BUILTIN_RULES: BuiltinRule[] = [
  {
    hosts: ['arxiv.org', 'ar5iv.labs.arxiv.org', 'ar5iv.org'],
    include: ['h1.title', 'blockquote.abstract', '.ltx_page_main', '.ltx_document'],
    exclude: ['.ltx_bibliography', '.ltx_equation', '.ltx_equationgroup', '.ltx_eqn_table', '.ltx_tag', '.ltx_listing', '.ltx_note_mark', '.ltx_authors', '.ltx_role_footnote .ltx_tag', '.descriptor', '.ltx_page_footer', '.ltx_TOC'],
  },
  {
    hosts: ['pubmed.ncbi.nlm.nih.gov'],
    include: ['h1.heading-title', '#eng-abstract', '#abstract', '.abstract-content', '#other-abstract', '.keywords'],
    exclude: ['#references', '.references', '.similar-articles', '.citedby-articles', '.full-view-snippet'],
  },
  {
    hosts: ['pmc.ncbi.nlm.nih.gov', 'www.ncbi.nlm.nih.gov', 'europepmc.org'],
    include: ['article', 'main', '#main-content'],
    exclude: ['.ref-list', 'section.ref-list', '#reference-list', '.fm-authors', '.pmc-sidenav', '.citation-list', '#article--references'],
  },
  {
    hosts: ['biorxiv.org', 'medrxiv.org'],
    include: ['.highwire-cite-title', '#abstract-1', '.abstract', '.article.fulltext-view'],
    exclude: ['.ref-list', '#ref-list-1', '.highwire-cite-authors', '.section.ref-list'],
  },
  {
    hosts: ['nature.com'],
    include: ['h1.c-article-title', '.c-article-section__content', '#Abs1-content', 'figcaption'],
    exclude: ['#references', '.c-article-references', '.c-article-references__links', '.c-article-author-list'],
  },
  {
    hosts: ['sciencedirect.com'],
    include: ['.title-text', '#abstracts', '#body', '.Abstracts', '.keywords-section'],
    exclude: ['.bibliography', '#references', '.References', '.author-group'],
  },
  {
    hosts: ['link.springer.com', 'springer.com', 'biomedcentral.com', 'springeropen.com'],
    include: ['h1.c-article-title', '.c-article-section__content', '#Abs1-content'],
    exclude: ['.c-article-references', '#Bib1-content', '.c-article-author-list'],
  },
  {
    hosts: ['onlinelibrary.wiley.com', 'febs.onlinelibrary.wiley.com'],
    include: ['.citation__title', '.article-section__content', '.article__body'],
    exclude: ['.article-section__references', '.references', '.loa-wrapper'],
  },
  {
    hosts: ['journals.plos.org'],
    include: ['#artTitle', '.abstract-content', '#artText'],
    exclude: ['.references', '.ref-tip'],
  },
  {
    hosts: ['frontiersin.org'],
    include: ['.JournalAbstract', '.JournalFullText', 'h1'],
    exclude: ['.References', '.authors'],
  },
  {
    hosts: ['mdpi.com'],
    include: ['h1.title', '.art-abstract', '.html-body', '.html-abstract'],
    exclude: ['.html-references', '.html-ref', '.art-authors'],
  },
  {
    hosts: ['cell.com', 'thelancet.com'],
    include: ['.article-header__title', '.section-paragraph', '.article__body', '.abstract'],
    exclude: ['.references', '.ref__title', '.article-header__authors'],
  },
  {
    hosts: ['elifesciences.org'],
    include: ['.content-header__title', '.article-section__body', '.abstract'],
    exclude: ['.reference-list', '.author-list'],
  },
  {
    hosts: ['wikipedia.org'],
    include: ['#firstHeading', '#mw-content-text'],
    exclude: ['.reference', '.reflist', '.mw-references-wrap', '.mw-editsection', '.navbox', '.metadata', '#toc', '.toc', '.mw-cite-backlink', '.hatnote'],
  },
  {
    hosts: ['twitter.com', 'x.com'],
    include: ['[data-testid="tweetText"]', '[data-testid="UserDescription"]', '[data-testid="card.layoutLarge.detail"]', '[data-testid="card.layoutSmall.detail"]'],
  },
  {
    hosts: ['reddit.com'],
    include: ['shreddit-post [slot="title"]', 'shreddit-post [slot="text-body"]', 'shreddit-comment [slot="comment"]', '[slot="post-media-container"] + div', '.md'],
  },
  {
    hosts: ['github.com'],
    include: ['.markdown-body', '.comment-body', '.js-issue-title', '.markdown-title', '[data-testid="issue-title"]', 'p.f4'],
    exclude: ['pre', '.blob-code', '.react-code-lines', '.highlight', '.commit-ref'],
  },
  {
    hosts: ['stackoverflow.com', 'stackexchange.com', 'superuser.com', 'serverfault.com', 'askubuntu.com'],
    include: ['#question-header h1', '.s-prose', '.comment-copy'],
    exclude: ['pre'],
  },
  {
    hosts: ['news.ycombinator.com'],
    include: ['.titleline > a', '.commtext', '.toptext'],
  },
  {
    hosts: ['youtube.com'],
    include: ['#content-text', 'h1.ytd-watch-metadata', '#description-inline-expander', '#title.ytd-video-primary-info-renderer'],
    exclude: ['#chat', '#guide', 'ytd-masthead'],
  },
  {
    hosts: ['scholar.google.com'],
    include: ['.gs_rt', '.gs_rs'],
  },
  {
    hosts: ['semanticscholar.org'],
    include: ['[data-test-id="paper-detail-title"]', '[data-test-id="abstract-text"]', '.tldr-abstract-replacement', '.cl-paper-abstract'],
  },
  {
    hosts: ['openreview.net'],
    include: ['.note_content_title', '.note-content-value', '.citation_title'],
  },
];

/** Tüm sitelerde atlanan öğeler. */
export const GLOBAL_EXCLUDE = [
  '[translate="no"]',
  '.notranslate',
  '[contenteditable=""]',
  '[contenteditable="true"]',
  '[contenteditable="plaintext-only"]',
  'math',
  'mjx-container',
  '.MathJax',
  '.MathJax_Display',
  '.MathJax_Preview',
  '.katex-display',
  '.mwe-math-element',
  '.ltx_equation',
  '.ltx_equationgroup',
];

/** Akademik "kaynakçayı çevirme" ayarı için kaynakça kapsayıcıları. */
export const REFERENCE_SELECTORS = [
  '#references',
  '#bibliography',
  '#reference-list',
  '#ref-list',
  '.references',
  '.bibliography',
  '.ref-list',
  '.reflist',
  '.ltx_bibliography',
  '.c-article-references',
  '.article-section__references',
  '.citation-list',
  '.html-references',
  '#Bib1-content',
  'ol.references',
  'ul.references',
  '[data-title="References"]',
  '[aria-labelledby="references"]',
];

export const REFERENCE_HEADING = /^\s*(\d+(\.\d+)*\.?\s*)?(references?|bibliography|literature cited|works cited|cited literature|kaynak[çc]a|kaynaklar|referanslar|literat[üu]r|bibliografie|literatur|références|bibliografía|referências)\s*:?\s*$/i;

function splitSelectors(s: string | undefined): string[] {
  return (s ?? '')
    .split(/[,\n]/)
    .map((x) => x.trim())
    .filter(Boolean);
}

export function resolveSiteRule(host: string, custom: SiteRule[]): ResolvedRule {
  const include: string[] = [];
  const exclude: string[] = [];
  for (const rule of BUILTIN_RULES) {
    if (hostMatches(host, rule.hosts)) {
      include.push(...(rule.include ?? []));
      exclude.push(...(rule.exclude ?? []));
    }
  }
  const user = custom.filter((r) => r.enabled && siteMatches(host, [r.pattern]));
  const userInclude = user.flatMap((r) => splitSelectors(r.include));
  const userExclude = user.flatMap((r) => splitSelectors(r.exclude));
  return {
    include: userInclude.length ? userInclude : include,
    exclude: [...exclude, ...userExclude],
  };
}

/** Geçersiz CSS seçicileri atarak güvenli bir birleşik seçici üretir. */
export function safeSelector(list: string[]): string | null {
  const valid = list.filter((sel) => {
    try {
      document.createDocumentFragment().querySelector(sel);
      return true;
    } catch {
      return false;
    }
  });
  return valid.length ? valid.join(',') : null;
}

import type { Expert, Settings } from './types';

/**
 * Yerleşik "AI uzmanları": alan bazlı çeviri istemleri.
 * İstemler İngilizce yazılır (modeller en tutarlı sonucu böyle verir); çıktı dili {{to}} ile belirlenir.
 * Biçim kuralları (parça etiketleri, yer tutucular) istem oluşturucuda ayrıca eklenir.
 */

const ACADEMIC_CORE = `Follow the conventions of scholarly writing in {{to}}:
- Use a formal, precise academic register. Keep the author's line of argument and sentence boundaries close to the original so the translation can be read side by side with the source.
- Translate technical terms with the equivalents established in the {{to}} scientific literature. Use the same translation for the same term throughout.
- Never change citations (e.g. [12], [3–5], (Smith et al., 2020)), cross-references (Fig. 2, Table 1, Eq. 3, Section 4.2), numbers, units, variables, symbols, mathematical expressions, chemical formulas or code.
- Keep abbreviations and acronyms as they appear in the source; expand them only if the source does.
- Preserve hedging and the strength of claims exactly (may, might, suggests, is associated with, significantly). Do not add, omit or interpret content.`;

export const BUILTIN_EXPERTS: Expert[] = [
  {
    id: 'academic',
    group: 'academic',
    icon: '🎓',
    name: 'Akademik (genel)',
    nameEn: 'Academic (general)',
    description: 'Makale, tez ve bildiriler için resmi ve terim odaklı çeviri.',
    descriptionEn: 'Formal, terminology-aware translation for papers, theses and proceedings.',
    temperature: 0.2,
    systemPrompt: `You are an expert academic translator who translates scholarly literature into {{to}} for researchers.
${ACADEMIC_CORE}`,
  },
  {
    id: 'biomedical',
    group: 'academic',
    icon: '🧬',
    name: 'Biyoloji ve yaşam bilimleri',
    nameEn: 'Biology & life sciences',
    description: 'Moleküler biyoloji, genetik, mikrobiyoloji, ekoloji ve biyoinformatik.',
    descriptionEn: 'Molecular biology, genetics, microbiology, ecology and bioinformatics.',
    temperature: 0.2,
    systemPrompt: `You are an expert translator of biological and biomedical research literature into {{to}}, with a PhD-level understanding of molecular biology, genetics, cell biology, microbiology, ecology and bioinformatics.
${ACADEMIC_CORE}
Domain rules:
- Keep gene and protein symbols exactly as written (e.g. TP53, BRCA1, p53, IL-6, CRISPR-Cas9), including their capitalisation.
- Keep binomial species names and strain designations unchanged (e.g. Escherichia coli K-12, Arabidopsis thaliana, Mus musculus).
- Keep accession numbers, sequence notation (5′→3′), restriction enzymes, cell line names, reagents, kits and software names unchanged.
- Use the standard {{to}} terminology for anatomy, physiology, cell biology and laboratory methods; keep statistics exactly as written (p < 0.05, 95% CI, n = 12, SD, SEM).`,
  },
  {
    id: 'medicine',
    group: 'academic',
    icon: '🩺',
    name: 'Tıp ve farmakoloji',
    nameEn: 'Medicine & pharmacology',
    description: 'Klinik çalışmalar, vaka raporları, farmakoloji ve halk sağlığı.',
    descriptionEn: 'Clinical studies, case reports, pharmacology and public health.',
    temperature: 0.2,
    systemPrompt: `You are an expert medical translator who translates clinical and pharmacological literature into {{to}} for physicians and researchers.
${ACADEMIC_CORE}
Domain rules:
- Use the official medical terminology of {{to}} for diseases, anatomy, symptoms, procedures and outcomes; follow ICD/MeSH concepts.
- Use International Nonproprietary Names (INN) for drugs; keep brand names, doses, routes and regimens exactly as written.
- Keep trial names, registration numbers, scales and scores (e.g. NYHA class II, HbA1c, GCS 13) and all statistics (HR, OR, RR, 95% CI, p values) unchanged.
- Preserve clinical nuance precisely (e.g. "associated with" vs "caused", "non-inferior", "statistically significant").`,
  },
  {
    id: 'chemistry',
    group: 'academic',
    icon: '⚗️',
    name: 'Kimya ve malzeme',
    nameEn: 'Chemistry & materials',
    description: 'Organik/anorganik kimya, biyokimya, malzeme bilimi.',
    descriptionEn: 'Organic/inorganic chemistry, biochemistry, materials science.',
    temperature: 0.2,
    systemPrompt: `You are an expert translator of chemistry, biochemistry and materials-science literature into {{to}}.
${ACADEMIC_CORE}
Domain rules:
- Keep chemical formulas, IUPAC names, CAS numbers, reaction schemes, oxidation states, isotopes and units unchanged.
- Use the IUPAC-recommended {{to}} nomenclature for common compounds and the established {{to}} terms for techniques (NMR, XRD, HPLC, DFT) and reaction types.`,
  },
  {
    id: 'physics-math',
    group: 'academic',
    icon: '∑',
    name: 'Fizik ve matematik',
    nameEn: 'Physics & mathematics',
    description: 'Teorem, ispat ve formül yoğun metinler.',
    descriptionEn: 'Theorem-, proof- and formula-heavy texts.',
    temperature: 0.1,
    systemPrompt: `You are an expert translator of physics and mathematics papers into {{to}}.
${ACADEMIC_CORE}
Domain rules:
- Leave every formula, LaTeX command, variable and symbol untouched; translate only the surrounding prose.
- Use the standard {{to}} terms for structures and statements (theorem, lemma, corollary, proof, without loss of generality, if and only if).`,
  },
  {
    id: 'cs-engineering',
    group: 'academic',
    icon: '💻',
    name: 'Bilgisayar ve mühendislik',
    nameEn: 'Computer science & engineering',
    description: 'Yapay zeka, yazılım, elektronik ve mühendislik makaleleri.',
    descriptionEn: 'AI, software, electronics and engineering papers.',
    temperature: 0.2,
    systemPrompt: `You are an expert translator of computer-science and engineering research into {{to}}.
${ACADEMIC_CORE}
Domain rules:
- Keep code, identifiers, API and library names, model names, dataset names, metrics (F1, BLEU, AUC) and hyperparameters unchanged.
- Where the {{to}} community commonly uses the English term (e.g. benchmark, fine-tuning, transformer), prefer the established usage.`,
  },
  {
    id: 'social-humanities',
    group: 'academic',
    icon: '📚',
    name: 'Sosyal ve beşeri bilimler',
    nameEn: 'Social sciences & humanities',
    description: 'Psikoloji, sosyoloji, tarih, felsefe, eğitim, dilbilim.',
    descriptionEn: 'Psychology, sociology, history, philosophy, education, linguistics.',
    temperature: 0.3,
    systemPrompt: `You are an expert translator of social-science and humanities scholarship into {{to}}.
${ACADEMIC_CORE}
Domain rules:
- Render theoretical concepts with the terms used in the {{to}} literature of the discipline; when a concept is usually cited in its original language, keep it and add the {{to}} gloss in parentheses on first use.
- Keep quotations faithful in meaning and style; keep the names of authors, schools of thought and works as they are usually written in {{to}}.`,
  },
  {
    id: 'economics',
    group: 'academic',
    icon: '📈',
    name: 'Ekonomi ve işletme',
    nameEn: 'Economics & business',
    description: 'İktisat, finans, yönetim ve pazarlama literatürü.',
    descriptionEn: 'Economics, finance, management and marketing literature.',
    temperature: 0.2,
    systemPrompt: `You are an expert translator of economics, finance and management research into {{to}}.
${ACADEMIC_CORE}
Domain rules:
- Use the standard {{to}} terms for economic and econometric concepts; keep model names, indices, tickers, currencies and figures unchanged.`,
  },
  {
    id: 'law',
    group: 'academic',
    icon: '⚖️',
    name: 'Hukuk',
    nameEn: 'Law',
    description: 'Mevzuat, içtihat ve hukuk makaleleri.',
    descriptionEn: 'Legislation, case law and legal scholarship.',
    temperature: 0.1,
    systemPrompt: `You are an expert legal translator who translates legal texts and scholarship into {{to}}.
${ACADEMIC_CORE}
Domain rules:
- Use the precise legal terminology of the {{to}} legal system; where concepts have no equivalent, keep the original term in parentheses.
- Keep case names, statute and article numbers, and citations unchanged.`,
  },
  {
    id: 'academic-writing',
    group: 'academic',
    icon: '✍️',
    name: 'Akademik yazım (makale dili)',
    nameEn: 'Academic writing (journal style)',
    description: 'Kendi taslağınızı dergiye uygun akademik dile çevirir (giriş kutusunda kullanışlı).',
    descriptionEn: 'Turns your own draft into journal-ready academic prose (handy in input boxes).',
    temperature: 0.3,
    systemPrompt: `You are an experienced academic editor and translator. Translate the author's draft into clear, idiomatic academic {{to}} suitable for submission to a peer-reviewed journal.
- Improve clarity, concision and flow while keeping the author's meaning, claims and level of certainty unchanged.
- Use precise discipline-appropriate terminology and avoid colloquial expressions.
- Keep citations, numbers, units, symbols and abbreviations unchanged.`,
  },
  {
    id: 'general',
    group: 'general',
    icon: '🌐',
    name: 'Genel',
    nameEn: 'General',
    description: 'Her tür metin için akıcı ve doğal çeviri.',
    descriptionEn: 'Fluent, natural translation for any text.',
    temperature: 0.3,
    systemPrompt: `You are a professional translator. Translate the text into fluent, natural {{to}} that reads as if it had been written in {{to}}, while preserving the meaning, tone and every fact of the original exactly.`,
  },
  {
    id: 'tech',
    group: 'general',
    icon: '🛠️',
    name: 'Teknik doküman',
    nameEn: 'Technical docs',
    description: 'Yazılım dokümanları, README, API ve geliştirici içerikleri.',
    descriptionEn: 'Software docs, READMEs, APIs and developer content.',
    temperature: 0.2,
    systemPrompt: `You are a senior technical translator for software documentation. Translate into {{to}} using the terminology {{to}}-speaking developers actually use. Keep code, commands, identifiers, file paths, product names and UI labels in backticks unchanged.`,
  },
  {
    id: 'news',
    group: 'general',
    icon: '📰',
    name: 'Haber',
    nameEn: 'News',
    description: 'Gazetecilik üslubunda, tarafsız ve net çeviri.',
    descriptionEn: 'Neutral, clear translation in journalistic style.',
    temperature: 0.3,
    systemPrompt: `You are a news translator. Translate into {{to}} in a neutral, concise journalistic style. Keep names, quotes, figures and dates accurate; use the forms of place and person names customary in {{to}} media.`,
  },
  {
    id: 'social',
    group: 'general',
    icon: '💬',
    name: 'Sosyal medya',
    nameEn: 'Social media',
    description: 'Twitter/X, Reddit gibi platformlar için gündelik ve doğal dil.',
    descriptionEn: 'Casual, natural language for platforms like X and Reddit.',
    temperature: 0.5,
    systemPrompt: `You translate social-media posts into {{to}} the way a native {{to}} speaker would write them online: casual, natural and concise. Render slang, memes and idioms with natural {{to}} equivalents. Keep @mentions, #hashtags, URLs, emoji and tickers unchanged.`,
  },
  {
    id: 'literature',
    group: 'general',
    icon: '📖',
    name: 'Edebiyat',
    nameEn: 'Literature',
    description: 'Roman ve öykülerde üslubu ve duyguyu koruyan çeviri.',
    descriptionEn: 'Preserves style and emotion in fiction.',
    temperature: 0.6,
    systemPrompt: `You are a literary translator. Translate into {{to}}, preserving the author's voice, rhythm, imagery and emotional tone. Prefer natural literary {{to}} over word-for-word renderings, but never add or remove content.`,
  },
  {
    id: 'subtitle',
    group: 'general',
    icon: '🎬',
    name: 'Video altyazısı',
    nameEn: 'Video subtitles',
    description: 'Kısa, konuşma diline uygun altyazı çevirisi.',
    descriptionEn: 'Short, spoken-style subtitle translation.',
    temperature: 0.3,
    systemPrompt: `You translate video subtitles into {{to}}. Keep each line short and easy to read, use natural spoken {{to}}, and keep the meaning of every line. Do not merge or split lines.`,
  },
  {
    id: 'simple',
    group: 'general',
    icon: '🌱',
    name: 'Sade dil',
    nameEn: 'Plain language',
    description: 'Karmaşık metinleri anlaşılır, sade bir dille çevirir.',
    descriptionEn: 'Translates complex text into plain, easy language.',
    temperature: 0.4,
    systemPrompt: `Translate into plain, easy-to-understand {{to}}. Use short sentences and everyday words, and briefly explain unavoidable technical terms in parentheses, without losing any information.`,
  },
];

const BUILTIN_BY_ID = new Map(BUILTIN_EXPERTS.map((e) => [e.id, e]));

/** Yerleşik + kullanıcı uzmanlarını (yerleşik değişiklikler uygulanmış olarak) döner. */
export function getAllExperts(settings: Pick<Settings, 'experts' | 'expertOverrides'>): Expert[] {
  const builtins = BUILTIN_EXPERTS.map((e) => {
    const o = settings.expertOverrides?.[e.id];
    return o ? { ...e, ...o, builtin: true } : { ...e, builtin: true };
  });
  return [...builtins, ...(settings.experts ?? []).map((e) => ({ ...e, builtin: false }))];
}

export function getExpert(settings: Pick<Settings, 'experts' | 'expertOverrides'>, id: string | undefined): Expert {
  const all = getAllExperts(settings);
  return all.find((e) => e.id === id) ?? all.find((e) => e.id === 'general')!;
}

export function isBuiltinExpert(id: string): boolean {
  return BUILTIN_BY_ID.has(id);
}

export function expertName(expert: Expert, ui: 'tr' | 'en'): string {
  return ui === 'en' && expert.nameEn ? expert.nameEn : expert.name;
}

export function expertDescription(expert: Expert, ui: 'tr' | 'en'): string {
  return (ui === 'en' && expert.descriptionEn ? expert.descriptionEn : expert.description) ?? '';
}

/** Akademik yayıncı ve veri tabanı alan adları. */
export const ACADEMIC_HOSTS = [
  'arxiv.org',
  'ar5iv.labs.arxiv.org',
  'ar5iv.org',
  'alphaxiv.org',
  'pubmed.ncbi.nlm.nih.gov',
  'pmc.ncbi.nlm.nih.gov',
  'ncbi.nlm.nih.gov',
  'europepmc.org',
  'biorxiv.org',
  'medrxiv.org',
  'chemrxiv.org',
  'ssrn.com',
  'researchsquare.com',
  'preprints.org',
  'nature.com',
  'science.org',
  'cell.com',
  'thelancet.com',
  'nejm.org',
  'jamanetwork.com',
  'bmj.com',
  'sciencedirect.com',
  'link.springer.com',
  'springer.com',
  'onlinelibrary.wiley.com',
  'tandfonline.com',
  'journals.sagepub.com',
  'academic.oup.com',
  'journals.plos.org',
  'mdpi.com',
  'frontiersin.org',
  'pubs.acs.org',
  'pubs.rsc.org',
  'ieeexplore.ieee.org',
  'dl.acm.org',
  'aclanthology.org',
  'openreview.net',
  'proceedings.neurips.cc',
  'proceedings.mlr.press',
  'jstor.org',
  'cambridge.org',
  'elifesciences.org',
  'embopress.org',
  'pnas.org',
  'royalsocietypublishing.org',
  'iopscience.iop.org',
  'journals.aps.org',
  'semanticscholar.org',
  'scholar.google.com',
  'researchgate.net',
  'dergipark.org.tr',
  'tez.yok.gov.tr',
];

const SOCIAL_HOSTS = ['twitter.com', 'x.com', 'reddit.com', 'facebook.com', 'instagram.com', 'threads.net', 'bsky.app', 'mastodon.social', 'linkedin.com', 'quora.com'];
const TECH_HOSTS = ['github.com', 'gitlab.com', 'stackoverflow.com', 'stackexchange.com', 'developer.mozilla.org', 'readthedocs.io', 'docs.python.org', 'npmjs.com', 'pypi.org', 'huggingface.co', 'dev.to', 'medium.com'];
const NEWS_HOSTS = ['bbc.com', 'bbc.co.uk', 'nytimes.com', 'cnn.com', 'reuters.com', 'theguardian.com', 'apnews.com', 'washingtonpost.com', 'bloomberg.com', 'ft.com', 'economist.com', 'aljazeera.com', 'dw.com', 'lemonde.fr', 'spiegel.de'];

export function hostMatches(host: string, domains: string[]): boolean {
  const h = host.toLowerCase().replace(/^www\./, '');
  return domains.some((d) => h === d || h.endsWith('.' + d));
}

/**
 * Siteye göre önerilen uzman. Sayfada makale meta verisi (citation_title) varsa akademik sayılır.
 * Dönüş undefined ise kullanıcının varsayılan uzmanı kullanılır.
 */
export function suggestExpertForPage(
  host: string,
  isPaper: boolean,
  academicExpertId: string,
): string | undefined {
  if (isPaper || hostMatches(host, ACADEMIC_HOSTS)) return academicExpertId;
  if (hostMatches(host, SOCIAL_HOSTS)) return 'social';
  if (hostMatches(host, TECH_HOSTS)) return 'tech';
  if (hostMatches(host, NEWS_HOSTS)) return 'news';
  return undefined;
}

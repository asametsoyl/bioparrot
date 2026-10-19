/** Uygulama genelindeki paylaşılan tipler (eklenti ve Zotero eklentisi ortak kullanır). */

export type UiLocalePref = 'auto' | 'tr' | 'en';
export type UiLocale = 'tr' | 'en';
export type ThemePref = 'system' | 'light' | 'dark';
export type DisplayMode = 'bilingual' | 'translation';
export type TranslationScope = 'main' | 'full';
export type HoverKey = 'Shift' | 'Alt' | 'Control' | 'Meta' | 'off';
export type SelectionMode = 'icon' | 'direct' | 'off';
export type SelectionDisplay = 'inline' | 'card';
export type InputTrigger = 'triple-space' | 'double-space' | 'off';

export type EngineType =
  | 'google'
  | 'microsoft'
  | 'deepl'
  | 'deeplx'
  | 'openai'
  | 'anthropic'
  | 'gemini';

export interface EngineConfig {
  id: string;
  type: EngineType;
  name: string;
  apiKey?: string;
  baseUrl?: string;
  model?: string;
  /** Microsoft (Azure) Translator bölgesi. */
  region?: string;
  temperature?: number;
  /** Claude düşünme derinliği (low/medium/high); boşsa gönderilmez. */
  effort?: '' | 'low' | 'medium' | 'high';
  /** Aynı anda en fazla kaç istek. */
  maxConcurrent: number;
  /** Saniyedeki en fazla istek. */
  maxRps: number;
  /** Tek istekte en fazla kaç parça. */
  batchSize: number;
  /** Tek istekte en fazla kaç karakter. */
  batchChars: number;
  /** İstek gövdesine eklenecek JSON (yapay zeka motorları için). */
  extraBody?: string;
  /** Ek HTTP başlıkları (JSON). */
  customHeaders?: string;
  /** OpenAI uyumlu hazır ayar kimliği. */
  preset?: string;
  /** Silinemeyen yerleşik motor. */
  builtin?: boolean;
}

export interface Expert {
  id: string;
  name: string;
  nameEn?: string;
  description?: string;
  descriptionEn?: string;
  /** {{to}}, {{from}}, {{title}}, {{summary}} değişkenlerini destekleyen sistem istemi. */
  systemPrompt: string;
  temperature?: number;
  icon?: string;
  group?: 'academic' | 'general';
  builtin?: boolean;
}

export interface GlossaryEntry {
  id: string;
  source: string;
  target: string;
  caseSensitive?: boolean;
  note?: string;
}

export interface SiteRule {
  id: string;
  /** Alan adı veya joker karakterli desen (örn. *.example.com). */
  pattern: string;
  /** Yalnızca bu seçicilerin içi çevrilir (virgülle ayrılmış CSS). */
  include?: string;
  /** Bu seçiciler hiç çevrilmez. */
  exclude?: string;
  enabled: boolean;
}

export type StyleId =
  | 'none'
  | 'underline'
  | 'dashed'
  | 'wavy'
  | 'highlight'
  | 'marker'
  | 'blur'
  | 'faded'
  | 'italic'
  | 'bold'
  | 'border'
  | 'quote'
  | 'card';

export interface TranslationStyle {
  id: StyleId;
  /** Vurgu rengi (alt çizgi, işaretleyici vb.). */
  color: string;
  /** Çeviri metni rengi; boşsa sayfanınkini kullanır. */
  textColor: string;
  /** Yazı boyutu çarpanı (1 = aynı). */
  fontScale: number;
}

export interface Settings {
  version: number;
  uiLanguage: UiLocalePref;
  theme: ThemePref;
  sourceLanguage: string;
  targetLanguage: string;
  engineId: string;
  expertId: string;
  displayMode: DisplayMode;
  style: TranslationStyle;
  rules: {
    alwaysSites: string[];
    neverSites: string[];
    alwaysLanguages: string[];
    neverLanguages: string[];
    custom: SiteRule[];
  };
  page: {
    /** Hedef dilden farklı dildeki tüm sayfaları açılınca otomatik çevir (varsayılan kapalı). */
    autoTranslateAll: boolean;
    scope: TranslationScope;
    /** Yalnızca ekranda görünen paragrafları çevir. */
    lazy: boolean;
    skipCode: boolean;
    /** Kısa metinleri (menü, düğme) satır içinde göster. */
    inlineShort: boolean;
    /** Sayfa başlığı ve özetini yapay zekaya bağlam olarak gönder. */
    contextAware: boolean;
    translateIframes: boolean;
  };
  academic: {
    /** Akademik sitelerde akademik uzmanı otomatik seç. */
    autoExpert: boolean;
    /** Otomatik seçilecek akademik uzman. */
    expertId: string;
    /** Kaynakça bölümlerini çevirme. */
    skipReferences: boolean;
    /** Yerleşik karşılığı olmayan terimlerin aslını parantez içinde koru. */
    keepTerms: boolean;
    /** Uzun belgelerde terimleri tutarlı çevirmek için önce terim listesi çıkar. */
    termConsistency: boolean;
  };
  floatingBall: {
    enabled: boolean;
    side: 'left' | 'right';
    /** Ekran yüksekliğine göre yüzde. */
    top: number;
    hiddenSites: string[];
  };
  hover: {
    key: HoverKey;
  };
  input: {
    trigger: InputTrigger;
    targetLanguage: string;
    /** Metin zaten hedef dildeyse ana dile çevir. */
    smartReverse: boolean;
    /** Boşsa varsayılan motor. */
    engineId: string;
    expertId: string;
  };
  selection: {
    mode: SelectionMode;
    /** Çevirinin gösterileceği yer: seçilen paragrafın altı veya açılır kart. */
    display: SelectionDisplay;
    engineId: string;
    inInputs: boolean;
  };
  engines: Record<string, EngineConfig>;
  /** Kullanıcının eklediği uzmanlar. */
  experts: Expert[];
  /** Yerleşik uzmanlarda yapılan değişiklikler. */
  expertOverrides: Record<string, Partial<Pick<Expert, 'systemPrompt' | 'temperature'>>>;
  glossary: GlossaryEntry[];
  cache: {
    enabled: boolean;
    ttlDays: number;
  };
  pdf: {
    /** Tarayıcıda açılan PDF'leri BioParrot PDF çevirmeninde aç. */
    autoOpen: boolean;
  };
}

/** Yapay zekaya gönderilen sayfa/belge bağlamı. */
export interface TranslationContext {
  title?: string;
  summary?: string;
  url?: string;
  /** Belgeye özel terim listesi (terim tutarlılığı). */
  terms?: { source: string; target: string }[];
}

/** Bir çeviri isteği. */
export interface TranslateRequest {
  texts: string[];
  from: string;
  to: string;
  engineId?: string;
  expertId?: string;
  /** Metinler <a0>…</a0> biçiminde yer tutucu etiket içeriyor. */
  rich?: boolean;
  context?: TranslationContext;
  /** Önbelleği atla. */
  noCache?: boolean;
}

export interface TranslateResponse {
  translations: string[];
  /** Çeviriler zengin biçimde mi (etiketli, HTML kaçışlı) yoksa düz metin mi. */
  rich: boolean;
  /** Motorun algıladığı kaynak dil (metin başına). */
  detected?: (string | undefined)[];
  engineId: string;
}

export class TranslateError extends Error {
  constructor(
    message: string,
    public code: 'auth' | 'quota' | 'network' | 'config' | 'parse' | 'server' | 'refused' | 'unknown' = 'unknown',
    public status?: number,
    /** 429 yanıtındaki Retry-After (ms). */
    public retryAfterMs?: number,
  ) {
    super(message);
    this.name = 'TranslateError';
  }
}

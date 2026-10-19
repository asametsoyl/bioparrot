/**
 * Desteklenen diller. `code` uygulama içi standart koddur; motorlara özel
 * kodlar `google` / `microsoft` / `deepl` alanlarında tutulur (yoksa `code`).
 */
export interface Language {
  code: string;
  /** İngilizce adı (yapay zeka istemlerinde de kullanılır). */
  en: string;
  tr: string;
  native: string;
  google?: string;
  microsoft?: string;
  /** DeepL hedef dil kodu (kaynak kodu bundan türetilir). */
  deepl?: string;
  /** Sağdan sola yazılan dil. */
  rtl?: boolean;
}

export const LANGUAGES: Language[] = [
  { code: 'en', en: 'English', tr: 'İngilizce', native: 'English', deepl: 'EN-US' },
  { code: 'tr', en: 'Turkish', tr: 'Türkçe', native: 'Türkçe' },
  { code: 'de', en: 'German', tr: 'Almanca', native: 'Deutsch' },
  { code: 'fr', en: 'French', tr: 'Fransızca', native: 'Français' },
  { code: 'es', en: 'Spanish', tr: 'İspanyolca', native: 'Español' },
  { code: 'it', en: 'Italian', tr: 'İtalyanca', native: 'Italiano' },
  { code: 'pt', en: 'Brazilian Portuguese', tr: 'Portekizce (Brezilya)', native: 'Português (Brasil)', deepl: 'PT-BR' },
  { code: 'pt-PT', en: 'European Portuguese', tr: 'Portekizce (Portekiz)', native: 'Português (Portugal)', microsoft: 'pt-pt', deepl: 'PT-PT' },
  { code: 'ru', en: 'Russian', tr: 'Rusça', native: 'Русский' },
  { code: 'uk', en: 'Ukrainian', tr: 'Ukraynaca', native: 'Українська' },
  { code: 'ar', en: 'Arabic', tr: 'Arapça', native: 'العربية', rtl: true },
  { code: 'fa', en: 'Persian', tr: 'Farsça', native: 'فارسی', rtl: true },
  { code: 'zh-CN', en: 'Simplified Chinese', tr: 'Çince (Basitleştirilmiş)', native: '简体中文', microsoft: 'zh-Hans', deepl: 'ZH-HANS' },
  { code: 'zh-TW', en: 'Traditional Chinese', tr: 'Çince (Geleneksel)', native: '繁體中文', microsoft: 'zh-Hant', deepl: 'ZH-HANT' },
  { code: 'ja', en: 'Japanese', tr: 'Japonca', native: '日本語' },
  { code: 'ko', en: 'Korean', tr: 'Korece', native: '한국어' },
  { code: 'hi', en: 'Hindi', tr: 'Hintçe', native: 'हिन्दी' },
  { code: 'bn', en: 'Bengali', tr: 'Bengalce', native: 'বাংলা' },
  { code: 'ur', en: 'Urdu', tr: 'Urduca', native: 'اردو', rtl: true },
  { code: 'id', en: 'Indonesian', tr: 'Endonezce', native: 'Bahasa Indonesia' },
  { code: 'ms', en: 'Malay', tr: 'Malayca', native: 'Bahasa Melayu' },
  { code: 'vi', en: 'Vietnamese', tr: 'Vietnamca', native: 'Tiếng Việt' },
  { code: 'th', en: 'Thai', tr: 'Tayca', native: 'ไทย' },
  { code: 'nl', en: 'Dutch', tr: 'Felemenkçe', native: 'Nederlands' },
  { code: 'sv', en: 'Swedish', tr: 'İsveççe', native: 'Svenska' },
  { code: 'no', en: 'Norwegian', tr: 'Norveççe', native: 'Norsk', microsoft: 'nb', deepl: 'NB' },
  { code: 'da', en: 'Danish', tr: 'Danca', native: 'Dansk' },
  { code: 'fi', en: 'Finnish', tr: 'Fince', native: 'Suomi' },
  { code: 'pl', en: 'Polish', tr: 'Lehçe', native: 'Polski' },
  { code: 'cs', en: 'Czech', tr: 'Çekçe', native: 'Čeština' },
  { code: 'sk', en: 'Slovak', tr: 'Slovakça', native: 'Slovenčina' },
  { code: 'hu', en: 'Hungarian', tr: 'Macarca', native: 'Magyar' },
  { code: 'ro', en: 'Romanian', tr: 'Rumence', native: 'Română' },
  { code: 'bg', en: 'Bulgarian', tr: 'Bulgarca', native: 'Български' },
  { code: 'el', en: 'Greek', tr: 'Yunanca', native: 'Ελληνικά' },
  { code: 'he', en: 'Hebrew', tr: 'İbranice', native: 'עברית', google: 'iw', rtl: true },
  { code: 'sr', en: 'Serbian', tr: 'Sırpça', native: 'Српски', microsoft: 'sr-Cyrl' },
  { code: 'hr', en: 'Croatian', tr: 'Hırvatça', native: 'Hrvatski' },
  { code: 'bs', en: 'Bosnian', tr: 'Boşnakça', native: 'Bosanski' },
  { code: 'sl', en: 'Slovenian', tr: 'Slovence', native: 'Slovenščina' },
  { code: 'mk', en: 'Macedonian', tr: 'Makedonca', native: 'Македонски' },
  { code: 'sq', en: 'Albanian', tr: 'Arnavutça', native: 'Shqip' },
  { code: 'be', en: 'Belarusian', tr: 'Belarusça', native: 'Беларуская' },
  { code: 'et', en: 'Estonian', tr: 'Estonca', native: 'Eesti' },
  { code: 'lv', en: 'Latvian', tr: 'Letonca', native: 'Latviešu' },
  { code: 'lt', en: 'Lithuanian', tr: 'Litvanca', native: 'Lietuvių' },
  { code: 'az', en: 'Azerbaijani', tr: 'Azerbaycan Türkçesi', native: 'Azərbaycanca' },
  { code: 'ka', en: 'Georgian', tr: 'Gürcüce', native: 'ქართული' },
  { code: 'hy', en: 'Armenian', tr: 'Ermenice', native: 'Հայերեն' },
  { code: 'kk', en: 'Kazakh', tr: 'Kazakça', native: 'Қазақ тілі' },
  { code: 'uz', en: 'Uzbek', tr: 'Özbekçe', native: 'Oʻzbekcha' },
  { code: 'ky', en: 'Kyrgyz', tr: 'Kırgızca', native: 'Кыргызча' },
  { code: 'tk', en: 'Turkmen', tr: 'Türkmence', native: 'Türkmençe' },
  { code: 'tt', en: 'Tatar', tr: 'Tatarca', native: 'Татарча' },
  { code: 'ba', en: 'Bashkir', tr: 'Başkurtça', native: 'Башҡортса' },
  { code: 'cv', en: 'Chuvash', tr: 'Çuvaşça', native: 'Чӑвашла' },
  { code: 'crh', en: 'Crimean Tatar', tr: 'Kırım Tatarcası', native: 'Qırımtatarca' },
  { code: 'ug', en: 'Uyghur', tr: 'Uygurca', native: 'ئۇيغۇرچە', rtl: true },
  { code: 'ku', en: 'Kurdish (Kurmanji)', tr: 'Kürtçe (Kurmanci)', native: 'Kurdî', microsoft: 'kmr' },
  { code: 'ckb', en: 'Kurdish (Sorani)', tr: 'Kürtçe (Sorani)', native: 'کوردی', microsoft: 'ku', rtl: true },
  { code: 'ps', en: 'Pashto', tr: 'Peştuca', native: 'پښتو', rtl: true },
  { code: 'tg', en: 'Tajik', tr: 'Tacikçe', native: 'Тоҷикӣ' },
  { code: 'mn', en: 'Mongolian', tr: 'Moğolca', native: 'Монгол', microsoft: 'mn-Cyrl' },
  { code: 'ne', en: 'Nepali', tr: 'Nepalce', native: 'नेपाली' },
  { code: 'si', en: 'Sinhala', tr: 'Seylanca', native: 'සිංහල' },
  { code: 'ta', en: 'Tamil', tr: 'Tamilce', native: 'தமிழ்' },
  { code: 'te', en: 'Telugu', tr: 'Telugu dili', native: 'తెలుగు' },
  { code: 'kn', en: 'Kannada', tr: 'Kannada dili', native: 'ಕನ್ನಡ' },
  { code: 'ml', en: 'Malayalam', tr: 'Malayalam dili', native: 'മലയാളം' },
  { code: 'mr', en: 'Marathi', tr: 'Marathi dili', native: 'मराठी' },
  { code: 'gu', en: 'Gujarati', tr: 'Guceratça', native: 'ગુજરાતી' },
  { code: 'pa', en: 'Punjabi', tr: 'Pencapça', native: 'ਪੰਜਾਬੀ' },
  { code: 'or', en: 'Odia', tr: 'Oriya dili', native: 'ଓଡ଼ିଆ' },
  { code: 'as', en: 'Assamese', tr: 'Assamca', native: 'অসমীয়া' },
  { code: 'sd', en: 'Sindhi', tr: 'Sindhi dili', native: 'سنڌي', rtl: true },
  { code: 'dv', en: 'Dhivehi', tr: 'Divehi dili', native: 'ދިވެހި', rtl: true },
  { code: 'bo', en: 'Tibetan', tr: 'Tibetçe', native: 'བོད་ཡིག' },
  { code: 'my', en: 'Burmese', tr: 'Birmanca', native: 'မြန်မာ' },
  { code: 'km', en: 'Khmer', tr: 'Kmerce', native: 'ខ្មែរ' },
  { code: 'lo', en: 'Lao', tr: 'Laoca', native: 'ລາວ' },
  { code: 'fil', en: 'Filipino', tr: 'Filipince', native: 'Filipino', google: 'tl' },
  { code: 'jv', en: 'Javanese', tr: 'Cava dili', native: 'Basa Jawa', google: 'jw' },
  { code: 'su', en: 'Sundanese', tr: 'Sunda dili', native: 'Basa Sunda' },
  { code: 'sw', en: 'Swahili', tr: 'Svahili', native: 'Kiswahili' },
  { code: 'am', en: 'Amharic', tr: 'Amharca', native: 'አማርኛ' },
  { code: 'ti', en: 'Tigrinya', tr: 'Tigrinya dili', native: 'ትግርኛ' },
  { code: 'ha', en: 'Hausa', tr: 'Hausa dili', native: 'Hausa' },
  { code: 'yo', en: 'Yoruba', tr: 'Yorubaca', native: 'Yorùbá' },
  { code: 'ig', en: 'Igbo', tr: 'İbo dili', native: 'Igbo' },
  { code: 'zu', en: 'Zulu', tr: 'Zuluca', native: 'isiZulu' },
  { code: 'xh', en: 'Xhosa', tr: 'Kosa dili', native: 'isiXhosa' },
  { code: 'af', en: 'Afrikaans', tr: 'Afrikaanca', native: 'Afrikaans' },
  { code: 'so', en: 'Somali', tr: 'Somalice', native: 'Soomaali' },
  { code: 'mg', en: 'Malagasy', tr: 'Malgaşça', native: 'Malagasy' },
  { code: 'rw', en: 'Kinyarwanda', tr: 'Kinyarwanda dili', native: 'Ikinyarwanda' },
  { code: 'ny', en: 'Chichewa', tr: 'Çiçeva dili', native: 'Chichewa' },
  { code: 'sn', en: 'Shona', tr: 'Shona dili', native: 'chiShona' },
  { code: 'st', en: 'Sesotho', tr: 'Sesotho dili', native: 'Sesotho' },
  { code: 'is', en: 'Icelandic', tr: 'İzlandaca', native: 'Íslenska' },
  { code: 'ga', en: 'Irish', tr: 'İrlandaca', native: 'Gaeilge' },
  { code: 'cy', en: 'Welsh', tr: 'Galce', native: 'Cymraeg' },
  { code: 'gd', en: 'Scottish Gaelic', tr: 'İskoç Galcesi', native: 'Gàidhlig' },
  { code: 'eu', en: 'Basque', tr: 'Baskça', native: 'Euskara' },
  { code: 'ca', en: 'Catalan', tr: 'Katalanca', native: 'Català' },
  { code: 'gl', en: 'Galician', tr: 'Galiçyaca', native: 'Galego' },
  { code: 'mt', en: 'Maltese', tr: 'Maltaca', native: 'Malti' },
  { code: 'lb', en: 'Luxembourgish', tr: 'Lüksemburgca', native: 'Lëtzebuergesch' },
  { code: 'fy', en: 'Frisian', tr: 'Frizce', native: 'Frysk' },
  { code: 'co', en: 'Corsican', tr: 'Korsikaca', native: 'Corsu' },
  { code: 'la', en: 'Latin', tr: 'Latince', native: 'Latina' },
  { code: 'eo', en: 'Esperanto', tr: 'Esperanto', native: 'Esperanto' },
  { code: 'yi', en: 'Yiddish', tr: 'Yidiş', native: 'ייִדיש', rtl: true },
  { code: 'ht', en: 'Haitian Creole', tr: 'Haiti Kreolü', native: 'Kreyòl Ayisyen' },
  { code: 'haw', en: 'Hawaiian', tr: 'Havai dili', native: 'ʻŌlelo Hawaiʻi' },
  { code: 'mi', en: 'Maori', tr: 'Maori dili', native: 'Te Reo Māori' },
  { code: 'sm', en: 'Samoan', tr: 'Samoa dili', native: 'Gagana Sāmoa' },
];

const BY_CODE = new Map(LANGUAGES.map((l) => [l.code, l]));

/** Farklı kaynaklardan gelen dil kodlarını uygulamanın standart koduna çevirir. */
const ALIASES: Record<string, string> = {
  iw: 'he',
  jw: 'jv',
  tl: 'fil',
  nb: 'no',
  nn: 'no',
  in: 'id',
  mo: 'ro',
  zh: 'zh-CN',
  'zh-cn': 'zh-CN',
  'zh-sg': 'zh-CN',
  'zh-hans': 'zh-CN',
  'zh-tw': 'zh-TW',
  'zh-hk': 'zh-TW',
  'zh-mo': 'zh-TW',
  'zh-hant': 'zh-TW',
  'pt-br': 'pt',
  'pt-pt': 'pt-PT',
  'sr-cyrl': 'sr',
  'sr-latn': 'sr',
  'mn-cyrl': 'mn',
  kmr: 'ku',
};

/** Herhangi bir dil kodunu (en-US, zh-Hant, iw…) standart koda çevirir. Bilinmiyorsa temel kodu döner. */
export function normalizeLang(code: string | undefined | null): string {
  if (!code) return '';
  const raw = code.trim().replace('_', '-');
  const lower = raw.toLowerCase();
  if (lower === 'auto' || lower === 'und') return lower === 'auto' ? 'auto' : '';
  if (ALIASES[lower]) return ALIASES[lower];
  if (BY_CODE.has(raw)) return raw;
  const base = lower.split('-')[0];
  if (ALIASES[base]) return ALIASES[base];
  return base;
}

export function getLanguage(code: string): Language | undefined {
  return BY_CODE.get(normalizeLang(code));
}

/** İki dil kodunun aynı dili gösterip göstermediği (Çince yazı sistemleri ayrı sayılır). */
export function isSameLanguage(a: string | undefined, b: string | undefined): boolean {
  const x = normalizeLang(a);
  const y = normalizeLang(b);
  if (!x || !y) return false;
  if (x === y) return true;
  if (x.startsWith('zh') || y.startsWith('zh')) return false;
  return x.split('-')[0] === y.split('-')[0];
}

export function languageName(code: string, ui: 'tr' | 'en'): string {
  if (code === 'auto') return ui === 'tr' ? 'Otomatik algıla' : 'Auto detect';
  const lang = getLanguage(code);
  if (!lang) return code;
  return ui === 'tr' ? lang.tr : lang.en;
}

/** Yapay zeka istemlerinde kullanılan İngilizce ad. */
export function englishName(code: string): string {
  return getLanguage(code)?.en ?? code;
}

export function isRtl(code: string): boolean {
  return !!getLanguage(code)?.rtl;
}

export function toGoogleCode(code: string): string {
  if (code === 'auto') return 'auto';
  const lang = getLanguage(code);
  return lang?.google ?? lang?.code ?? code;
}

export function toMicrosoftCode(code: string): string {
  const lang = getLanguage(code);
  return lang?.microsoft ?? lang?.code ?? code;
}

export function toDeeplTarget(code: string): string {
  const lang = getLanguage(code);
  return (lang?.deepl ?? lang?.code ?? code).toUpperCase();
}

export function toDeeplSource(code: string): string | undefined {
  if (code === 'auto') return undefined;
  return toDeeplTarget(code).split('-')[0];
}

/** Motorların döndürdüğü algılanan dil kodunu standartlaştırır. */
export function fromEngineCode(code: string | undefined, engine: 'google' | 'microsoft' | 'deepl' | 'other'): string {
  if (!code) return '';
  if (engine === 'microsoft') {
    const lower = code.toLowerCase();
    if (lower === 'ku') return 'ckb';
    if (lower === 'kmr') return 'ku';
  }
  return normalizeLang(code);
}

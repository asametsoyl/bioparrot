/**
 * Akademik makalelerin standart bölüm başlıkları için kesin karşılıklar.
 * Klasik çeviri motorları bu tek kelimelik başlıkları bağlamsız çevirince hata yapar
 * ("Abstract" → "Soyut"); bu sözlük başlıkları doğrudan doğru karşılıkla gösterir.
 */
const HEADINGS: Record<string, Record<string, string>> = {
  tr: {
    abstract: 'Özet',
    summary: 'Özet',
    'graphical abstract': 'Grafik Özet',
    highlights: 'Öne Çıkanlar',
    keywords: 'Anahtar Kelimeler',
    'key words': 'Anahtar Kelimeler',
    introduction: 'Giriş',
    background: 'Arka Plan',
    'related work': 'İlgili Çalışmalar',
    'literature review': 'Literatür Taraması',
    methods: 'Yöntemler',
    method: 'Yöntem',
    methodology: 'Yöntem',
    'materials and methods': 'Gereç ve Yöntemler',
    'material and methods': 'Gereç ve Yöntemler',
    'methods and materials': 'Yöntem ve Gereçler',
    'experimental procedures': 'Deneysel Yöntemler',
    'experimental section': 'Deneysel Bölüm',
    experiments: 'Deneyler',
    'study design': 'Çalışma Tasarımı',
    participants: 'Katılımcılar',
    patients: 'Hastalar',
    'statistical analysis': 'İstatistiksel Analiz',
    'statistical analyses': 'İstatistiksel Analizler',
    results: 'Bulgular',
    findings: 'Bulgular',
    'results and discussion': 'Bulgular ve Tartışma',
    discussion: 'Tartışma',
    limitations: 'Sınırlılıklar',
    'strengths and limitations': 'Güçlü Yönler ve Sınırlılıklar',
    conclusion: 'Sonuç',
    conclusions: 'Sonuçlar',
    'conclusion and future work': 'Sonuç ve Gelecek Çalışmalar',
    'future work': 'Gelecek Çalışmalar',
    'case report': 'Olgu Sunumu',
    'case presentation': 'Olgu Sunumu',
    acknowledgments: 'Teşekkür',
    acknowledgements: 'Teşekkür',
    acknowledgment: 'Teşekkür',
    acknowledgement: 'Teşekkür',
    funding: 'Finansman',
    'author contributions': 'Yazar Katkıları',
    'competing interests': 'Çıkar Çatışması',
    'conflict of interest': 'Çıkar Çatışması',
    'conflicts of interest': 'Çıkar Çatışması',
    'declaration of competing interest': 'Çıkar Çatışması Beyanı',
    'data availability': 'Veri Erişilebilirliği',
    'data availability statement': 'Veri Erişilebilirliği Beyanı',
    'code availability': 'Kod Erişilebilirliği',
    'ethics statement': 'Etik Beyanı',
    'ethics approval': 'Etik Onay',
    references: 'Kaynaklar',
    bibliography: 'Kaynakça',
    'supplementary material': 'Ek Materyal',
    'supplementary materials': 'Ek Materyaller',
    'supplementary information': 'Ek Bilgiler',
    'supplementary data': 'Ek Veriler',
    appendix: 'Ek',
    appendices: 'Ekler',
    figures: 'Şekiller',
    tables: 'Tablolar',
    abbreviations: 'Kısaltmalar',
  },
};

const NUMBERING = /^((?:\d+(?:\.\d+)*|[IVXLC]+|[A-Z])[.)]?\s+)/;

/** Bölüm başlığı sözlükte varsa hedef dildeki karşılığını döner (numaralandırma ve büyük harf korunur). */
export function lookupSectionHeading(text: string, target: string): string | undefined {
  const dict = HEADINGS[target.split('-')[0]];
  if (!dict) return undefined;
  const trimmed = text.replace(/\s+/g, ' ').trim();
  if (trimmed.length > 60) return undefined;
  const num = NUMBERING.exec(trimmed)?.[1] ?? '';
  const body = trimmed.slice(num.length).replace(/[:.]$/, '').trim();
  const hit = dict[body.toLowerCase()];
  if (!hit) return undefined;
  const upper = body.length > 2 && body === body.toUpperCase();
  return num + (upper ? hit.toLocaleUpperCase(target) : hit);
}

/** Kaynak büyük harfle başlıyorsa çevirinin ilk harfini de büyütür (Türkçe i/İ dahil). */
export function matchInitialCase(source: string, translation: string, target: string): string {
  const s = source.trimStart();
  const firstSrc = s.match(/\p{L}/u)?.[0];
  if (!firstSrc || firstSrc !== firstSrc.toLocaleUpperCase() || firstSrc === firstSrc.toLocaleLowerCase()) return translation;
  return translation.replace(/\p{L}/u, (ch) => ch.toLocaleUpperCase(target));
}

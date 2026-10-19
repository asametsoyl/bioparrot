/** Uygulamanın görünen adı. Yeniden adlandırmak için burayı, public/_locales ve zotero/addon dosyalarını değiştirin. */
export const APP_NAME = 'BioParrot';

/** Sayfaya eklenen tüm öğelerde kullanılan CSS/öznitelik öneki. */
export const PREFIX = 'bp';

export const CLS = {
  wrap: `${PREFIX}-wrap`,
  block: `${PREFIX}-block`,
  inline: `${PREFIX}-inline`,
  text: `${PREFIX}-text`,
  loading: `${PREFIX}-loading`,
  error: `${PREFIX}-error`,
  only: `${PREFIX}-only`,
  ui: `${PREFIX}-ui`,
  hoverTarget: `${PREFIX}-hover-target`,
  origHidden: `${PREFIX}-orig-hidden`,
} as const;

export const ATTR = {
  unit: `data-${PREFIX}-unit`,
} as const;

export const STORAGE_KEYS = {
  settings: 'settings',
  usage: 'usage',
} as const;

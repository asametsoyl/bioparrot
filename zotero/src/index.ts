/**
 * BioParrot Zotero eklentisi ana betiği. bootstrap.js tarafından yüklenir ve
 * Zotero.BioParrot nesnesini oluşturur.
 */
import { setLocale, t } from '@/core/i18n';
import { installGlobals } from './env';
import { openBilingualPdf, translateMetadata } from './items';
import { initPrefsPane } from './prefs-ui';
import { getPref } from './prefs';
import { onTextSelectionPopup, translationFor } from './reader';
import { maybeRunSelfTest } from './selftest';
import { getService, translateText } from './service';

declare const Zotero: any;
declare const _globalThis: any;
declare const rootURI: string;
declare const pluginID: string;

installGlobals(typeof _globalThis !== 'undefined' ? _globalThis : globalThis);

const MENU_IDS = ['bioparrot-itemmenu-sep', 'bioparrot-itemmenu-translate', 'bioparrot-itemmenu-pdf'];
let notifierID: string | null = null;

function selectedItems(win: any): any[] {
  try {
    return win.ZoteroPane?.getSelectedItems?.() ?? [];
  } catch {
    return [];
  }
}

function addMenus(win: any): void {
  const doc = win.document;
  const menu = doc.getElementById('zotero-itemmenu');
  if (!menu || doc.getElementById(MENU_IDS[1])) return;
  const icon = rootURI + 'content/icons/favicon@0.5x.png';
  const sep = doc.createXULElement('menuseparator');
  sep.id = MENU_IDS[0];
  const translate = doc.createXULElement('menuitem');
  translate.id = MENU_IDS[1];
  translate.setAttribute('label', t('zot.menuTranslateMeta'));
  translate.classList.add('menuitem-iconic');
  translate.setAttribute('image', icon);
  translate.addEventListener('command', () => void translateMetadata(selectedItems(win)));
  const pdf = doc.createXULElement('menuitem');
  pdf.id = MENU_IDS[2];
  pdf.setAttribute('label', t('zot.menuBilingualPdf'));
  pdf.classList.add('menuitem-iconic');
  pdf.setAttribute('image', icon);
  pdf.addEventListener('command', () => {
    const item = selectedItems(win)[0];
    if (item) void openBilingualPdf(win, item, rootURI);
  });
  menu.append(sep, translate, pdf);
  menu.addEventListener('popupshowing', () => {
    const items = selectedItems(win);
    const regular = items.some((i) => i?.isRegularItem?.());
    const pdfCapable = items.length === 1 && (items[0]?.isPDFAttachment?.() || items[0]?.isRegularItem?.());
    translate.hidden = !regular;
    pdf.hidden = !pdfCapable;
    sep.hidden = translate.hidden && pdf.hidden;
  });
}

function removeMenus(win: any): void {
  for (const id of MENU_IDS) win.document.getElementById(id)?.remove();
}

const hooks = {
  async onStartup() {
    setLocale(String(Zotero.locale ?? '').toLowerCase().startsWith('tr') ? 'tr' : 'en');
    const win = Zotero.getMainWindow?.();
    getService(win?.fetch ? win.fetch.bind(win) : undefined);

    try {
      Zotero.PreferencePanes.register({
        pluginID,
        src: rootURI + 'content/preferences.xhtml',
        scripts: [rootURI + 'content/preferences.js'],
        label: 'BioParrot',
        image: rootURI + 'content/icons/favicon.png',
      });
    } catch (err) {
      Zotero.logError(err);
    }

    Zotero.Reader.registerEventListener('renderTextSelectionPopup', onTextSelectionPopup, pluginID);

    // Vurgu oluşturulunca, aynı metnin çevirisini yorum olarak ekle (isteğe bağlı).
    notifierID = Zotero.Notifier.registerObserver(
      {
        notify: async (event: string, type: string, ids: number[]) => {
          if (event !== 'add' || type !== 'item' || !getPref<boolean>('annotationComment', false)) return;
          for (const id of ids) {
            const item = Zotero.Items.get(id);
            if (!item?.isAnnotation?.() || item.annotationComment || !item.annotationText) continue;
            const translation = translationFor(item.annotationText);
            if (!translation) continue;
            item.annotationComment = translation;
            await item.saveTx();
          }
        },
      },
      ['item'],
      'bioparrot',
    );

    for (const w of Zotero.getMainWindows?.() ?? []) hooks.onMainWindowLoad(w);
    // Yalnızca geliştirici öz testi tercihi ayarlıysa çalışır.
    void maybeRunSelfTest(rootURI);
  },

  onMainWindowLoad(win: any) {
    try {
      addMenus(win);
    } catch (err) {
      Zotero.logError(err);
    }
  },

  onMainWindowUnload(win: any) {
    removeMenus(win);
  },

  onShutdown() {
    try {
      Zotero.Reader.unregisterEventListener('renderTextSelectionPopup', onTextSelectionPopup);
    } catch {
      /* zaten kaldırılmış */
    }
    if (notifierID) Zotero.Notifier.unregisterObserver(notifierID);
    for (const w of Zotero.getMainWindows?.() ?? []) removeMenus(w);
  },
};

Zotero.BioParrot = {
  hooks,
  version: '0.1.0',
  initPrefsPane,
  // Komut satırından/testlerden kullanım için
  translateText,
  translateMetadata,
  openBilingualPdf: (win: any, item: any) => openBilingualPdf(win, item, rootURI),
  onTextSelectionPopup,
};

/**
 * Geliştirici öz testi: yalnızca `extensions.bioparrot.selfTest` tercihi bir PDF yoluna ayarlandığında çalışır.
 * Sonuçları `extensions.bioparrot.selfTestOut` dosyasına JSON olarak yazar ve Zotero'yu kapatır.
 * Normal kullanımda hiçbir etkisi yoktur.
 */
import { getPref } from './prefs';
import { openBilingualPdf, translateMetadata } from './items';
import { onTextSelectionPopup } from './reader';
import { translateText } from './service';

declare const Zotero: any;
declare const Services: any;
declare const IOUtils: any;

const sleep = (ms: number) => new Promise((r) => setTimeout(r, ms));

export async function maybeRunSelfTest(rootURI: string): Promise<void> {
  const pdfPath = getPref<string>('selfTest', '');
  const outPath = getPref<string>('selfTestOut', '');
  if (!pdfPath || !outPath) return;
  const results: { step: string; ok: boolean; detail: unknown }[] = [];
  const step = async (name: string, fn: () => Promise<unknown>) => {
    try {
      const detail = await fn();
      results.push({ step: name, ok: true, detail });
    } catch (err) {
      results.push({ step: name, ok: false, detail: String((err as Error)?.stack ?? err) });
    }
    await IOUtils.writeUTF8(outPath, JSON.stringify(results, null, 2));
  };

  await sleep(3000);
  await step('plugin', async () => ({ zotero: Zotero.version, locale: Zotero.locale }));
  await step('menus', async () => {
    const win = Zotero.getMainWindow();
    return ['bioparrot-itemmenu-translate', 'bioparrot-itemmenu-pdf'].map((id) => !!win?.document.getElementById(id));
  });
  await step('translateText', async () => translateText('The cell membrane is a selectively permeable barrier.'));

  let parentID = 0;
  await step('metadataNote', async () => {
    const item = new Zotero.Item('journalArticle');
    item.setField('title', 'CRISPR-Cas9 editing of TP53 reveals a role in organoid growth');
    item.setField('abstractNote', 'Mitochondria are the powerhouse of the cell. Here we show that loss of TP53 in mouse organoids increases growth.');
    await item.saveTx();
    parentID = item.id;
    await translateMetadata([item]);
    return item.getNotes().map((id: number) => Zotero.Items.get(id).getNote().replace(/<[^>]+>/g, ' ').replace(/\s+/g, ' ').trim());
  });

  await step('importPdf', async () => {
    const att = await Zotero.Attachments.importFromFile({ file: pdfPath, parentItemID: parentID });
    return { id: att.id, file: att.attachmentFilename };
  });

  await step('readerPopup', async () => {
    const att = Zotero.Items.get(Zotero.Items.get(parentID).getAttachments()[0]);
    const reader = await Zotero.Reader.open(att.id);
    for (let i = 0; i < 60 && !reader?._iframeWindow?.document?.body; i++) await sleep(500);
    const doc = reader._iframeWindow.document;
    let box: any = null;
    onTextSelectionPopup({
      reader,
      doc,
      params: { annotation: { text: 'The dominant sequence trans-\nduction models are based on complex recurrent networks.' } },
      append: (el: Element) => {
        box = el;
        doc.body.appendChild(el);
      },
    });
    for (let i = 0; i < 40 && box && /Çevriliyor|Translating/.test(box.textContent); i++) await sleep(500);
    return box?.textContent?.slice(0, 300);
  });

  await step('bilingualPdf', async () => {
    const parent = Zotero.Items.get(parentID);
    const before = parent.getAttachments().length;
    const att = Zotero.Items.get(parent.getAttachments()[0]);
    const win = Zotero.getMainWindow();
    win.alert = () => {};
    await openBilingualPdf(win, att, rootURI);
    let dlg: any = null;
    for (let i = 0; i < 60 && !dlg; i++) {
      await sleep(500);
      dlg = [...Services.wm.getEnumerator(null)].find((w: any) => String(w.location?.href).includes('bioparrot/content/pdf.html'));
    }
    if (!dlg) throw new Error('PDF penceresi açılmadı');
    dlg.alert = () => {};
    let status = '';
    for (let i = 0; i < 300; i++) {
      await sleep(500);
      const doc = dlg.document;
      status = doc.querySelector('.pdf-status')?.textContent ?? doc.querySelector('.pdf-error')?.textContent ?? doc.body?.textContent?.slice(0, 160) ?? '';
      if (/tamamlandı|complete|açılamadı|could not/i.test(status)) break;
    }
    const blocks = dlg.document.querySelectorAll('.pv-block').length;
    dlg.document.querySelector('.export .btn-primary')?.click();
    await sleep(400);
    dlg.document.querySelector('.export .menu button')?.click();
    let after = before;
    for (let i = 0; i < 180 && after === before; i++) {
      await sleep(500);
      after = Zotero.Items.get(parentID).getAttachments().length;
    }
    const ids = Zotero.Items.get(parentID).getAttachments();
    const newest = Zotero.Items.get(ids[ids.length - 1]);
    const error = dlg.document.querySelector('.pdf-error')?.textContent;
    dlg.close();
    return { status: status.trim().slice(0, 120), blocks, before, after, newest: newest?.attachmentFilename, error };
  });

  await step('prefsPane', async () => (Zotero.PreferencePanes.pluginPanes ?? []).map((p: any) => p.label ?? p.id));
  await sleep(500);
  Services.startup.quit(Services.startup.eForceQuit);
}

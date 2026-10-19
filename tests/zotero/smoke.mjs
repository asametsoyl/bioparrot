// Zotero içinde BioParrot eklentisinin duman testi (geçici profil, -marionette ile başlatılmış Zotero gerekir).
import { Marionette } from './marionette.mjs';

const m = new Marionette(2829);
const hello = await m.connect(90000);
console.log('Marionette:', JSON.stringify(hello));
await m.send('WebDriver:NewSession', { capabilities: {} });
await m.send('Marionette:SetContext', { value: 'chrome' });

const step = async (name, fn, args = [], timeout = 120000) => {
  try {
    const v = await m.exec(fn, args, timeout);
    console.log(`✓ ${name}:`, typeof v === 'string' ? v : JSON.stringify(v));
    return v;
  } catch (err) {
    console.log(`✗ ${name}:`, String(err.message).split('\n')[0]);
    return null;
  }
};

await step('Zotero hazır ve eklenti yüklendi', async () => {
  for (let i = 0; i < 120 && !(Zotero.initialized && Zotero.BioParrot); i++) await new Promise((r) => setTimeout(r, 500));
  return { version: Zotero.version, plugin: !!Zotero.BioParrot, pluginVersion: Zotero.BioParrot?.version, locale: Zotero.locale };
});

await step('Menüler eklendi', async () => {
  const win = Zotero.getMainWindow();
  return ['bioparrot-itemmenu-translate', 'bioparrot-itemmenu-pdf'].map((id) => !!win.document.getElementById(id));
});

await step('Eklenti içinden çeviri (Google)', async () => {
  return await Zotero.BioParrot.translateText('The cell membrane is a selectively permeable barrier.');
});

const itemID = await step('Öğe oluşturma + başlık/özet çevirisi notu', async () => {
  const item = new Zotero.Item('journalArticle');
  item.setField('title', 'CRISPR-Cas9 editing of TP53 reveals a role in organoid growth');
  item.setField('abstractNote', 'Mitochondria are the powerhouse of the cell. Here we show that loss of TP53 in mouse organoids increases growth.');
  await item.saveTx();
  await Zotero.BioParrot.translateMetadata([item]);
  const notes = item.getNotes().map((id) => Zotero.Items.get(id).getNote());
  return { id: item.id, notes: notes.map((n) => n.replace(/<[^>]+>/g, ' ').replace(/\s+/g, ' ').slice(0, 220)) };
});

if (process.argv[2]) {
  await step('PDF eki içe aktarma', async (pdfPath, parentID) => {
    const att = await Zotero.Attachments.importFromFile({ file: pdfPath, parentItemID: parentID });
    return att.id;
  }, [process.argv[2], itemID?.id]);

  await step('Okuyucu seçim penceresi (doğrudan işleyici)', async (parentID) => {
    const parent = Zotero.Items.get(parentID);
    const att = Zotero.Items.get(parent.getAttachments()[0]);
    const reader = await Zotero.Reader.open(att.id);
    for (let i = 0; i < 60 && !reader?._iframeWindow?.document?.body; i++) await new Promise((r) => setTimeout(r, 500));
    const doc = reader._iframeWindow.document;
    let box = null;
    Zotero.BioParrot.onTextSelectionPopup({ reader, doc, params: { annotation: { text: 'The dominant sequence trans-\nduction models are based on complex recurrent networks.' } }, append: (el) => { box = el; doc.body.appendChild(el); } });
    for (let i = 0; i < 40 && box && /Çevriliyor|Translating/.test(box.textContent); i++) await new Promise((r) => setTimeout(r, 500));
    return box?.textContent?.slice(0, 200);
  }, [itemID?.id]);

  await step('İki dilli PDF penceresi + dışa aktarma', async (parentID) => {
    const parent = Zotero.Items.get(parentID);
    const before = parent.getAttachments().length;
    const att = Zotero.Items.get(parent.getAttachments()[0]);
    const win = Zotero.getMainWindow();
    win.alert = () => {};
    await Zotero.BioParrot.openBilingualPdf(win, att);
    let dlg = null;
    for (let i = 0; i < 60 && !dlg; i++) {
      await new Promise((r) => setTimeout(r, 500));
      dlg = [...Services.wm.getEnumerator(null)].find((w) => String(w.location?.href).includes('bioparrot/content/pdf.html'));
    }
    if (!dlg) return 'pencere açılmadı';
    dlg.alert = () => {};
    let status = '';
    for (let i = 0; i < 240; i++) {
      await new Promise((r) => setTimeout(r, 500));
      status = dlg.document.querySelector('.pdf-status')?.textContent ?? dlg.document.querySelector('.pdf-error')?.textContent ?? dlg.document.body?.textContent?.slice(0, 120) ?? '';
      if (/tamamlandı|complete|hata|error|açılamadı/i.test(status)) break;
    }
    const blocks = dlg.document.querySelectorAll('.pv-block').length;
    dlg.document.querySelector('.export .btn-primary')?.click();
    await new Promise((r) => setTimeout(r, 300));
    dlg.document.querySelector('.export .menu button')?.click();
    let after = before;
    for (let i = 0; i < 120 && after === before; i++) {
      await new Promise((r) => setTimeout(r, 500));
      after = Zotero.Items.get(parentID).getAttachments().length;
    }
    const newest = Zotero.Items.get(Zotero.Items.get(parentID).getAttachments().at(-1));
    dlg.close();
    return { status: status.trim().slice(0, 80), blocks, attachmentsBefore: before, attachmentsAfter: after, newest: newest?.attachmentFilename };
  }, [itemID?.id], 240000);
}

await step('Tercih paneli kayıtlı', async () => {
  return Zotero.PreferencePanes.pluginPanes?.map?.((p) => p.label ?? p.id) ?? 'bilinmiyor';
});

m.close();

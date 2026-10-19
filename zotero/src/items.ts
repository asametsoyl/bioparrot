import { languageName } from '@/core/languages';
import { escapeHtml } from '@/core/markup';
import { getLocale, t } from '@/core/i18n';
import { getService } from './service';
import { loadSettings } from './prefs';

declare const Zotero: any;

function progress(text: string): any {
  const win = new Zotero.ProgressWindow({ closeOnClick: true });
  win.changeHeadline('BioParrot');
  win.addDescription(text);
  win.show();
  return win;
}

/** Seçili öğelerin başlık ve özetini çevirip alt not olarak ekler. */
export async function translateMetadata(items: any[]): Promise<void> {
  const regular = items.filter((i) => i?.isRegularItem?.());
  if (!regular.length) return;
  const pw = progress(t('zot.translatingItems'));
  const s = loadSettings();
  let done = 0;
  try {
    for (const item of regular) {
      const title: string = item.getField('title') ?? '';
      const abstract: string = item.getField('abstractNote') ?? '';
      const texts = [title, abstract].filter((x) => x.trim());
      if (!texts.length) continue;
      const res = await getService().translate({
        texts,
        from: s.sourceLanguage,
        to: s.targetLanguage,
        engineId: s.engineId,
        expertId: s.expertId,
        rich: false,
        context: { title, summary: abstract },
      });
      const trTitle = title ? res.translations[0] : '';
      const trAbstract = abstract ? res.translations[title ? 1 : 0] : '';
      const lang = languageName(s.targetLanguage, getLocale());
      const html = [
        `<h1>${escapeHtml(t('zot.noteTitle', { lang }))}</h1>`,
        trTitle ? `<p><strong>${escapeHtml(t('zot.noteTitleLabel'))}:</strong> ${escapeHtml(trTitle)}</p>` : '',
        trAbstract ? `<p><strong>${escapeHtml(t('zot.noteAbstractLabel'))}:</strong></p><p>${escapeHtml(trAbstract).replace(/\n+/g, '</p><p>')}</p>` : '',
        `<p><em>${escapeHtml(s.engines[s.engineId]?.name ?? s.engineId)}</em></p>`,
      ].join('');
      const note = new Zotero.Item('note');
      note.libraryID = item.libraryID;
      note.parentID = item.id;
      note.setNote(html);
      await note.saveTx();
      done++;
    }
    pw.addDescription(t('zot.itemsDone', { n: done }));
  } catch (err) {
    pw.addDescription(t('err.generic', { msg: String((err as Error)?.message ?? err) }));
  }
  pw.startCloseTimer(5000);
}

/** Seçili öğenin PDF ekini iki dilli PDF penceresinde açar. */
export async function openBilingualPdf(win: any, item: any, rootURI: string): Promise<void> {
  let attachment = item;
  if (!item?.isAttachment?.() || !item.isPDFAttachment?.()) attachment = await item?.getBestAttachment?.();
  if (!attachment || !attachment.isPDFAttachment?.()) {
    win.alert(t('zot.noPdf'));
    return;
  }
  const path = await attachment.getFilePathAsync();
  if (!path) {
    win.alert(t('zot.fileMissing'));
    return;
  }
  win.openDialog(
    'chrome://bioparrot/content/pdf.html',
    `bioparrot-pdf-${attachment.id}`,
    'chrome,dialog=no,resizable,centerscreen,width=1400,height=920',
    {
      Zotero,
      rootURI,
      attachmentID: attachment.id,
      parentID: attachment.parentID ?? null,
      path,
      name: (attachment.attachmentFilename as string) || `${attachment.getField('title') || 'document'}.pdf`,
    },
  );
}

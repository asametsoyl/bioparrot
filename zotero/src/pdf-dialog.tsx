/**
 * Zotero içinde açılan iki dilli PDF penceresi. Tarayıcı eklentisindeki PDF çevirmeninin aynısını
 * kullanır; çeviriler doğrudan Zotero tercihlerindeki motorla yapılır ve sonuç öğeye ek olarak kaydedilir.
 */
import { useCallback, useState } from 'react';
import { createRoot } from 'react-dom/client';
import { setFetch } from '@/core/http';
import { setLocale, t } from '@/core/i18n';
import { TranslationService } from '@/core/service';
import type { Settings } from '@/core/types';
import { PdfApp, type PdfAppHost, type PdfFile } from '@/pdf-app/PdfApp';
import type { UpdateSettings } from '@/ui/hooks';
import '@/ui/theme.css';
import '@/pdf-app/pdf-app.css';
import { loadSettings, saveBasicSettings } from './prefs';

declare const IOUtils: { read(path: string): Promise<Uint8Array>; write(path: string, data: Uint8Array): Promise<number>; remove(path: string): Promise<void> };
declare const PathUtils: { join(...parts: string[]): string; tempDir: string };

interface DialogArgs {
  Zotero: any;
  rootURI: string;
  attachmentID: number;
  parentID: number | null;
  path: string;
  name: string;
}

const args: DialogArgs = (window as any).arguments?.[0] ?? (window as any).opener?.__bioparrotPdfArgs ?? {};
const Zotero = args.Zotero ?? (window as any).opener?.Zotero;
(globalThis as any).Zotero = Zotero;
setLocale(String(Zotero?.locale ?? navigator.language).toLowerCase().startsWith('tr') ? 'tr' : 'en');
setFetch((input, init) => window.fetch(input, init));

const service = new TranslationService({ getSettings: () => loadSettings() });
const BASE = 'chrome://bioparrot/content/pdf/';

const host: PdfAppHost = {
  workerSrc: BASE + 'pdf.worker.min.mjs',
  fontUrls: { regular: BASE + 'fonts/NotoSans-Regular.ttf', bold: BASE + 'fonts/NotoSans-Bold.ttf', math: BASE + 'fonts/NotoSansMath-Regular.ttf' },
  translate: (req) => service.translate(req),
  extractTerms: (text, to, engineId, expertId) => service.extractTerms(text, to, engineId, expertId),
  async loadUrl(url: string): Promise<PdfFile> {
    const res = await fetch(url);
    return { name: url.split('/').pop() || 'document.pdf', data: new Uint8Array(await res.arrayBuffer()) };
  },
  async save(bytes, filename) {
    const tmp = PathUtils.join(PathUtils.tempDir, filename.replace(/[\\/:*?"<>|]/g, '_'));
    await IOUtils.write(tmp, bytes);
    try {
      const attachment = Zotero.Items.get(args.attachmentID);
      await Zotero.Attachments.importFromFile({
        file: tmp,
        parentItemID: args.parentID ?? undefined,
        libraryID: args.parentID ? undefined : attachment?.libraryID,
        title: filename,
      });
      window.alert(t('zot.saved', { name: filename }));
    } finally {
      await IOUtils.remove(tmp).catch(() => {});
    }
  },
  openSettings: () => {
    try {
      Zotero.Utilities.Internal.openPreferences('bioparrot@bioparrot.app');
    } catch {
      /* eski sürüm */
    }
  },
};

function App({ file }: { file: PdfFile }) {
  const [settings, setSettings] = useState<Settings>(() => loadSettings());
  const update = useCallback<UpdateSettings>(async (mutator) => {
    setSettings((prev) => {
      const draft = structuredClone(prev);
      mutator(draft);
      saveBasicSettings(draft);
      return draft;
    });
  }, []);
  return <PdfApp host={host} settings={settings} update={update} initialFile={file} />;
}

async function main() {
  const root = createRoot(document.getElementById('root')!);
  try {
    const data = await IOUtils.read(args.path);
    document.title = `${args.name} – BioParrot`;
    root.render(<App file={{ name: args.name, data }} />);
  } catch (err) {
    root.render(<p style={{ padding: 24 }}>{t('pdf.loadFailed', { msg: String((err as Error)?.message ?? err) })}</p>);
  }
}

void main();

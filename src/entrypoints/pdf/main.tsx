import { createRoot } from 'react-dom/client';
import { browser } from 'wxt/browser';
import pdfWorkerUrl from 'pdfjs-dist/build/pdf.worker.min.mjs?url';
import notoRegularUrl from '@expo-google-fonts/noto-sans/400Regular/NotoSans_400Regular.ttf?url';
import notoBoldUrl from '@expo-google-fonts/noto-sans/700Bold/NotoSans_700Bold.ttf?url';
import notoMathUrl from '@expo-google-fonts/noto-sans-math/400Regular/NotoSansMath_400Regular.ttf?url';
import { useEffect, useMemo } from 'react';
import { sendToBackground } from '@/lib/messaging';
import { SelectionPopup } from '@/content/selection';
import { PdfApp, type PdfAppHost, type PdfFile } from '@/pdf-app/PdfApp';
import { useSettings, useTheme } from '@/ui/hooks';
import '@/ui/theme.css';
import '@/pdf-app/pdf-app.css';

function fileNameFromUrl(url: string, res?: Response): string {
  const cd = res?.headers.get('content-disposition') ?? '';
  const m = /filename\*?=(?:UTF-8'')?"?([^";]+)"?/i.exec(cd);
  if (m) return decodeURIComponent(m[1]);
  try {
    const path = new URL(url).pathname.split('/').filter(Boolean).pop() ?? 'document.pdf';
    const name = decodeURIComponent(path);
    return /\.pdf$/i.test(name) ? name : `${name}.pdf`;
  } catch {
    return 'document.pdf';
  }
}

const host: PdfAppHost = {
  workerSrc: pdfWorkerUrl,
  fontUrls: { regular: notoRegularUrl, bold: notoBoldUrl, math: notoMathUrl },
  translate: (req) => sendToBackground({ type: 'translate', req }),
  extractTerms: (text, to, engineId, expertId) => sendToBackground({ type: 'extractTerms', text, to, engineId, expertId }),
  async loadUrl(url: string): Promise<PdfFile> {
    const res = await fetch(url, { credentials: 'include' });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const data = new Uint8Array(await res.arrayBuffer());
    if (!(data[0] === 0x25 && data[1] === 0x50 && data[2] === 0x44 && data[3] === 0x46)) {
      throw new Error('Bu adres bir PDF dosyası döndürmedi (giriş gerektiren bir sayfa olabilir).');
    }
    return { name: fileNameFromUrl(url, res), data };
  },
  async save(bytes, filename) {
    const blob = new Blob([bytes as BlobPart], { type: 'application/pdf' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 30_000);
  },
  openSettings: () => void browser.tabs.create({ url: browser.runtime.getURL('/options.html#documents') }),
};

function App() {
  const { settings, update } = useSettings();
  useTheme(settings?.theme);
  const initialUrl = useMemo(() => new URLSearchParams(location.search).get('file'), []);

  // PDF sayfasında da metin seçince çeviri kartı açılsın.
  useEffect(() => {
    if (!settings) return;
    const popup = new SelectionPopup(() => settings, { title: document.title, summary: '', isPaper: true, isPdf: true });
    popup.attach();
    return () => popup.detach();
  }, [settings]);

  if (!settings) return null;
  return <PdfApp host={host} settings={settings} update={update} initialUrl={initialUrl} />;
}

createRoot(document.getElementById('root')!).render(<App />);

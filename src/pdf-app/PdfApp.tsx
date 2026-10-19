import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import { isLlmEngine } from '@/core/engines/meta';
import { t } from '@/core/i18n';
import { buildTranslatedPdf, type ExportMode } from '@/core/pdf/export';
import { PdfTranslationJob, type JobState, type PdfJobHost } from '@/core/pdf/job';
import { configurePdfWorker, openPdf, type PDFDocumentProxy } from '@/core/pdf/pdfjs';
import type { Settings } from '@/core/types';
import { EngineOptions, ExpertOptions, LanguageSelect, Logo, Segmented, Select } from '@/ui/components';
import type { UpdateSettings } from '@/ui/hooks';
import { PageView, type ViewMode } from './PageView';

export interface PdfFile {
  name: string;
  data: Uint8Array;
}

export interface PdfAppHost extends PdfJobHost {
  workerSrc: string;
  fontUrls: { regular: string; bold: string; math: string };
  loadUrl(url: string): Promise<PdfFile>;
  save(bytes: Uint8Array, filename: string): Promise<void>;
  openSettings?(): void;
}

interface Props {
  host: PdfAppHost;
  settings: Settings;
  update: UpdateSettings;
  initialUrl?: string | null;
  initialFile?: PdfFile | null;
}

interface Loaded {
  file: PdfFile;
  doc: PDFDocumentProxy;
  sizes: { width: number; height: number }[];
}

async function fetchBytes(url: string): Promise<Uint8Array> {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  return new Uint8Array(await res.arrayBuffer());
}

export function PdfApp({ host, settings, update, initialUrl, initialFile }: Props) {
  const [loaded, setLoaded] = useState<Loaded | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState('');
  const [job, setJob] = useState<PdfTranslationJob | null>(null);
  const [jobState, setJobState] = useState<JobState | null>(null);
  const [version, setVersion] = useState(0);
  const [view, setView] = useState<ViewMode>('dual');
  const [zoom, setZoom] = useState(1);
  const [containerWidth, setContainerWidth] = useState(1200);
  const [exportMenu, setExportMenu] = useState(false);
  const [exporting, setExporting] = useState<string | null>(null);
  const [url, setUrl] = useState('');
  const [dragging, setDragging] = useState(false);
  const mainRef = useRef<HTMLDivElement>(null);
  const fileInput = useRef<HTMLInputElement>(null);
  const jobKey = useRef('');

  useEffect(() => configurePdfWorker(host.workerSrc), [host.workerSrc]);

  useEffect(() => {
    const el = mainRef.current;
    if (!el) return;
    const ro = new ResizeObserver(() => setContainerWidth(el.clientWidth));
    ro.observe(el);
    setContainerWidth(el.clientWidth);
    return () => ro.disconnect();
  }, [loaded]);

  const load = useCallback(async (file: PdfFile) => {
    setLoading(true);
    setError('');
    try {
      const doc = await openPdf(file.data);
      const sizes = [];
      for (let i = 1; i <= doc.numPages; i++) {
        const page = await doc.getPage(i);
        const vp = page.getViewport({ scale: 1 });
        sizes.push({ width: vp.width, height: vp.height });
      }
      setLoaded({ file, doc, sizes });
      document.title = `${file.name} – BioParrot`;
    } catch (err) {
      setError(t('pdf.loadFailed', { msg: String((err as Error)?.message ?? err) }));
    } finally {
      setLoading(false);
    }
  }, []);

  useEffect(() => {
    if (initialFile) void load(initialFile);
    else if (initialUrl) {
      setLoading(true);
      host
        .loadUrl(initialUrl)
        .then(load)
        .catch((err) => {
          setLoading(false);
          const msg = String((err as Error)?.message ?? err);
          setError(t('pdf.loadFailed', { msg }) + (initialUrl.startsWith('file:') ? `\n${t('pdf.fileAccess')}` : ''));
        });
    }
  }, [initialFile, initialUrl, host, load]);

  // Belge veya çeviri ayarları değişince işi (yeniden) başlat.
  const settingsKey = `${settings.targetLanguage}|${settings.sourceLanguage}|${settings.engineId}|${settings.expertId}|${settings.academic.skipReferences}`;
  useEffect(() => {
    if (!loaded) return;
    const key = `${loaded.file.name}|${loaded.doc.numPages}|${settingsKey}`;
    if (jobKey.current === key) return;
    jobKey.current = key;
    job?.stop();
    const next = new PdfTranslationJob(loaded.doc, host, settings, (j) => {
      setJobState({ ...j.state });
      setVersion((v) => v + 1);
    });
    setJob(next);
    setJobState({ ...next.state });
    void next.run();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [loaded, settingsKey]);

  const maxWidth = useMemo(() => Math.max(...(loaded?.sizes.map((s) => s.width) ?? [600])), [loaded]);
  const fitScale = Math.max(0.3, (containerWidth - 48 - (view === 'dual' ? 24 : 0)) / (maxWidth * (view === 'dual' ? 2 : 1)));
  const scale = Math.min(fitScale * zoom, 4);

  const onFiles = async (files: FileList | null) => {
    const f = files?.[0];
    if (!f) return;
    await load({ name: f.name, data: new Uint8Array(await f.arrayBuffer()) });
  };

  const doExport = async (mode: ExportMode) => {
    if (!loaded || !job) return;
    setExportMenu(false);
    setExporting(t('pdf.exporting', { done: 0, total: loaded.doc.numPages }));
    try {
      const [regular, bold, math] = await Promise.all([fetchBytes(host.fontUrls.regular), fetchBytes(host.fontUrls.bold), fetchBytes(host.fontUrls.math)]);
      const bytes = await buildTranslatedPdf({
        original: loaded.file.data,
        layouts: job.layouts,
        translations: job.translations,
        rich: jobState?.rich ?? true,
        mode,
        fonts: { regular, bold, math },
        credit: `BioParrot (${settings.engines[settings.engineId]?.name ?? ''} → ${settings.targetLanguage})`,
        onProgress: (done, total) => setExporting(t('pdf.exporting', { done, total })),
      });
      const base = loaded.file.name.replace(/\.pdf$/i, '');
      const suffix = mode === 'dual' ? 'bilingual' : mode === 'interleaved' ? 'interleaved' : settings.targetLanguage;
      await host.save(bytes, `${base}.${suffix}.pdf`);
    } catch (err) {
      setError(t('pdf.exportFailed', { msg: String((err as Error)?.message ?? err) }));
    } finally {
      setExporting(null);
    }
  };

  const llm = isLlmEngine(settings.engines[settings.engineId]);
  const st = jobState;
  const running = st?.phase === 'analyzing' || st?.phase === 'translating';
  const failedPages = st?.pageStatus.filter((p) => p === 'error').length ?? 0;
  const progress = st ? (st.phase === 'analyzing' ? st.analyzed / Math.max(st.total, 1) : st.pagesDone / Math.max(st.total, 1)) : 0;
  const statusText = !st
    ? ''
    : st.phase === 'analyzing'
      ? t('pdf.analyzing', { done: st.analyzed, total: st.total })
      : st.phase === 'translating'
        ? t('pdf.translating', { done: st.pagesDone, total: st.total })
        : st.phase === 'done'
          ? t('pdf.done')
          : st.phase === 'stopped'
            ? t('pdf.stopped')
            : st.error ?? '';

  return (
    <div
      className={`pdf-app ${dragging ? 'dragging' : ''}`}
      onDragOver={(e) => {
        e.preventDefault();
        setDragging(true);
      }}
      onDragLeave={() => setDragging(false)}
      onDrop={(e) => {
        e.preventDefault();
        setDragging(false);
        void onFiles(e.dataTransfer.files);
      }}
    >
      <header className="pdf-bar">
        <Logo size={26} withText={!loaded} sub={t('pdf.title')} />
        {loaded && (
          <div className="pdf-file" title={loaded.file.name}>
            {loaded.file.name}
            <span className="faint"> · {loaded.doc.numPages} s.</span>
          </div>
        )}
        <div className="pdf-controls">
          <LanguageSelect value={settings.targetLanguage} onChange={(v) => update((d) => void (d.targetLanguage = v))} ariaLabel={t('lang.target')} />
          <Select value={settings.engineId} onChange={(v) => update((d) => void (d.engineId = v))} ariaLabel={t('engine.label')}>
            <EngineOptions settings={settings} />
          </Select>
          {llm && (
            <Select
              value={settings.academic.autoExpert && settings.expertId === 'general' ? settings.academic.expertId : settings.expertId}
              onChange={(v) => update((d) => void (d.expertId = v))}
              ariaLabel={t('expert.label')}
            >
              <ExpertOptions settings={settings} />
            </Select>
          )}
        </div>
        <div className="spacer" />
        {loaded && (
          <>
            <Segmented
              value={view}
              onChange={setView}
              options={[
                { value: 'dual', label: t('pdf.view.dual') },
                { value: 'translated', label: t('pdf.view.translated') },
                { value: 'original', label: t('pdf.view.original') },
              ]}
            />
            <div className="zoom">
              <button className="btn btn-ghost btn-icon" onClick={() => setZoom((z) => Math.max(0.4, z / 1.15))} aria-label="zoom out">
                −
              </button>
              <span className="small muted">{Math.round(zoom * 100)}%</span>
              <button className="btn btn-ghost btn-icon" onClick={() => setZoom((z) => Math.min(3, z * 1.15))} aria-label="zoom in">
                +
              </button>
            </div>
            {running ? (
              <button className="btn" onClick={() => job?.stop()}>
                {t('pdf.stop')}
              </button>
            ) : (
              <button className="btn" onClick={() => void job?.run()}>
                {st?.phase === 'done' ? t('pdf.restart') : t('pdf.start')}
              </button>
            )}
            <div className="export">
              <button className="btn btn-primary" disabled={!!exporting || !job?.translations.size} onClick={() => setExportMenu((v) => !v)}>
                ⬇ {t('pdf.download')}
              </button>
              {exportMenu && (
                <div className="menu">
                  <button onClick={() => void doExport('dual')}>{t('pdf.export.dual')}</button>
                  <button onClick={() => void doExport('translated')}>{t('pdf.export.translated')}</button>
                  <button onClick={() => void doExport('interleaved')}>{t('pdf.export.interleaved')}</button>
                </div>
              )}
            </div>
          </>
        )}
        {host.openSettings && (
          <button className="btn btn-ghost btn-icon" onClick={host.openSettings} title={t('common.settings')}>
            ⚙
          </button>
        )}
      </header>

      {loaded && st && (
        <div className="pdf-status">
          <div className="bar">
            <div className="fill" style={{ width: `${Math.round(progress * 100)}%` }} />
          </div>
          <span className="small">
            {running && <span className="spinner" />} {statusText}
            {job?.context.terms?.length ? ` · ${t('pdf.terms', { n: job.context.terms.length })}` : ''}
            {exporting ? ` · ${exporting}` : ''}
          </span>
          {failedPages > 0 && !running && <span className="small err-text">{t('pdf.pageErrors', { n: failedPages })}</span>}
        </div>
      )}

      <main className="pdf-main" ref={mainRef}>
        {error && <div className="pdf-error">{error}</div>}
        {!loaded && !loading && (
          <div className="dropzone">
            <div className="dz-icon">📄</div>
            <h2>{t('pdf.drop')}</h2>
            <p className="muted">{t('pdf.dropHint')}</p>
            <button className="btn btn-primary" onClick={() => fileInput.current?.click()}>
              {t('pdf.choose')}
            </button>
            <div className="input-group" style={{ width: 'min(520px, 90%)' }}>
              <input className="input" value={url} placeholder={t('pdf.urlPlaceholder')} onChange={(e) => setUrl(e.target.value)} />
              <button
                className="btn"
                disabled={!url.trim()}
                onClick={() => {
                  setLoading(true);
                  host
                    .loadUrl(url.trim())
                    .then(load)
                    .catch((err) => {
                      setLoading(false);
                      setError(t('pdf.loadFailed', { msg: String((err as Error)?.message ?? err) }));
                    });
                }}
              >
                {t('pdf.open')}
              </button>
            </div>
          </div>
        )}
        {loading && (
          <div className="dropzone">
            <span className="spinner" /> {t('pdf.loading')}
          </div>
        )}
        {loaded && (
          <div className="pv-list">
            {loaded.sizes.map((size, i) => (
              <PageView
                key={`${loaded.file.name}-${i}`}
                doc={loaded.doc}
                index={i}
                width={size.width}
                height={size.height}
                scale={scale}
                view={view}
                layout={job?.layouts[i]}
                translations={job?.translations ?? new Map()}
                rich={st?.rich ?? true}
                status={st?.pageStatus[i] ?? 'pending'}
                version={version}
              />
            ))}
          </div>
        )}
        <input ref={fileInput} type="file" accept="application/pdf,.pdf" hidden onChange={(e) => void onFiles(e.target.files)} />
      </main>
    </div>
  );
}

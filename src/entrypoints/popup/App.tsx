import { useCallback, useEffect, useState } from 'react';
import { browser, type Browser } from 'wxt/browser';
import { isLlmEngine } from '@/core/engines/meta';
import { getAllExperts } from '@/core/experts';
import { getLocale, t } from '@/core/i18n';
import { languageName } from '@/core/languages';
import { sendTabAction, type PageState } from '@/lib/messaging';
import { siteMatches } from '@/lib/settings';
import { EngineOptions, ExpertOptions, LanguageSelect, Logo, Segmented, Select, Switch } from '@/ui/components';
import { useSettings, useTheme } from '@/ui/hooks';

type PageStatus = 'loading' | 'unsupported' | PageState;

function hostOf(url?: string): string {
  try {
    return url ? new URL(url).hostname : '';
  } catch {
    return '';
  }
}

function looksLikePdf(url?: string): boolean {
  return !!url && /\.pdf($|[?#])|arxiv\.org\/pdf\//i.test(url);
}

export default function App() {
  const { settings, update } = useSettings();
  useTheme(settings?.theme);
  const [tab, setTab] = useState<Browser.tabs.Tab | null>(null);
  const [status, setStatus] = useState<PageStatus>('loading');
  const [shortcut, setShortcut] = useState('Alt+A');

  const refresh = useCallback(async (target: Browser.tabs.Tab | null, inject = false) => {
    if (!target?.id) return setStatus('unsupported');
    let state = await sendTabAction<PageState>(target.id, { type: 'getState' }, true);
    if (!state && inject && /^(https?|file):/.test(target.url ?? '')) {
      try {
        await browser.scripting.executeScript({ target: { tabId: target.id, allFrames: true }, files: ['/content-scripts/content.js'] });
        await new Promise((r) => setTimeout(r, 200));
        state = await sendTabAction<PageState>(target.id, { type: 'getState' }, true);
      } catch {
        /* korumalı sayfa */
      }
    }
    setStatus(state ?? 'unsupported');
  }, []);

  useEffect(() => {
    void (async () => {
      const [active] = await browser.tabs.query({ active: true, currentWindow: true });
      setTab(active ?? null);
      await refresh(active ?? null, true);
      const commands = await browser.commands.getAll().catch(() => []);
      const cmd = commands.find((c) => c.name === 'toggle-translate');
      setShortcut(cmd?.shortcut ?? '');
    })();
  }, [refresh]);

  // Çeviri sürerken durumu düzenli yenile.
  useEffect(() => {
    if (typeof status !== 'object') return;
    const timer = setInterval(() => void refresh(tab), 900);
    return () => clearInterval(timer);
  }, [status, tab, refresh]);

  if (!settings) return <div className="popup popup-loading" />;

  const ui = getLocale();
  const host = hostOf(tab?.url);
  const page = typeof status === 'object' ? status : null;
  const unsupported = status === 'unsupported';
  const isPdf = page?.isPdf || looksLikePdf(tab?.url);
  const engine = settings.engines[settings.engineId];
  const llm = isLlmEngine(engine);
  const autoExpert = page?.expertId && page.expertId !== settings.expertId ? getAllExperts(settings).find((e) => e.id === page.expertId) : undefined;

  const openOptions = (section?: string) => {
    void browser.tabs.create({ url: browser.runtime.getURL('/options.html') + (section ? `#${section}` : '') });
    window.close();
  };
  const openPdf = (file?: string) => {
    void browser.tabs.create({ url: browser.runtime.getURL('/pdf.html') + (file ? `?file=${encodeURIComponent(file)}` : '') });
    window.close();
  };

  const onTranslate = async () => {
    if (isPdf) return openPdf(tab?.url);
    if (!tab?.id || unsupported) return;
    await sendTabAction(tab.id, { type: 'toggle' }, true);
    setTimeout(() => void refresh(tab), 250);
  };

  const toggleList = (key: 'alwaysSites' | 'neverSites', on: boolean) =>
    update((d) => {
      const list = d.rules[key].filter((s) => !siteMatches(host, [s]));
      if (on) list.push(host);
      d.rules[key] = list;
      const other = key === 'alwaysSites' ? 'neverSites' : 'alwaysSites';
      if (on) d.rules[other] = d.rules[other].filter((s) => !siteMatches(host, [s]));
    });

  const pageLang = page?.pageLang;
  const showLangToggle = !!pageLang && pageLang !== settings.targetLanguage;

  return (
    <div className="popup">
      <header className="pop-head">
        <Logo size={28} sub={t('app.tagline')} />
        <div className="spacer" />
        <button className="btn btn-ghost btn-icon" title={t('common.settings')} onClick={() => openOptions()}>
          <svg viewBox="0 0 24 24" width="18" height="18" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
            <circle cx="12" cy="12" r="3" />
            <path d="M19.4 15a1.7 1.7 0 0 0 .3 1.8l.1.1a2 2 0 1 1-2.8 2.8l-.1-.1a1.7 1.7 0 0 0-1.8-.3 1.7 1.7 0 0 0-1 1.5V21a2 2 0 1 1-4 0v-.1a1.7 1.7 0 0 0-1.1-1.5 1.7 1.7 0 0 0-1.8.3l-.1.1a2 2 0 1 1-2.8-2.8l.1-.1a1.7 1.7 0 0 0 .3-1.8 1.7 1.7 0 0 0-1.5-1H3a2 2 0 1 1 0-4h.1a1.7 1.7 0 0 0 1.5-1.1 1.7 1.7 0 0 0-.3-1.8l-.1-.1a2 2 0 1 1 2.8-2.8l.1.1a1.7 1.7 0 0 0 1.8.3H9a1.7 1.7 0 0 0 1-1.5V3a2 2 0 1 1 4 0v.1a1.7 1.7 0 0 0 1 1.5 1.7 1.7 0 0 0 1.8-.3l.1-.1a2 2 0 1 1 2.8 2.8l-.1.1a1.7 1.7 0 0 0-.3 1.8V9a1.7 1.7 0 0 0 1.5 1H21a2 2 0 1 1 0 4h-.1a1.7 1.7 0 0 0-1.5 1z" />
          </svg>
        </button>
      </header>

      <main className="pop-body">
        <div className="lang-row">
          <div className="lang-col">
            <label>{t('lang.source')}</label>
            <LanguageSelect includeAuto value={settings.sourceLanguage} onChange={(v) => update((d) => void (d.sourceLanguage = v))} ariaLabel={t('lang.source')} />
          </div>
          <button
            className="swap"
            title={t('lang.swap')}
            disabled={settings.sourceLanguage === 'auto'}
            onClick={() =>
              update((d) => {
                const src = d.sourceLanguage;
                d.sourceLanguage = d.targetLanguage;
                d.targetLanguage = src;
              })
            }
          >
            ⇄
          </button>
          <div className="lang-col">
            <label>{t('lang.target')}</label>
            <LanguageSelect value={settings.targetLanguage} onChange={(v) => update((d) => void (d.targetLanguage = v))} ariaLabel={t('lang.target')} />
          </div>
        </div>

        <div className="field">
          <label>{t('engine.label')}</label>
          <Select value={settings.engineId} onChange={(v) => update((d) => void (d.engineId = v))} ariaLabel={t('engine.label')}>
            <EngineOptions settings={settings} />
          </Select>
        </div>

        {llm ? (
          <div className="field">
            <label>{t('expert.label')}</label>
            <Select value={settings.expertId} onChange={(v) => update((d) => void (d.expertId = v))} ariaLabel={t('expert.label')}>
              <ExpertOptions settings={settings} />
            </Select>
            {(page?.isPaper || autoExpert) && (
              <div className="paper-hint">
                🎓 {t('pop.paperDetected')}
                {autoExpert ? ` (${autoExpert.icon ?? ''} ${ui === 'en' && autoExpert.nameEn ? autoExpert.nameEn : autoExpert.name})` : ''}
              </div>
            )}
          </div>
        ) : (
          <button className="llm-hint" onClick={() => openOptions('engines')}>
            ✨ {t('pop.setupLlm')} →
          </button>
        )}

        <button className={`translate-btn ${page?.active ? 'active' : ''}`} disabled={unsupported && !isPdf} onClick={onTranslate}>
          {page && page.pending > 0 && <span className="spinner" />}
          <span>{isPdf ? t('pop.openPdf') : page?.active ? t('common.showOriginal') : t('pop.translatePage')}</span>
          {shortcut && !isPdf && <kbd>{shortcut}</kbd>}
        </button>
        {unsupported && !isPdf && <p className="unsupported">{t('pop.unsupported')}</p>}

        <div className="mode-row">
          <span className="muted small">{t('mode.label')}</span>
          <Segmented
            value={settings.displayMode}
            onChange={(v) => update((d) => void (d.displayMode = v))}
            options={[
              { value: 'bilingual', label: t('mode.bilingual') },
              { value: 'translation', label: t('mode.translation') },
            ]}
          />
        </div>

        {host && !unsupported && (
          <div className="toggles">
            <label className="toggle">
              <span>{t('pop.alwaysSite')}</span>
              <Switch checked={siteMatches(host, settings.rules.alwaysSites)} onChange={(v) => toggleList('alwaysSites', v)} />
            </label>
            <label className="toggle">
              <span>{t('pop.neverSite')}</span>
              <Switch checked={siteMatches(host, settings.rules.neverSites)} onChange={(v) => toggleList('neverSites', v)} />
            </label>
            {showLangToggle && (
              <label className="toggle">
                <span>{t('pop.alwaysLang', { lang: languageName(pageLang!, ui) })}</span>
                <Switch
                  checked={settings.rules.alwaysLanguages.includes(pageLang!)}
                  onChange={(v) =>
                    update((d) => {
                      d.rules.alwaysLanguages = d.rules.alwaysLanguages.filter((l) => l !== pageLang);
                      if (v) d.rules.alwaysLanguages.push(pageLang!);
                    })
                  }
                />
              </label>
            )}
            <label className="toggle">
              <span>{t('pop.autoAll')}</span>
              <Switch checked={settings.page.autoTranslateAll} onChange={(v) => update((d) => void (d.page.autoTranslateAll = v))} />
            </label>
            <label className="toggle">
              <span>{t('pop.floatingBall')}</span>
              <Switch
                checked={settings.floatingBall.enabled && !siteMatches(host, settings.floatingBall.hiddenSites)}
                onChange={(v) =>
                  update((d) => {
                    d.floatingBall.hiddenSites = d.floatingBall.hiddenSites.filter((s) => !siteMatches(host, [s]));
                    if (v) d.floatingBall.enabled = true;
                    else d.floatingBall.hiddenSites.push(host);
                  })
                }
              />
            </label>
          </div>
        )}
      </main>

      <footer className="pop-tools">
        <button className="tool" onClick={() => openPdf(isPdf ? tab?.url : undefined)}>
          <span className="tool-ico">📄</span>
          <span>{t('pop.pdf')}</span>
        </button>
        <button className="tool" onClick={() => openOptions('glossary')}>
          <span className="tool-ico">📚</span>
          <span>{t('opt.nav.glossary')}</span>
        </button>
        <button className="tool" onClick={() => openOptions('academic')}>
          <span className="tool-ico">🎓</span>
          <span>{t('opt.nav.academic')}</span>
        </button>
        <button className="tool" onClick={() => openOptions('zotero')}>
          <span className="tool-ico">🦜</span>
          <span>{t('pop.zotero')}</span>
        </button>
      </footer>
    </div>
  );
}

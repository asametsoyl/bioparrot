import { browser } from 'wxt/browser';
import { defineContentScript } from 'wxt/utils/define-content-script';
import { setLocale } from '@/core/i18n';
import { isSameLanguage, normalizeLang } from '@/core/languages';
import type { Settings } from '@/core/types';
import { FloatingBall } from '@/content/floating-ball';
import { HoverTranslator } from '@/content/hover';
import { InputTranslator } from '@/content/input-translate';
import { PageTranslator } from '@/content/page-translator';
import { getPageInfo } from '@/content/page-info';
import { SelectionPopup } from '@/content/selection';
import { showErrorToast } from '@/content/toast';
import { sendToBackground, type PageState, type TabAction, type TabActionMessage } from '@/lib/messaging';
import { getSettings, siteMatches, updateSettings, watchSettings } from '@/lib/settings';

declare global {
  interface Window {
    __bioparrot?: boolean;
  }
}

export default defineContentScript({
  matches: ['<all_urls>'],
  allFrames: true,
  matchAboutBlank: true,
  runAt: 'document_idle',
  async main(ctx) {
    // scripting.executeScript ile ikinci kez eklenirse tekrar başlatma.
    if (window.__bioparrot) return;
    window.__bioparrot = true;

    let settings: Settings = await getSettings();
    setLocale(settings.uiLanguage);
    const isTop = window.top === window;
    const pageInfo = getPageInfo();
    const host = location.hostname;
    if (!isTop && !settings.page.translateIframes) return;

    let ball: FloatingBall | null = null;
    const translator = new PageTranslator({
      settings: () => settings,
      pageInfo,
      onChange: () => updateBall(),
      onError: (err) => {
        if (!isTop) return;
        showErrorToast(err.message, err.code === 'auth' || err.code === 'config' || err.code === 'quota');
      },
    });

    // Testler için: yalıtılmış dünyada (sayfa betikleri erişemez) iç duruma erişim.
    (window as unknown as { __bioparrotDebug?: unknown }).__bioparrotDebug = { translator, pageInfo };

    const hover = new HoverTranslator(translator, () => settings);
    hover.attach();
    const input = new InputTranslator(() => settings);
    input.attach();
    const selection = isTop ? new SelectionPopup(() => settings, pageInfo, (range, text) => void translator.translateSelection(range, text)) : null;
    selection?.attach();

    function updateBall() {
      ball?.update({ active: translator.active, pending: translator.pending, isPdf: pageInfo.isPdf });
    }

    function shouldShowBall(s: Settings): boolean {
      return isTop && s.floatingBall.enabled && !siteMatches(host, s.floatingBall.hiddenSites) && document.contentType !== 'text/plain';
    }

    async function requestToggle() {
      if (pageInfo.isPdf) {
        await sendToBackground({ type: 'openPdf', url: location.href }).catch(() => {});
        return;
      }
      const next = !translator.active;
      await sendToBackground({ type: 'broadcast', action: { type: next ? 'translate' : 'restore' } }).catch(() => {
        // Arka planla bağlantı yoksa en azından bu çerçevede uygula.
        if (next) void translator.start();
        else translator.stop();
      });
    }

    function syncBall(s: Settings) {
      if (shouldShowBall(s)) {
        if (!ball) {
          ball = new FloatingBall(
            {
              toggle: () => void requestToggle(),
              toggleMode: () =>
                void updateSettings((d) => {
                  d.displayMode = d.displayMode === 'bilingual' ? 'translation' : 'bilingual';
                }),
              settings: () => void sendToBackground({ type: 'openOptions' }).catch(() => {}),
              hideSite: () =>
                void updateSettings((d) => {
                  if (!d.floatingBall.hiddenSites.includes(host)) d.floatingBall.hiddenSites.push(host);
                }),
              openPdf: () => void sendToBackground({ type: 'openPdf', url: location.href }).catch(() => {}),
              savePosition: (top) =>
                void updateSettings((d) => {
                  d.floatingBall.top = top;
                }),
            },
            s.floatingBall.side,
            s.floatingBall.top,
          );
          ball.mount();
          updateBall();
        } else {
          ball.setPosition(s.floatingBall.side, s.floatingBall.top);
          ball.applyLabels();
        }
      } else if (ball) {
        ball.unmount();
        ball = null;
      }
    }
    syncBall(settings);

    const unwatch = watchSettings((next) => {
      const prev = settings;
      settings = next;
      setLocale(next.uiLanguage);
      translator.onSettingsChanged(prev, next);
      syncBall(next);
      updateBall();
    });

    async function handleAction(action: TabAction): Promise<unknown> {
      switch (action.type) {
        case 'toggle':
          if (isTop) await requestToggle();
          return true;
        case 'translate':
          if (!pageInfo.isPdf) await translator.start();
          return true;
        case 'restore':
          translator.stop();
          return true;
        case 'toggleMode':
          if (isTop)
            await updateSettings((d) => {
              d.displayMode = d.displayMode === 'bilingual' ? 'translation' : 'bilingual';
            });
          return true;
        case 'setMode':
          if (isTop)
            await updateSettings((d) => {
              d.displayMode = action.mode;
            });
          return true;
        case 'translateSelection':
          selection?.openForSelection(action.text);
          return true;
        case 'toggleHover':
          hover.enabled = !hover.enabled;
          return true;
        case 'getState': {
          const state: PageState = {
            active: translator.active,
            mode: settings.displayMode,
            pageLang: await detectPageLanguage(),
            isPaper: pageInfo.isPaper,
            isPdf: pageInfo.isPdf,
            pending: translator.pending,
            expertId: translator.resolveExpert(settings),
            url: location.href,
          };
          return state;
        }
        case 'settingsHint':
          return true;
      }
    }

    const onMessage = (raw: unknown, _sender: unknown, sendResponse: (r: unknown) => void) => {
      const msg = raw as TabActionMessage;
      if (!msg || msg.type !== 'bp:tab') return;
      handleAction(msg.action).then(sendResponse, () => sendResponse(false));
      return true;
    };
    browser.runtime.onMessage.addListener(onMessage);

    let pageLangPromise: Promise<string> | null = null;
    function detectPageLanguage(): Promise<string> {
      pageLangPromise ??= (async () => {
        const sample = (document.body?.innerText ?? '').slice(0, 2000);
        if (sample.trim().length > 40) {
          try {
            const r = await browser.i18n.detectLanguage(sample);
            const top = r?.languages?.[0];
            if (top?.language && (r.isReliable || (top.percentage ?? 0) > 60)) return normalizeLang(top.language);
          } catch {
            /* yok */
          }
        }
        return normalizeLang(pageInfo.htmlLang);
      })();
      return pageLangPromise;
    }

    // Ayar açıksa tarayıcıda açılan PDF'ler BioParrot PDF çevirmenine yönlendirilir.
    if (isTop && pageInfo.isPdf && settings.pdf.autoOpen) {
      location.replace(browser.runtime.getURL('/pdf.html') + `?file=${encodeURIComponent(location.href)}`);
      return;
    }

    // Otomatik çeviri kararı yalnızca ana çerçevede verilir, tüm çerçevelere yayılır.
    if (isTop && !pageInfo.isPdf) {
      const s = settings;
      if (!siteMatches(host, s.rules.neverSites)) {
        let auto = siteMatches(host, s.rules.alwaysSites);
        if (!auto && (s.rules.alwaysLanguages.length || s.page.autoTranslateAll)) {
          const lang = await detectPageLanguage();
          const foreign = !!lang && !isSameLanguage(lang, s.targetLanguage) && !s.rules.neverLanguages.some((l) => isSameLanguage(l, lang));
          auto = foreign && (s.page.autoTranslateAll || s.rules.alwaysLanguages.some((l) => isSameLanguage(l, lang)));
        }
        if (auto) void sendToBackground({ type: 'broadcast', action: { type: 'translate' } }).catch(() => translator.start());
      }
    }

    ctx.onInvalidated(() => {
      unwatch();
      browser.runtime.onMessage.removeListener(onMessage);
      translator.stop();
      hover.detach();
      input.detach();
      selection?.detach();
      ball?.unmount();
      window.__bioparrot = false;
    });
  },
});

import { browser, type Browser } from 'wxt/browser';
import { defineBackground } from 'wxt/utils/define-background';
import { setLocale, t } from '@/core/i18n';
import { TranslationService } from '@/core/service';
import { IdbCache } from '@/lib/idb-cache';
import { errorResult, okResult, sendTabAction, type BackgroundMessage } from '@/lib/messaging';
import { getSettings, watchSettings } from '@/lib/settings';
import { StorageUsageRecorder } from '@/lib/usage';

const MENU = {
  page: 'bp-toggle-page',
  selection: 'bp-translate-selection',
  pdfLink: 'bp-open-pdf-link',
  pdfPage: 'bp-open-pdf-page',
} as const;

export default defineBackground(() => {
  const cache = new IdbCache();
  const usage = new StorageUsageRecorder();
  const service = new TranslationService({ getSettings, cache, usage });

  // ---- Mesajlar -------------------------------------------------------------
  browser.runtime.onMessage.addListener((raw: unknown, sender: Browser.runtime.MessageSender, sendResponse: (r: unknown) => void) => {
    const msg = raw as BackgroundMessage;
    if (!msg || typeof msg !== 'object' || !('type' in msg)) return;
    if ((msg as { type: string }).type === 'bp:tab') return;
    handle(msg, sender)
      .then((data) => sendResponse(okResult(data)))
      .catch((err) => sendResponse(errorResult(err)));
    return true;
  });

  async function handle(msg: BackgroundMessage, sender: Browser.runtime.MessageSender): Promise<unknown> {
    switch (msg.type) {
      case 'translate':
        return service.translate(msg.req);
      case 'extractTerms':
        return service.extractTerms(msg.text, msg.to, msg.engineId, msg.expertId);
      case 'explain':
        return service.explain(msg.text, msg.to, msg.context, msg.engineId, msg.expertId);
      case 'testEngine': {
        const start = Date.now();
        const res = await service.translate({ texts: [msg.text], from: 'auto', to: msg.to, engineId: msg.engineId, noCache: true });
        return { translation: res.translations[0], ms: Date.now() - start };
      }
      case 'openOptions': {
        const url = browser.runtime.getURL('/options.html') + (msg.section ? `#${msg.section}` : '');
        await browser.tabs.create({ url });
        return;
      }
      case 'openPdf': {
        const target = msg.url ?? sender.tab?.url;
        const url = browser.runtime.getURL('/pdf.html') + (target ? `?file=${encodeURIComponent(target)}` : '');
        await browser.tabs.create({ url, index: sender.tab ? sender.tab.index + 1 : undefined });
        return;
      }
      case 'broadcast': {
        const tabId = sender.tab?.id;
        if (tabId !== undefined) await sendTabAction(tabId, msg.action);
        return;
      }
      case 'cacheStats':
        return { entries: await cache.count() };
      case 'clearCache':
        await cache.clear();
        return;
    }
  }

  // ---- Sağ tık menüleri -------------------------------------------------------
  async function buildMenus() {
    const settings = await getSettings();
    setLocale(settings.uiLanguage);
    await browser.contextMenus.removeAll();
    browser.contextMenus.create({ id: MENU.page, title: `${t('common.translate')} / ${t('common.showOriginal')}`, contexts: ['page'] });
    browser.contextMenus.create({ id: MENU.selection, title: `${t('common.translate')}: “%s”`, contexts: ['selection'] });
    browser.contextMenus.create({
      id: MENU.pdfLink,
      title: t('pop.openPdf'),
      contexts: ['link'],
      targetUrlPatterns: ['*://*/*.pdf', '*://*/*.pdf?*', '*://*/*.PDF', '*://arxiv.org/pdf/*', '*://*/*/pdf/*'],
    });
    browser.contextMenus.create({
      id: MENU.pdfPage,
      title: t('pop.openPdf'),
      contexts: ['page'],
      documentUrlPatterns: ['*://*/*.pdf', '*://*/*.pdf?*', '*://*/*.PDF', 'file:///*.pdf', '*://arxiv.org/pdf/*'],
    });
  }

  browser.runtime.onInstalled.addListener(async (details) => {
    await buildMenus();
    const settings = await getSettings();
    await cache.prune(settings.cache.ttlDays).catch(() => {});
    if (details.reason === 'install') {
      await browser.tabs.create({ url: browser.runtime.getURL('/options.html') + '#welcome' });
    }
  });
  browser.runtime.onStartup?.addListener(async () => {
    await buildMenus();
    const settings = await getSettings();
    await cache.prune(settings.cache.ttlDays).catch(() => {});
  });
  watchSettings((s) => {
    setLocale(s.uiLanguage);
    void buildMenus();
  });

  browser.contextMenus.onClicked.addListener(async (info, tab) => {
    const tabId = tab?.id;
    if (tabId === undefined) return;
    switch (info.menuItemId) {
      case MENU.page:
        await ensureContentScript(tabId);
        await sendTabAction(tabId, { type: 'toggle' });
        break;
      case MENU.selection:
        await ensureContentScript(tabId);
        await sendTabAction(tabId, { type: 'translateSelection', text: info.selectionText }, true);
        break;
      case MENU.pdfLink:
      case MENU.pdfPage: {
        const target = info.menuItemId === MENU.pdfLink ? info.linkUrl : tab?.url;
        const url = browser.runtime.getURL('/pdf.html') + (target ? `?file=${encodeURIComponent(target)}` : '');
        await browser.tabs.create({ url, index: tab ? tab.index + 1 : undefined });
        break;
      }
    }
  });

  // ---- Klavye kısayolları -----------------------------------------------------
  browser.commands.onCommand.addListener(async (command, tab) => {
    const active = tab?.id !== undefined ? tab : (await browser.tabs.query({ active: true, currentWindow: true }))[0];
    const tabId = active?.id;
    if (tabId === undefined) return;
    await ensureContentScript(tabId);
    switch (command) {
      case 'toggle-translate':
        await sendTabAction(tabId, { type: 'toggle' });
        break;
      case 'toggle-display-mode':
        await sendTabAction(tabId, { type: 'toggleMode' });
        break;
      case 'translate-selection':
        await sendTabAction(tabId, { type: 'translateSelection' }, true);
        break;
      case 'toggle-hover':
        await sendTabAction(tabId, { type: 'toggleHover' }, true);
        break;
    }
  });

  // Uzantı yüklenmeden önce açılmış sekmelere içerik betiğini sonradan ekler.
  async function ensureContentScript(tabId: number) {
    const alive = await sendTabAction(tabId, { type: 'getState' }, true);
    if (alive) return;
    try {
      await browser.scripting.executeScript({ target: { tabId, allFrames: true }, files: ['/content-scripts/content.js'] });
    } catch {
      /* korumalı sayfa */
    }
  }
});

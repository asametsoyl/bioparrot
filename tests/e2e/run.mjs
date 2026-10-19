// Uçtan uca test: derlenmiş eklentiyi yerel Chrome'a (geçici profil) yükler ve örnek makale sayfasını çevirir.
// Kullanım: npm run build && node tests/e2e/run.mjs
import fs from 'node:fs';
import http from 'node:http';
import os from 'node:os';
import path from 'node:path';
import puppeteer from 'puppeteer-core';

const ROOT = path.resolve(import.meta.dirname, '../..');
const EXT = path.join(ROOT, '.output/chrome-mv3');
const OUT = path.join(ROOT, 'test-results');
const CHROME = process.env.CHROME_PATH ?? '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
fs.mkdirSync(OUT, { recursive: true });

const server = http
  .createServer((req, res) => {
    const name = req.url === '/' ? 'article.html' : req.url.slice(1).split('?')[0];
    fs.readFile(path.join(ROOT, 'tests/e2e/fixtures', name), (err, data) => {
      if (err) return res.writeHead(404).end();
      res.writeHead(200, { 'Content-Type': name.endsWith('.pdf') ? 'application/pdf' : 'text/html; charset=utf-8' }).end(data);
    });
  })
  .listen(0);
const port = server.address().port;
const base = `http://localhost:${port}`;

const results = [];
const check = (name, ok, detail = '') => {
  results.push({ name, ok, detail });
  console.log(`${ok ? '✓' : '✗'} ${name}${detail ? ` — ${detail}` : ''}`);
};
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

const userDataDir = fs.mkdtempSync(path.join(os.tmpdir(), 'bioparrot-e2e-'));
const browser = await puppeteer.launch({
  executablePath: CHROME,
  headless: true,
  pipe: true,
  enableExtensions: [EXT],
  userDataDir,
  args: ['--no-first-run', '--no-default-browser-check', '--lang=en-US', '--window-size=1280,900'],
  defaultViewport: { width: 1280, height: 900 },
});

try {
  const swTarget = await browser.waitForTarget((t) => t.type() === 'service_worker' && t.url().startsWith('chrome-extension://'), { timeout: 15000 });
  const extId = new URL(swTarget.url()).host;
  const sw = await swTarget.worker();
  check('Eklenti yüklendi', !!extId, extId);

  // Ayar: hedef dil Türkçe, motor Google (ücretsiz)
  await sw.evaluate(async () => {
    const { settings } = await chrome.storage.local.get('settings');
    await chrome.storage.local.set({ settings: { ...(settings ?? {}), targetLanguage: 'tr', engineId: 'google' } });
  });

  // Kurulumda açılan karşılama sekmesini kapat; test sekmesi ön planda olmalı (arka plan sekmeleri çizim yapmaz).
  for (const p of await browser.pages()) if (p.url().includes('options.html')) await p.close();
  const page = await browser.newPage();
  page.on('console', (m) => m.type() === 'error' && console.log('  [page]', m.text()));
  await page.goto(base + '/', { waitUntil: 'load' });
  await page.bringToFront();
  await page.waitForSelector('[data-bp-ui="floating-ball"]', { timeout: 10000 });
  check('Yüzen düğme göründü', true);

  await page.click('pierce/.bp-ball');
  await page.waitForSelector('font.bp-wrap .bp-text', { timeout: 30000 });
  await sleep(5000);

  const info = await page.evaluate(() => {
    const textOf = (sel) => document.querySelector(sel)?.querySelector('font.bp-wrap .bp-text')?.textContent ?? null;
    const wraps = Array.from(document.querySelectorAll('font.bp-wrap .bp-text'));
    return {
      count: wraps.length,
      h1: textOf('h1'),
      abstract: textOf('#abstract p'),
      abstractLinks: document.querySelectorAll('#abstract p font.bp-wrap a[href="#fig1"]').length,
      abstractBold: document.querySelectorAll('#abstract p font.bp-wrap b').length,
      qpcr: textOf('#p-qpcr'),
      qpcrCode: document.querySelectorAll('#p-qpcr font.bp-wrap code').length,
      qpcrSup: document.querySelectorAll('#p-qpcr font.bp-wrap sup').length,
      math: document.querySelectorAll('#p-math font.bp-wrap math').length,
      pre: document.querySelectorAll('pre font.bp-wrap').length,
      refs: document.querySelectorAll('#references li font.bp-wrap').length,
      refsHeading: textOf('#references h2'),
      nav: document.querySelectorAll('nav font.bp-wrap').length,
      footer: document.querySelectorAll('footer font.bp-wrap').length,
      li: document.querySelectorAll('li font.bp-wrap').length,
      td: textOf('td:last-child'),
      card: textOf('.card'),
      mixed: Array.from(document.querySelectorAll('.card:nth-child(2) font.bp-wrap .bp-text')).map((e) => e.textContent),
      dyn: textOf('#p-dyn'),
      errors: document.querySelectorAll('font.bp-error').length,
      loading: document.querySelectorAll('font.bp-loading').length,
    };
  });
  console.log(JSON.stringify(info, null, 2));
  check('Paragraflar çevrildi', info.count >= 10, `${info.count} çeviri`);
  check('Başlık çevrildi', !!info.h1 && info.h1 !== 'CRISPR-Cas9 editing of TP53 reveals a role in organoid growth', info.h1 ?? '');
  check('Özet çevrildi', !!info.abstract, info.abstract ?? '');
  check('Bağlantı çeviride korundu', info.abstractLinks === 1);
  check('Kalın yazı çeviride korundu', info.abstractBold >= 1);
  check('Satır içi kod korundu', info.qpcrCode === 1);
  check('Atıf üst simgesi korundu', info.qpcrSup === 1);
  check('Matematik ifadesi korundu', info.math === 1);
  check('Kod bloğu çevrilmedi', info.pre === 0);
  check('Kaynakça çevrilmedi', info.refs === 0);
  check('Gezinme menüsü atlandı (ana içerik modu)', info.nav === 0);
  check('Liste öğeleri çevrildi', info.li >= 2);
  check('Tablo hücresi çevrildi', !!info.td);
  check('Karma içerik ayrı ayrı çevrildi', info.mixed.length >= 2, JSON.stringify(info.mixed));
  check('Sonradan eklenen paragraf çevrildi', !!info.dyn, info.dyn ?? '');
  check('Hata yok', info.errors === 0 && info.loading === 0, `hata ${info.errors}, bekleyen ${info.loading}`);
  await page.screenshot({ path: path.join(OUT, '1-bilingual.png'), fullPage: true });

  // Yalnız çeviri modu
  await sw.evaluate(async () => {
    const { settings } = await chrome.storage.local.get('settings');
    await chrome.storage.local.set({ settings: { ...settings, displayMode: 'translation' } });
  });
  await sleep(600);
  const onlyInfo = await page.evaluate(() => {
    const p = document.querySelector('#abstract p');
    const li = document.querySelector('li');
    return {
      visibleText: p.innerText.slice(0, 160),
      liMarker: getComputedStyle(li).display,
      pMargin: getComputedStyle(p).marginTop,
    };
  });
  check(
    'Yalnız çeviri modu uygulandı',
    !onlyInfo.visibleText.includes('Mitochondria are') && onlyInfo.visibleText.includes('Mitokondri') && onlyInfo.pMargin !== '0px',
    JSON.stringify(onlyInfo),
  );
  await page.screenshot({ path: path.join(OUT, '2-translation-only.png'), fullPage: true });
  await sw.evaluate(async () => {
    const { settings } = await chrome.storage.local.get('settings');
    await chrome.storage.local.set({ settings: { ...settings, displayMode: 'bilingual' } });
  });
  await sleep(400);

  // Orijinali göster
  await page.click('pierce/.bp-ball');
  await sleep(600);
  const restored = await page.evaluate(() => document.querySelectorAll('font.bp-wrap').length);
  check('Orijinale dönüldü', restored === 0, `${restored} sarmalayıcı kaldı`);

  // Fareyle çeviri (Shift)
  const box = await page.$eval('#p-qpcr', (el) => {
    const r = el.getBoundingClientRect();
    return { x: r.left + 40, y: r.top + 8 };
  });
  await page.mouse.move(box.x, box.y);
  await page.keyboard.down('Shift');
  await page.keyboard.up('Shift');
  await page.waitForSelector('#p-qpcr font.bp-wrap .bp-text', { timeout: 15000 }).catch(() => {});
  const hovered = await page.evaluate(() => ({
    target: !!document.querySelector('#p-qpcr font.bp-wrap .bp-text'),
    others: document.querySelectorAll('font.bp-wrap').length,
  }));
  check('Fareyle tek paragraf çevrildi', hovered.target && hovered.others === 1, JSON.stringify(hovered));

  // Seçim çevirisi (varsayılan: küçük nokta → paragrafın altına)
  const selectText = async (selector, from, to) => {
    await page.evaluate(
      (sel, a, b) => {
        const el = document.querySelector(sel);
        const walker = document.createTreeWalker(el, NodeFilter.SHOW_TEXT);
        const first = walker.nextNode();
        const range = document.createRange();
        range.setStart(first, a);
        if (b === -1) range.setEndAfter(el.lastChild);
        else range.setEnd(first, b);
        getSelection().removeAllRanges();
        getSelection().addRange(range);
        const r = range.getBoundingClientRect();
        el.dispatchEvent(new MouseEvent('mouseup', { bubbles: true, button: 0, clientX: r.right, clientY: r.bottom }));
      },
      selector,
      from,
      to,
    );
    await sleep(300);
  };
  await page.evaluate(() => document.querySelectorAll('font.bp-wrap').forEach((w) => w.remove()));
  await selectText('#abstract p', 0, 'Mitochondria are the'.length);
  const dot = await page.$('pierce/.bp-sel-icon');
  check('Seçimde küçük nokta göründü', !!dot);
  if (dot) {
    await dot.click();
    await page.waitForSelector('#abstract p .bp-sel-wrap .bp-text', { timeout: 15000 }).catch(() => {});
    const inline = await page.evaluate(() => document.querySelector('#abstract p .bp-sel-wrap .bp-text')?.textContent ?? null);
    check('Seçimin çevirisi paragrafın altına eklendi', !!inline, inline ?? '');
    await page.screenshot({ path: path.join(OUT, '3-selection-inline.png') });
  }
  // Tüm paragraf seçilince biçimlendirme korunarak çevrilir
  await selectText('#p-qpcr', 0, -1);
  const dot2 = await page.$('pierce/.bp-sel-icon');
  if (dot2) await dot2.click();
  await page.waitForSelector('#p-qpcr font.bp-wrap .bp-text', { timeout: 15000 }).catch(() => {});
  const full = await page.evaluate(() => ({ text: document.querySelector('#p-qpcr font.bp-wrap .bp-text')?.textContent, code: document.querySelectorAll('#p-qpcr font.bp-wrap code').length }));
  check('Paragrafın tamamı biçimiyle çevrildi', !!full.text && full.code === 1, JSON.stringify(full));

  // Kart modu
  await sw.evaluate(async () => {
    const { settings } = await chrome.storage.local.get('settings');
    await chrome.storage.local.set({ settings: { ...settings, selection: { ...settings.selection, display: 'card' } } });
  });
  await sleep(400);
  await selectText('#p-dyn', 0, 'This paragraph was added'.length);
  const dot3 = await page.$('pierce/.bp-sel-icon');
  if (dot3) {
    await dot3.click();
    await page
      .waitForFunction(
        () => {
          const r = document.querySelector('[data-bp-ui="selection"]')?.shadowRoot?.querySelector('.bp-result');
          return r && r.textContent && !r.querySelector('.bp-spin');
        },
        { timeout: 15000 },
      )
      .catch(() => {});
    const selText = await page.evaluate(() => document.querySelector('[data-bp-ui="selection"]')?.shadowRoot?.querySelector('.bp-result')?.textContent ?? null);
    check('Kart modunda çeviri', !!selText, selText ?? '');
    await page.screenshot({ path: path.join(OUT, '3-selection-card.png') });
  } else check('Kart modunda çeviri', false, 'nokta yok');

  // Giriş kutusu: üç kez boşluk
  await page.keyboard.press('Escape');
  await page.click('#ta');
  await page.keyboard.type('Merhaba dünya, bu bir denemedir.');
  await page.keyboard.press('Space');
  await page.keyboard.press('Space');
  await page.keyboard.press('Space');
  await page.waitForFunction(() => !/Merhaba/.test(document.getElementById('ta').value), { timeout: 15000 }).catch(() => {});
  const ta = await page.$eval('#ta', (el) => el.value);
  check('Giriş kutusu çevirisi (TR → EN)', /hello|world|test/i.test(ta), ta);

  // Otomatik çeviri seçeneği: açıkken sayfa kendiliğinden çevrilir
  await sw.evaluate(async () => {
    const { settings } = await chrome.storage.local.get('settings');
    await chrome.storage.local.set({ settings: { ...settings, page: { ...settings.page, autoTranslateAll: true } } });
  });
  await page.reload({ waitUntil: 'load' });
  await page.bringToFront();
  await page.waitForSelector('font.bp-wrap .bp-text', { timeout: 20000 }).catch(() => {});
  await sleep(4000);
  const autoCount = await page.evaluate(() => document.querySelectorAll('font.bp-wrap .bp-text').length);
  check('“Tüm sayfaları otomatik çevir” açıkken sayfa kendiliğinden çevrildi', autoCount > 3, `${autoCount} çeviri`);
  await sw.evaluate(async () => {
    const { settings } = await chrome.storage.local.get('settings');
    await chrome.storage.local.set({ settings: { ...settings, page: { ...settings.page, autoTranslateAll: false } } });
  });

  // Popup ve ayarlar sayfası ekran görüntüleri
  await page.bringToFront();
  const popup = await browser.newPage();
  await popup.setViewport({ width: 360, height: 640 });
  await popup.goto(`chrome-extension://${extId}/popup.html`);
  await sleep(800);
  await popup.screenshot({ path: path.join(OUT, '4-popup.png') });
  const options = await browser.newPage();
  await options.goto(`chrome-extension://${extId}/options.html#engines`);
  await sleep(800);
  await options.screenshot({ path: path.join(OUT, '5-options-engines.png'), fullPage: false });
  await options.goto(`chrome-extension://${extId}/options.html#style`);
  await sleep(800);
  await options.screenshot({ path: path.join(OUT, '6-options-style.png'), fullPage: false });
} catch (err) {
  check('Beklenmeyen hata', false, String(err?.stack ?? err));
} finally {
  await browser.close();
  server.close();
  fs.rmSync(userDataDir, { recursive: true, force: true });
}

const failed = results.filter((r) => !r.ok);
console.log(`\n${results.length - failed.length}/${results.length} kontrol geçti.`);
process.exit(failed.length ? 1 : 0);

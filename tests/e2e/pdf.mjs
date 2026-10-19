// PDF çevirmeni uçtan uca testi: node tests/e2e/pdf.mjs <pdf-dosyası>
import fs from 'node:fs';
import http from 'node:http';
import os from 'node:os';
import path from 'node:path';
import puppeteer from 'puppeteer-core';

const ROOT = path.resolve(import.meta.dirname, '../..');
const EXT = path.join(ROOT, '.output/chrome-mv3');
const OUT = path.join(ROOT, 'test-results');
const PDF = path.resolve(process.argv[2] ?? path.join(ROOT, 'tests/e2e/fixtures/sample-paper.pdf'));
const CHROME = process.env.CHROME_PATH ?? '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
if (!fs.existsSync(PDF)) {
  console.error('PDF dosyası bulunamadı:', PDF);
  process.exit(2);
}
fs.mkdirSync(OUT, { recursive: true });
const downloads = fs.mkdtempSync(path.join(os.tmpdir(), 'bioparrot-dl-'));

const server = http
  .createServer((req, res) => {
    const file = req.url.startsWith('/out/') ? path.join(downloads, decodeURIComponent(req.url.slice(5))) : PDF;
    fs.readFile(file, (err, data) => (err ? res.writeHead(404).end() : res.writeHead(200, { 'Content-Type': 'application/pdf' }).end(data)));
  })
  .listen(0);
const base = `http://localhost:${server.address().port}`;
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
const userDataDir = fs.mkdtempSync(path.join(os.tmpdir(), 'bioparrot-pdf-'));
const browser = await puppeteer.launch({
  executablePath: CHROME,
  headless: true,
  pipe: true,
  enableExtensions: [EXT],
  userDataDir,
  args: ['--no-first-run', '--lang=en-US'],
  defaultViewport: { width: 1600, height: 1000 },
  protocolTimeout: 120000,
});
let ok = true;
try {
  const swTarget = await browser.waitForTarget((t) => t.type() === 'service_worker' && t.url().startsWith('chrome-extension://'));
  const extId = new URL(swTarget.url()).host;
  const sw = await swTarget.worker();
  await sw.evaluate(async () => {
    const { settings } = await chrome.storage.local.get('settings');
    await chrome.storage.local.set({ settings: { ...(settings ?? {}), targetLanguage: 'tr', engineId: 'google' } });
  });
  for (const p of await browser.pages()) if (p.url().includes('options.html')) await p.close();
  const page = await browser.newPage();
  page.on('console', (m) => ['error', 'warn'].includes(m.type()) && console.log('  [page]', m.text().slice(0, 200)));
  page.on('pageerror', (e) => console.log('  [pageerror]', e.message));
  const cdp = await page.target().createCDPSession();
  await cdp.send('Browser.setDownloadBehavior', { behavior: 'allow', downloadPath: downloads });
  const t0 = Date.now();
  await page.goto(`chrome-extension://${extId}/pdf.html?file=${encodeURIComponent(base + '/paper.pdf')}`);
  await page.bringToFront();
  await page.waitForFunction(() => /tamamlandı|complete/i.test(document.querySelector('.pdf-status')?.textContent ?? ''), { timeout: 180000, polling: 1000 });
  console.log(`✓ Çeviri tamamlandı (${((Date.now() - t0) / 1000).toFixed(1)} sn)`);
  const stats = await page.evaluate(() => ({
    status: document.querySelector('.pdf-status')?.textContent,
    blocks: document.querySelectorAll('.pv-block').length,
    sample: Array.from(document.querySelectorAll('.pv-block')).slice(0, 6).map((b) => b.textContent.slice(0, 120)),
  }));
  console.log(JSON.stringify(stats, null, 2));
  if (stats.blocks < 20) ok = false;
  await page.evaluate(() => document.querySelector('.pdf-main').scrollTo(0, 0));
  await sleep(1500);
  await page.screenshot({ path: path.join(OUT, 'pdf-1-page1.png') });
  await page.evaluate(() => {
    const row = document.querySelector('.pv-row[data-page="3"]');
    row?.scrollIntoView();
  });
  await sleep(2000);
  await page.screenshot({ path: path.join(OUT, 'pdf-2-page3.png') });

  // Dışa aktar: iki dilli PDF
  await page.click('.export .btn-primary');
  await page.click('.export .menu button');
  let exported = null;
  for (let i = 0; i < 120 && !exported; i++) {
    await sleep(1000);
    exported = fs.readdirSync(downloads).find((f) => f.endsWith('.pdf'));
  }
  if (exported) {
    const size = fs.statSync(path.join(downloads, exported)).size;
    console.log(`✓ Dışa aktarıldı: ${exported} (${(size / 1024).toFixed(0)} KB)`);
    fs.copyFileSync(path.join(downloads, exported), path.join(OUT, exported));
    // Dışa aktarılan PDF'i görüntüleyicide aç ve ilk sayfanın görüntüsünü al
    const viewer = await browser.newPage();
    await viewer.goto(`chrome-extension://${extId}/pdf.html?file=${encodeURIComponent(`${base}/out/${exported}`)}`);
    await viewer.bringToFront();
    await viewer.waitForSelector('.pv-page canvas', { timeout: 60000 });
    await viewer.evaluate(() => document.querySelector('.seg button:nth-child(3)')?.click());
    await sleep(4000);
    await viewer.screenshot({ path: path.join(OUT, 'pdf-3-exported.png') });
  } else {
    console.log('✗ Dışa aktarılan dosya bulunamadı');
    ok = false;
  }
} catch (err) {
  console.log('✗ Hata:', err?.stack ?? err);
  ok = false;
} finally {
  await browser.close();
  server.close();
}
process.exit(ok ? 0 : 1);

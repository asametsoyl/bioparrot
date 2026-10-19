// Test ve README için örnek iki sütunlu makale PDF'i üretir: node scripts/make-sample-pdf.mjs
import path from 'node:path';
import puppeteer from 'puppeteer-core';

const ROOT = path.resolve(import.meta.dirname, '..');
const CHROME = process.env.CHROME_PATH ?? '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome';
const browser = await puppeteer.launch({ executablePath: CHROME, headless: true });
const page = await browser.newPage();
await page.goto('file://' + path.join(ROOT, 'tests/e2e/fixtures/sample-paper.html'));
await page.pdf({ path: path.join(ROOT, 'tests/e2e/fixtures/sample-paper.pdf'), format: 'A4', printBackground: true, preferCSSPageSize: true });
await browser.close();
console.log('✔ tests/e2e/fixtures/sample-paper.pdf');

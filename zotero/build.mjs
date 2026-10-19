// Zotero eklentisini derler ve .xpi paketi oluşturur: npm run build:zotero
import { execFileSync } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import react from '@vitejs/plugin-react';
import { build } from 'vite';

const ROOT = path.resolve(import.meta.dirname, '..');
const ADDON = path.join(ROOT, 'zotero/addon');
const CONTENT = path.join(ADDON, 'content');
const pkg = JSON.parse(fs.readFileSync(path.join(ROOT, 'package.json'), 'utf8'));

const shared = {
  configFile: false,
  root: ROOT,
  publicDir: false,
  logLevel: 'error',
  resolve: {
    alias: {
      '@': path.join(ROOT, 'src'),
      // Zotero 7 (Firefox 115) ile uyum için pdf.js'in eski tarayıcı sürümü
      'pdfjs-dist': path.join(ROOT, 'node_modules/pdfjs-dist/legacy/build/pdf.mjs'),
    },
  },
  define: { 'process.env.NODE_ENV': '"production"' },
};

// 1) Ana eklenti betiği (bootstrap.js tarafından yüklenir)
await build({
  ...shared,
  build: {
    outDir: path.join(CONTENT, 'scripts'),
    emptyOutDir: true,
    target: 'firefox115',
    minify: false,
    lib: { entry: path.join(ROOT, 'zotero/src/index.ts'), formats: ['iife'], name: '__bioparrot', fileName: () => 'bioparrot.js' },
  },
});

// 2) İki dilli PDF penceresi (React)
await build({
  ...shared,
  plugins: [react()],
  build: {
    outDir: path.join(CONTENT, 'pdf'),
    emptyOutDir: true,
    target: 'firefox115',
    cssCodeSplit: false,
    lib: {
      entry: path.join(ROOT, 'zotero/src/pdf-dialog.tsx'),
      formats: ['iife'],
      name: '__bioparrotPdf',
      fileName: () => 'pdf-dialog.js',
      cssFileName: 'pdf-dialog',
    },
  },
});

// 3) pdf.js işçisi ve yazı tipleri
const pdfDir = path.join(CONTENT, 'pdf');
fs.mkdirSync(path.join(pdfDir, 'fonts'), { recursive: true });
fs.copyFileSync(path.join(ROOT, 'node_modules/pdfjs-dist/legacy/build/pdf.worker.min.mjs'), path.join(pdfDir, 'pdf.worker.min.mjs'));
const fonts = {
  'NotoSans-Regular.ttf': '@expo-google-fonts/noto-sans/400Regular/NotoSans_400Regular.ttf',
  'NotoSans-Bold.ttf': '@expo-google-fonts/noto-sans/700Bold/NotoSans_700Bold.ttf',
  'NotoSansMath-Regular.ttf': '@expo-google-fonts/noto-sans-math/400Regular/NotoSansMath_400Regular.ttf',
};
for (const [name, src] of Object.entries(fonts)) fs.copyFileSync(path.join(ROOT, 'node_modules', src), path.join(pdfDir, 'fonts', name));

// 4) Sürüm numarasını eşitle ve paketle
const manifestPath = path.join(ADDON, 'manifest.json');
const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
manifest.version = pkg.version;
fs.writeFileSync(manifestPath, JSON.stringify(manifest, null, 2) + '\n');

// exFAT disklerde oluşan ._* dosyalarını temizle
const clean = (dir) => {
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, entry.name);
    if (entry.name.startsWith('._') || entry.name === '.DS_Store') fs.rmSync(p, { force: true });
    else if (entry.isDirectory()) clean(p);
  }
};
clean(ADDON);

const outDir = path.join(ROOT, '.output');
fs.mkdirSync(outDir, { recursive: true });
const xpi = path.join(outDir, `bioparrot-zotero-${pkg.version}.xpi`);
fs.rmSync(xpi, { force: true });
execFileSync('zip', ['-r', '-X', '-q', xpi, '.', '-x', '._*', '*/._*', '.DS_Store'], { cwd: ADDON });
console.log(`✔ Zotero eklentisi: ${path.relative(ROOT, xpi)} (${(fs.statSync(xpi).size / 1024 / 1024).toFixed(1)} MB)`);

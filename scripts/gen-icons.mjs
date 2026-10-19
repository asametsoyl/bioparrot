// Logo SVG'sinden eklenti simgelerini üretir: node scripts/gen-icons.mjs
import { mkdir, writeFile } from 'node:fs/promises';
import sharp from 'sharp';

const brand = await import('../src/core/brand.ts').catch(() => null);
if (!brand) throw new Error('src/core/brand.ts yüklenemedi (Node 22.6+ gerekir)');

await mkdir('public/icon', { recursive: true });
await writeFile('public/logo.svg', brand.LOGO_SVG);
for (const size of [16, 32, 48, 96, 128]) {
  await sharp(Buffer.from(brand.LOGO_SVG), { density: 512 }).resize(size, size).png().toFile(`public/icon/${size}.png`);
}
await sharp(Buffer.from(brand.LOGO_SVG), { density: 512 }).resize(512, 512).png().toFile('docs/logo-512.png').catch(async () => {
  await mkdir('docs', { recursive: true });
  await sharp(Buffer.from(brand.LOGO_SVG), { density: 512 }).resize(512, 512).png().toFile('docs/logo-512.png');
});
console.log('Simgeler üretildi.');

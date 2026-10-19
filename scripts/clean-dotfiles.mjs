// macOS'un exFAT/FAT disklerde oluşturduğu "._*" ve .DS_Store dosyalarını siler (diğer disklerde etkisizdir).
import fs from 'node:fs';
import path from 'node:path';

function clean(dir) {
  if (!fs.existsSync(dir)) return 0;
  let removed = 0;
  for (const entry of fs.readdirSync(dir, { withFileTypes: true })) {
    const p = path.join(dir, entry.name);
    if (entry.name.startsWith('._') || entry.name === '.DS_Store') {
      fs.rmSync(p, { force: true, recursive: true });
      removed++;
    } else if (entry.isDirectory() && entry.name !== 'node_modules' && entry.name !== '.git') removed += clean(p);
  }
  return removed;
}

const dirs = process.argv.slice(2);
let total = 0;
for (const d of dirs.length ? dirs : ['.output']) total += clean(path.resolve(d));
if (total) console.log(`clean-dotfiles: ${total} dosya silindi`);

# Katkıda bulunma / Contributing

Teşekkürler! / Thank you!

1. Depoyu çatallayın ve bir dal açın / Fork the repo and create a branch.
2. `npm install`, ardından değişikliklerinizi yapın / then make your changes.
3. Göndermeden önce / Before sending a pull request:
   ```bash
   npm run compile && npm test && npm run build
   ```
   Arayüzü veya sayfa çevirisini etkileyen değişikliklerde `npm run test:e2e` ve PDF ile ilgili değişikliklerde `npm run test:pdf` çalıştırın (yerel Chrome gerekir).
4. Arayüz metinleri `src/core/i18n/tr.ts` ve `en.ts` dosyalarındadır; yeni metinleri iki dile de ekleyin.
   UI strings live in `src/core/i18n/tr.ts` and `en.ts`; add new strings to both.
5. Yeni bir çeviri motoru eklemek için `src/core/engines/` altına bir modül yazıp `index.ts` ve `meta.ts` dosyalarına kaydedin.
   To add a translation service, write a module in `src/core/engines/` and register it in `index.ts` and `meta.ts`.

PDF düzen analizini ayıklamak için / To debug PDF layout analysis:
```bash
BP_PDF=/yol/makale.pdf BP_PAGES=3 BP_SHOW=1,2 node node_modules/vitest/vitest.mjs run --reporter=verbose
```

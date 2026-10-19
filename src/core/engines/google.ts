import { httpJson } from '../http';
import { fromEngineCode, toGoogleCode } from '../languages';
import { markupToPlain, plainToMarkup } from '../markup';
import { TranslateError } from '../types';
import type { Engine } from './types';

/**
 * Google Translate (anahtarsız web uç noktası). HTML biçimini destekler ve
 * çoklu `q` parametresiyle toplu çeviri yapar.
 * Yanıt: kaynak "auto" ise [[metin, dil], …], aksi halde [metin, …].
 */
export const googleEngine: Engine = {
  async translate(input) {
    const sl = toGoogleCode(input.from);
    const tl = toGoogleCode(input.to);
    const body = new URLSearchParams();
    for (const text of input.texts) body.append('q', input.rich ? text : plainToMarkup(text));

    const url = `https://translate.googleapis.com/translate_a/t?client=gtx&sl=${encodeURIComponent(sl)}&tl=${encodeURIComponent(tl)}&format=html`;
    const data = await httpJson<unknown>(url, {
      method: 'POST',
      headers: { 'Content-Type': 'application/x-www-form-urlencoded;charset=UTF-8' },
      body: body.toString(),
      signal: input.signal,
      timeoutMs: 30_000,
    });

    const items = Array.isArray(data) ? data : [data];
    if (items.length !== input.texts.length) {
      throw new TranslateError(`Google beklenmeyen sayıda sonuç döndürdü (${items.length}/${input.texts.length}).`, 'parse');
    }
    const translations: string[] = [];
    const detected: (string | undefined)[] = [];
    for (const item of items) {
      let text: unknown = item;
      let lang: string | undefined;
      if (Array.isArray(item)) {
        text = item[0];
        lang = typeof item[1] === 'string' ? item[1] : undefined;
      }
      if (typeof text !== 'string') throw new TranslateError('Google yanıtı çözümlenemedi.', 'parse');
      translations.push(input.rich ? text : markupToPlain(text));
      detected.push(lang ? fromEngineCode(lang, 'google') : undefined);
    }
    return { translations, detected };
  },
};

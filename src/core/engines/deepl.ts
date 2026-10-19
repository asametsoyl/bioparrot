import { httpJson } from '../http';
import { fromEngineCode, toDeeplSource, toDeeplTarget } from '../languages';
import { markupToPlain, plainToMarkup } from '../markup';
import { TranslateError } from '../types';
import type { Engine } from './types';

interface DeeplResponse {
  translations?: { detected_source_language?: string; text: string }[];
}

/** DeepL API (Free anahtarları ":fx" ile biter). */
export const deeplEngine: Engine = {
  async translate(input, config) {
    const key = config.apiKey?.trim();
    if (!key) throw new TranslateError('DeepL için API anahtarı girin.', 'config');
    const host = key.endsWith(':fx') ? 'https://api-free.deepl.com' : 'https://api.deepl.com';

    const body: Record<string, unknown> = {
      text: input.texts.map((t) => (input.rich ? t : plainToMarkup(t))),
      target_lang: toDeeplTarget(input.to),
      tag_handling: 'html',
      preserve_formatting: true,
    };
    const source = toDeeplSource(input.from);
    if (source) body.source_lang = source;
    const context = [input.context?.title, input.context?.summary].filter(Boolean).join('\n').slice(0, 2000);
    if (context) body.context = context;

    const data = await httpJson<DeeplResponse>(`${host}/v2/translate`, {
      method: 'POST',
      headers: { Authorization: `DeepL-Auth-Key ${key}` },
      json: body,
      signal: input.signal,
      timeoutMs: 45_000,
    });
    const list = data.translations;
    if (!Array.isArray(list)) throw new TranslateError('DeepL yanıtı çözümlenemedi.', 'parse');
    return {
      translations: list.map((t) => (input.rich ? t.text : markupToPlain(t.text))),
      detected: list.map((t) => fromEngineCode(t.detected_source_language, 'deepl') || undefined),
    };
  },
};

/** DeepLX: DeepL'i kendi sunucunuz üzerinden kullanan açık kaynak proxy. Düz metin, istek başına tek metin. */
export const deeplxEngine: Engine = {
  async translate(input, config) {
    const url = config.baseUrl?.trim();
    if (!url) throw new TranslateError('DeepLX sunucu adresini girin (örn. http://localhost:1188/translate).', 'config');
    const headers: Record<string, string> = {};
    if (config.apiKey?.trim()) headers.Authorization = `Bearer ${config.apiKey.trim()}`;
    const translations: (string | undefined)[] = [];
    for (const text of input.texts) {
      const data = await httpJson<{ code?: number; data?: string; message?: string }>(url, {
        method: 'POST',
        headers,
        json: {
          text: input.rich ? markupToPlain(text) : text,
          source_lang: toDeeplSource(input.from) ?? 'auto',
          target_lang: toDeeplTarget(input.to),
        },
        signal: input.signal,
        timeoutMs: 30_000,
      });
      if (data.code !== undefined && data.code !== 200) {
        throw new TranslateError(`DeepLX hatası: ${data.message ?? data.code}`, 'server');
      }
      translations.push(typeof data.data === 'string' ? (input.rich ? plainToMarkup(data.data) : data.data) : undefined);
    }
    return { translations };
  },
};

import { httpJson } from '../http';
import { fromEngineCode, toMicrosoftCode } from '../languages';
import { markupToPlain, plainToMarkup } from '../markup';
import { TranslateError } from '../types';
import type { Engine } from './types';

interface MsItem {
  detectedLanguage?: { language: string; score: number };
  translations?: { text: string; to: string }[];
}

/** Microsoft Translator (Azure AI Translator, resmi API anahtarıyla). */
export const microsoftEngine: Engine = {
  async translate(input, config) {
    const key = config.apiKey?.trim();
    if (!key) throw new TranslateError('Microsoft Translator için API anahtarı girin.', 'config');
    const base = (config.baseUrl?.trim() || 'https://api.cognitive.microsofttranslator.com').replace(/\/+$/, '');
    const params = new URLSearchParams({ 'api-version': '3.0', to: toMicrosoftCode(input.to), textType: 'html' });
    if (input.from !== 'auto') params.set('from', toMicrosoftCode(input.from));

    const headers: Record<string, string> = { 'Ocp-Apim-Subscription-Key': key };
    if (config.region?.trim()) headers['Ocp-Apim-Subscription-Region'] = config.region.trim();

    const data = await httpJson<MsItem[]>(`${base}/translate?${params}`, {
      method: 'POST',
      headers,
      json: input.texts.map((t) => ({ Text: input.rich ? t : plainToMarkup(t) })),
      signal: input.signal,
      timeoutMs: 30_000,
    });
    if (!Array.isArray(data)) throw new TranslateError('Microsoft yanıtı çözümlenemedi.', 'parse');
    return {
      translations: data.map((item) => {
        const text = item.translations?.[0]?.text;
        return text === undefined ? undefined : input.rich ? text : markupToPlain(text);
      }),
      detected: data.map((item) => fromEngineCode(item.detectedLanguage?.language, 'microsoft') || undefined),
    };
  },
};

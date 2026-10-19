import { httpJson, parseJsonObject } from '../http';
import { TranslateError, type EngineConfig } from '../types';
import { deepMerge, llmTranslate, resolveTemperature } from './llm';
import type { Engine, EngineUsage } from './types';

const SAFETY_OFF = ['HARM_CATEGORY_HARASSMENT', 'HARM_CATEGORY_HATE_SPEECH', 'HARM_CATEGORY_SEXUALLY_EXPLICIT', 'HARM_CATEGORY_DANGEROUS_CONTENT'].map(
  (category) => ({ category, threshold: 'BLOCK_NONE' }),
);

async function generate(config: EngineConfig, system: string, user: string, temperature: number | undefined, signal?: AbortSignal) {
  const key = config.apiKey?.trim();
  if (!key) throw new TranslateError('Gemini için API anahtarı girin.', 'config');
  const base = (config.baseUrl?.trim() || 'https://generativelanguage.googleapis.com/v1beta').replace(/\/+$/, '');
  const model = (config.model?.trim() || 'gemini-flash-latest').replace(/^models\//, '');

  const generationConfig: Record<string, any> = {};
  if (temperature !== undefined) generationConfig.temperature = temperature;
  // Gemini 2.5 Flash'ta düşünmeyi kapatmak çeviriyi belirgin şekilde hızlandırır.
  if (/2\.5-flash/.test(model)) generationConfig.thinkingConfig = { thinkingBudget: 0 };

  let body: Record<string, any> = {
    systemInstruction: { parts: [{ text: system }] },
    contents: [{ role: 'user', parts: [{ text: user }] }],
    generationConfig,
    // Akademik/tıbbi metinlerin gereksiz yere engellenmemesi için.
    safetySettings: SAFETY_OFF,
  };
  body = deepMerge(body, parseJsonObject(config.extraBody, 'Ek parametreler'));

  const data = await httpJson<any>(`${base}/models/${model}:generateContent`, {
    method: 'POST',
    headers: { 'x-goog-api-key': key },
    json: body,
    signal,
    timeoutMs: 180_000,
  });
  if (data?.promptFeedback?.blockReason) {
    throw new TranslateError(`Gemini isteği engelledi (${data.promptFeedback.blockReason}).`, 'refused');
  }
  const candidate = data?.candidates?.[0];
  const parts: any[] = candidate?.content?.parts ?? [];
  const text = parts
    .filter((p) => typeof p?.text === 'string' && !p.thought)
    .map((p) => p.text)
    .join('');
  if (!text) {
    if (candidate?.finishReason && candidate.finishReason !== 'STOP') {
      throw new TranslateError(`Gemini yanıt üretmedi (${candidate.finishReason}).`, candidate.finishReason === 'SAFETY' ? 'refused' : 'parse');
    }
    throw new TranslateError('Gemini boş yanıt döndürdü.', 'parse');
  }
  const usage: EngineUsage = {
    input: data?.usageMetadata?.promptTokenCount,
    output: data?.usageMetadata?.candidatesTokenCount,
  };
  return { text, usage };
}

/** Google Gemini API. */
export const geminiEngine: Engine = {
  translate(input, config) {
    return llmTranslate(input, config, (prompt, cfg, signal) =>
      generate(cfg, prompt.system, prompt.user, resolveTemperature(cfg, input.expert.temperature), signal),
    );
  },
  complete(system, user, config, signal) {
    return generate(config, system, user, resolveTemperature(config, 0.2), signal);
  },
};

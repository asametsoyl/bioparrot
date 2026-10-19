import { httpJson, parseJsonObject } from '../http';
import { TranslateError, type EngineConfig } from '../types';
import { getPreset } from './meta';
import { deepMerge, llmTranslate, resolveTemperature } from './llm';
import type { Engine, EngineUsage } from './types';

function endpoint(config: EngineConfig): string {
  const base = (config.baseUrl ?? '').trim().replace(/\/+$/, '');
  if (!base) throw new TranslateError('API adresi (Base URL) girin.', 'config');
  return /\/chat\/completions/.test(base) ? base : `${base}/chat/completions`;
}

/** Sıcaklık parametresini desteklemeyen akıl yürütme modelleri. */
function isReasoningModel(model: string): boolean {
  const name = model.replace(/^.*\//, '').toLowerCase();
  return /^(o\d|gpt-5)/.test(name) || name.includes('reasoner');
}

async function chat(config: EngineConfig, system: string, user: string, temperature: number | undefined, signal?: AbortSignal) {
  const model = config.model?.trim();
  if (!model) throw new TranslateError('Model adını girin.', 'config');
  const headers: Record<string, string> = {};
  const key = config.apiKey?.trim();
  if (key) {
    if (/openai\.azure\.com|cognitiveservices\.azure\.com/.test(config.baseUrl ?? '')) headers['api-key'] = key;
    else headers.Authorization = `Bearer ${key}`;
  }
  Object.assign(headers, getPreset(config.preset)?.headers, parseJsonObject(config.customHeaders, 'Ek başlıklar'));

  let body: Record<string, any> = {
    model,
    messages: [
      { role: 'system', content: system },
      { role: 'user', content: user },
    ],
    stream: false,
  };
  if (temperature !== undefined && !isReasoningModel(model)) body.temperature = temperature;
  body = deepMerge(body, parseJsonObject(config.extraBody, 'Ek parametreler'));

  const data = await httpJson<any>(endpoint(config), { method: 'POST', headers, json: body, signal, timeoutMs: 180_000 });
  const choice = data?.choices?.[0];
  let content: unknown = choice?.message?.content;
  if (Array.isArray(content)) content = content.map((p: any) => (typeof p === 'string' ? p : p?.text ?? '')).join('');
  if (typeof content !== 'string' || (!content && choice?.finish_reason === 'content_filter')) {
    if (choice?.finish_reason === 'content_filter') throw new TranslateError('Model içerik filtresi nedeniyle yanıt vermedi.', 'refused');
    throw new TranslateError(`Model yanıtı çözümlenemedi: ${JSON.stringify(data).slice(0, 200)}`, 'parse');
  }
  const usage: EngineUsage = { input: data?.usage?.prompt_tokens, output: data?.usage?.completion_tokens };
  return { text: content, usage };
}

/** OpenAI ve OpenAI uyumlu tüm servisler (DeepSeek, OpenRouter, Ollama, Groq, Azure OpenAI…). */
export const openaiEngine: Engine = {
  translate(input, config) {
    return llmTranslate(input, config, (prompt, cfg, signal) =>
      chat(cfg, prompt.system, prompt.user, resolveTemperature(cfg, input.expert.temperature), signal),
    );
  },
  complete(system, user, config, signal) {
    return chat(config, system, user, resolveTemperature(config, 0.2), signal);
  },
};

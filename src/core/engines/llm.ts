import { buildPrompt, type BuiltPrompt } from '../prompt';
import { parseSegments } from '../segments';
import type { EngineConfig } from '../types';
import type { EngineInput, EngineOutput, EngineUsage } from './types';

export type CompleteFn = (prompt: BuiltPrompt, config: EngineConfig, signal?: AbortSignal) => Promise<{ text: string; usage?: EngineUsage }>;

/** Tüm yapay zeka motorları için ortak çeviri akışı: istem oluştur → çağır → parçaları ayrıştır. */
export async function llmTranslate(input: EngineInput, config: EngineConfig, complete: CompleteFn): Promise<EngineOutput> {
  const prompt = buildPrompt({
    texts: input.texts,
    from: input.from,
    to: input.to,
    rich: input.rich,
    expert: input.expert,
    context: input.context,
    glossary: input.glossary,
    academic: input.academic,
  });
  const { text, usage } = await complete(prompt, config, input.signal);
  return { translations: parseSegments(text, input.texts.length), usage };
}

export function resolveTemperature(config: EngineConfig, fallback?: number): number | undefined {
  const t = config.temperature;
  if (typeof t === 'number' && Number.isFinite(t)) return t;
  return fallback;
}

export function deepMerge<T extends Record<string, any>>(base: T, extra: Record<string, any>): T {
  const out: Record<string, any> = { ...base };
  for (const [k, v] of Object.entries(extra)) {
    if (v && typeof v === 'object' && !Array.isArray(v) && out[k] && typeof out[k] === 'object' && !Array.isArray(out[k])) {
      out[k] = deepMerge(out[k], v);
    } else {
      out[k] = v;
    }
  }
  return out as T;
}

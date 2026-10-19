import type { EngineType } from '../types';
import { anthropicEngine } from './anthropic';
import { deeplEngine, deeplxEngine } from './deepl';
import { geminiEngine } from './gemini';
import { googleEngine } from './google';
import { microsoftEngine } from './microsoft';
import { openaiEngine } from './openai';
import type { Engine } from './types';

const ENGINES: Record<EngineType, Engine> = {
  google: googleEngine,
  microsoft: microsoftEngine,
  deepl: deeplEngine,
  deeplx: deeplxEngine,
  openai: openaiEngine,
  anthropic: anthropicEngine,
  gemini: geminiEngine,
};

export function getEngine(type: EngineType): Engine {
  const engine = ENGINES[type];
  if (!engine) throw new Error(`Bilinmeyen çeviri motoru: ${type}`);
  return engine;
}

export * from './meta';
export type * from './types';

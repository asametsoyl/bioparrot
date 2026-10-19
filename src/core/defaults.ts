import { DEFAULT_CLAUDE_MODEL } from './engines/anthropic';
import { getPreset } from './engines/meta';
import { normalizeLang, getLanguage } from './languages';
import type { EngineConfig, Settings } from './types';

export const SETTINGS_VERSION = 1;

const MT_LIMITS = { maxConcurrent: 4, maxRps: 6, batchSize: 40, batchChars: 4000 };
const LLM_LIMITS = { maxConcurrent: 4, maxRps: 4, batchSize: 12, batchChars: 2400 };

function openaiPreset(id: string, presetId: string, name?: string, extra: Partial<EngineConfig> = {}): EngineConfig {
  const preset = getPreset(presetId)!;
  return {
    id,
    type: 'openai',
    name: name ?? preset.label,
    preset: presetId,
    baseUrl: preset.baseUrl,
    model: preset.model,
    apiKey: '',
    builtin: true,
    ...LLM_LIMITS,
    ...extra,
  };
}

export function defaultEngines(): Record<string, EngineConfig> {
  const list: EngineConfig[] = [
    { id: 'google', type: 'google', name: 'Google Translate', builtin: true, ...MT_LIMITS, maxConcurrent: 6, maxRps: 8 },
    { id: 'microsoft', type: 'microsoft', name: 'Microsoft Translator', apiKey: '', region: '', builtin: true, ...MT_LIMITS, batchSize: 50, batchChars: 9000 },
    { id: 'deepl', type: 'deepl', name: 'DeepL', apiKey: '', builtin: true, ...MT_LIMITS, maxConcurrent: 3, maxRps: 3, batchChars: 9000 },
    { id: 'deeplx', type: 'deeplx', name: 'DeepLX', baseUrl: '', builtin: true, ...MT_LIMITS, maxConcurrent: 2, maxRps: 3, batchSize: 1 },
    openaiPreset('openai', 'openai', 'OpenAI', { temperature: undefined }),
    {
      id: 'anthropic',
      type: 'anthropic',
      name: 'Claude',
      apiKey: '',
      model: DEFAULT_CLAUDE_MODEL,
      effort: 'low',
      builtin: true,
      ...LLM_LIMITS,
    },
    { id: 'gemini', type: 'gemini', name: 'Gemini', apiKey: '', model: 'gemini-flash-latest', builtin: true, ...LLM_LIMITS },
    openaiPreset('deepseek', 'deepseek'),
    openaiPreset('openrouter', 'openrouter'),
    openaiPreset('ollama', 'ollama', 'Ollama', { maxConcurrent: 1, maxRps: 2, batchSize: 8 }),
  ];
  return Object.fromEntries(list.map((e) => [e.id, e]));
}

function browserLanguage(): string {
  try {
    const nav = (globalThis as any).navigator;
    const lang = normalizeLang(nav?.language ?? 'en');
    return getLanguage(lang) ? lang : 'en';
  } catch {
    return 'en';
  }
}

export function createDefaultSettings(): Settings {
  const native = browserLanguage();
  return {
    version: SETTINGS_VERSION,
    uiLanguage: 'auto',
    theme: 'system',
    sourceLanguage: 'auto',
    targetLanguage: native === 'en' ? 'tr' : native,
    engineId: 'google',
    expertId: 'general',
    displayMode: 'bilingual',
    style: { id: 'underline', color: '#22c55e', textColor: '', fontScale: 1 },
    rules: {
      alwaysSites: [],
      neverSites: [],
      alwaysLanguages: [],
      neverLanguages: [],
      custom: [],
    },
    page: {
      autoTranslateAll: false,
      scope: 'main',
      lazy: true,
      skipCode: true,
      inlineShort: true,
      contextAware: true,
      translateIframes: true,
    },
    academic: {
      autoExpert: true,
      expertId: 'academic',
      skipReferences: true,
      keepTerms: true,
      termConsistency: true,
    },
    floatingBall: { enabled: true, side: 'right', top: 62, hiddenSites: [] },
    hover: { key: 'Shift' },
    input: { trigger: 'triple-space', targetLanguage: 'en', smartReverse: true, engineId: '', expertId: '' },
    selection: { mode: 'icon', display: 'inline', engineId: '', inInputs: false },
    engines: defaultEngines(),
    experts: [],
    expertOverrides: {},
    glossary: [],
    cache: { enabled: true, ttlDays: 30 },
    pdf: { autoOpen: false },
  };
}

function isPlainObject(v: unknown): v is Record<string, any> {
  return !!v && typeof v === 'object' && !Array.isArray(v);
}

/** Varsayılanların üzerine kayıtlı ayarları birleştirir (diziler değiştirilir, nesneler birleştirilir). */
export function mergeSettings(defaults: Settings, stored: unknown): Settings {
  if (!isPlainObject(stored)) return defaults;
  const merge = (base: any, over: any): any => {
    if (!isPlainObject(base) || !isPlainObject(over)) return over === undefined ? base : over;
    const out: Record<string, any> = { ...base };
    for (const [k, v] of Object.entries(over)) {
      if (v === undefined) continue;
      out[k] = isPlainObject(base[k]) && isPlainObject(v) ? merge(base[k], v) : v;
    }
    return out;
  };
  const merged = merge(defaults, stored) as Settings;
  // Motorlar: varsayılan alanlar yeni sürümlerde eklendiğinde kayıtlı motorlara da gelsin.
  const engines: Record<string, EngineConfig> = {};
  for (const [id, cfg] of Object.entries({ ...defaults.engines, ...(merged.engines ?? {}) })) {
    const base = defaults.engines[id];
    engines[id] = base ? { ...base, ...cfg, builtin: true } : (cfg as EngineConfig);
  }
  merged.engines = engines;
  if (!merged.engines[merged.engineId]) merged.engineId = 'google';
  return merged;
}

import type { EngineConfig, EngineType } from '../types';
import type { EngineMeta } from './types';

export const ENGINE_META: Record<EngineType, EngineMeta> = {
  google: {
    type: 'google',
    kind: 'mt',
    label: 'Google Translate',
    free: true,
    requiresKey: false,
    richText: 'html',
    fields: [],
  },
  microsoft: {
    type: 'microsoft',
    kind: 'mt',
    label: 'Microsoft Translator (Azure)',
    requiresKey: true,
    richText: 'html',
    fields: ['apiKey', 'region', 'baseUrl'],
    keyUrl: 'https://portal.azure.com/#create/Microsoft.CognitiveServicesTextTranslation',
    docsUrl: 'https://learn.microsoft.com/azure/ai-services/translator/text-translation/overview',
  },
  deepl: {
    type: 'deepl',
    kind: 'mt',
    label: 'DeepL',
    requiresKey: true,
    richText: 'html',
    fields: ['apiKey'],
    keyUrl: 'https://www.deepl.com/your-account/keys',
    docsUrl: 'https://developers.deepl.com/docs',
  },
  deeplx: {
    type: 'deeplx',
    kind: 'mt',
    label: 'DeepLX (kendi sunucunuz)',
    requiresKey: false,
    richText: 'none',
    fields: ['baseUrl', 'apiKey'],
    docsUrl: 'https://github.com/OwO-Network/DeepLX',
  },
  openai: {
    type: 'openai',
    kind: 'llm',
    label: 'OpenAI uyumlu',
    requiresKey: true,
    richText: 'html',
    fields: ['apiKey', 'baseUrl', 'model', 'temperature', 'extraBody', 'customHeaders'],
  },
  anthropic: {
    type: 'anthropic',
    kind: 'llm',
    label: 'Anthropic Claude',
    requiresKey: true,
    richText: 'html',
    fields: ['apiKey', 'model', 'effort', 'baseUrl'],
    keyUrl: 'https://console.anthropic.com/settings/keys',
    docsUrl: 'https://platform.claude.com/docs',
    models: ['claude-opus-5-5', 'claude-sonnet-5-5', 'claude-haiku-5-5', 'claude-fable-5-1'],
  },
  gemini: {
    type: 'gemini',
    kind: 'llm',
    label: 'Google Gemini',
    requiresKey: true,
    richText: 'html',
    fields: ['apiKey', 'model', 'temperature', 'extraBody'],
    keyUrl: 'https://aistudio.google.com/apikey',
    docsUrl: 'https://ai.google.dev/gemini-api/docs',
    models: ['gemini-flash-latest', 'gemini-flash-lite-latest', 'gemini-pro-latest', 'gemini-2.5-flash', 'gemini-2.5-pro'],
  },
};

export interface OpenAIPreset {
  id: string;
  label: string;
  baseUrl: string;
  model: string;
  models?: string[];
  keyUrl?: string;
  requiresKey?: boolean;
  headers?: Record<string, string>;
  note?: string;
}

/** OpenAI uyumlu API sunan servisler için hazır ayarlar. */
export const OPENAI_PRESETS: OpenAIPreset[] = [
  {
    id: 'openai',
    label: 'OpenAI',
    baseUrl: 'https://api.openai.com/v1',
    model: 'gpt-4.1-mini',
    models: ['gpt-4.1-mini', 'gpt-4.1', 'gpt-4.1-nano', 'gpt-4o-mini', 'gpt-5-mini', 'gpt-5-nano', 'gpt-5'],
    keyUrl: 'https://platform.openai.com/api-keys',
  },
  {
    id: 'deepseek',
    label: 'DeepSeek',
    baseUrl: 'https://api.deepseek.com',
    model: 'deepseek-chat',
    models: ['deepseek-chat', 'deepseek-reasoner'],
    keyUrl: 'https://platform.deepseek.com/api_keys',
  },
  {
    id: 'openrouter',
    label: 'OpenRouter',
    baseUrl: 'https://openrouter.ai/api/v1',
    model: 'openai/gpt-4.1-mini',
    models: ['openai/gpt-4.1-mini', 'google/gemini-2.5-flash', 'deepseek/deepseek-chat', 'qwen/qwen-2.5-72b-instruct', 'meta-llama/llama-3.3-70b-instruct'],
    keyUrl: 'https://openrouter.ai/keys',
    headers: { 'HTTP-Referer': 'https://github.com/asametsoyl/bioparrot', 'X-Title': 'BioParrot' },
  },
  {
    id: 'groq',
    label: 'Groq',
    baseUrl: 'https://api.groq.com/openai/v1',
    model: 'llama-3.3-70b-versatile',
    keyUrl: 'https://console.groq.com/keys',
  },
  {
    id: 'mistral',
    label: 'Mistral',
    baseUrl: 'https://api.mistral.ai/v1',
    model: 'mistral-small-latest',
    models: ['mistral-small-latest', 'mistral-medium-latest', 'mistral-large-latest'],
    keyUrl: 'https://console.mistral.ai/api-keys',
  },
  {
    id: 'xai',
    label: 'xAI Grok',
    baseUrl: 'https://api.x.ai/v1',
    model: 'grok-3-mini',
    keyUrl: 'https://console.x.ai',
  },
  {
    id: 'qwen',
    label: 'Qwen (DashScope)',
    baseUrl: 'https://dashscope-intl.aliyuncs.com/compatible-mode/v1',
    model: 'qwen-plus',
    models: ['qwen-plus', 'qwen-turbo', 'qwen-max', 'qwen-mt-turbo'],
    keyUrl: 'https://modelstudio.console.alibabacloud.com',
  },
  {
    id: 'siliconflow',
    label: 'SiliconFlow',
    baseUrl: 'https://api.siliconflow.com/v1',
    model: 'Qwen/Qwen2.5-72B-Instruct',
    keyUrl: 'https://cloud.siliconflow.com/account/ak',
  },
  {
    id: 'together',
    label: 'Together AI',
    baseUrl: 'https://api.together.xyz/v1',
    model: 'meta-llama/Llama-3.3-70B-Instruct-Turbo',
    keyUrl: 'https://api.together.ai/settings/api-keys',
  },
  {
    id: 'ollama',
    label: 'Ollama (yerel)',
    baseUrl: 'http://localhost:11434/v1',
    model: 'qwen2.5:7b',
    requiresKey: false,
    note: 'Ollama eklenti isteklerini reddederse OLLAMA_ORIGINS="*" ortam değişkeniyle başlatın.',
  },
  {
    id: 'lmstudio',
    label: 'LM Studio (yerel)',
    baseUrl: 'http://localhost:1234/v1',
    model: 'local-model',
    requiresKey: false,
  },
  {
    id: 'azure-openai',
    label: 'Azure OpenAI',
    baseUrl: 'https://KAYNAK.openai.azure.com/openai/deployments/DAGITIM/chat/completions?api-version=2024-10-21',
    model: 'gpt-4.1-mini',
    note: 'Tam dağıtım adresini yazın; anahtar "api-key" başlığıyla gönderilir.',
  },
  {
    id: 'custom',
    label: 'Özel (OpenAI uyumlu)',
    baseUrl: 'https://',
    model: '',
  },
];

export function getPreset(id: string | undefined): OpenAIPreset | undefined {
  return OPENAI_PRESETS.find((p) => p.id === id);
}

export function engineRequiresKey(config: EngineConfig): boolean {
  if (config.type === 'openai') {
    const preset = getPreset(config.preset);
    if (preset?.requiresKey === false) return false;
    if (/localhost|127\.0\.0\.1|\[::1\]/.test(config.baseUrl ?? '')) return false;
    return true;
  }
  return ENGINE_META[config.type].requiresKey;
}

/** Motorun kullanıma hazır olup olmadığı (anahtar vb. girilmiş mi). */
export function isEngineReady(config: EngineConfig | undefined): boolean {
  if (!config) return false;
  if (config.type === 'deeplx') return !!config.baseUrl?.trim();
  if (engineRequiresKey(config) && !config.apiKey?.trim()) return false;
  if (config.type === 'openai' && (!config.baseUrl?.trim() || !config.model?.trim())) return false;
  return true;
}

export function isLlmEngine(config: EngineConfig | undefined): boolean {
  return !!config && ENGINE_META[config.type].kind === 'llm';
}

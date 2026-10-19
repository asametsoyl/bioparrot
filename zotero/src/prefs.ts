import { createDefaultSettings, mergeSettings } from '@/core/defaults';
import type { Settings } from '@/core/types';

declare const Zotero: any;

export const PREF = 'extensions.bioparrot.';

export function getPref<T = string | boolean>(key: string, fallback?: T): T {
  const v = Zotero.Prefs.get(PREF + key, true);
  return (v === undefined || v === null ? fallback : v) as T;
}

export function setPref(key: string, value: string | boolean | number): void {
  Zotero.Prefs.set(PREF + key, value, true);
}

/**
 * Zotero tercihlerinden BioParrot ayarlarını oluşturur:
 * önce tarayıcı eklentisinden içe aktarılan tam ayarlar (JSON), sonra panelde girilen tekil tercihler.
 */
export function loadSettings(): Settings {
  let imported: unknown = {};
  const json = getPref<string>('settingsJson', '');
  if (json) {
    try {
      const parsed = JSON.parse(json);
      imported = parsed?.settings ?? parsed;
    } catch {
      /* geçersiz JSON yok sayılır */
    }
  }
  const s = mergeSettings(createDefaultSettings(), imported);
  const target = getPref<string>('targetLanguage', '');
  if (target) s.targetLanguage = target;
  const engine = getPref<string>('engineId', '');
  if (engine && s.engines[engine]) s.engineId = engine;
  const expert = getPref<string>('expertId', '');
  if (expert) {
    s.expertId = expert;
    s.academic.expertId = expert;
  }
  const keys: [string, string][] = [
    ['openaiKey', 'openai'],
    ['anthropicKey', 'anthropic'],
    ['geminiKey', 'gemini'],
    ['deeplKey', 'deepl'],
    ['deepseekKey', 'deepseek'],
    ['openrouterKey', 'openrouter'],
  ];
  for (const [pref, id] of keys) {
    const key = getPref<string>(pref, '').trim();
    if (key && s.engines[id]) s.engines[id].apiKey = key;
  }
  const models: [string, string][] = [
    ['openaiModel', 'openai'],
    ['anthropicModel', 'anthropic'],
    ['geminiModel', 'gemini'],
  ];
  for (const [pref, id] of models) {
    const model = getPref<string>(pref, '').trim();
    if (model && s.engines[id]) s.engines[id].model = model;
  }
  // Zotero'da sayfa bağlamı yoktur; akademik uzman varsayılandır.
  if (s.expertId === 'general') s.expertId = s.academic.expertId;
  return s;
}

/** PDF penceresindeki seçimleri (dil, motor, uzman) tercihlere geri yazar. */
export function saveBasicSettings(s: Settings): void {
  setPref('targetLanguage', s.targetLanguage);
  setPref('engineId', s.engineId);
  setPref('expertId', s.expertId);
}

import { browser } from 'wxt/browser';
import { STORAGE_KEYS } from '@/core/constants';
import { createDefaultSettings, mergeSettings } from '@/core/defaults';
import type { Settings } from '@/core/types';

/** Ayarları tarayıcı deposundan okur (varsayılanlarla birleştirilmiş olarak). */
export async function getSettings(): Promise<Settings> {
  const data = await browser.storage.local.get(STORAGE_KEYS.settings);
  return mergeSettings(createDefaultSettings(), data[STORAGE_KEYS.settings]);
}

export async function saveSettings(settings: Settings): Promise<void> {
  await browser.storage.local.set({ [STORAGE_KEYS.settings]: settings });
}

/** Ayarları bir fonksiyonla günceller ve kaydeder. */
export async function updateSettings(mutator: (draft: Settings) => Settings | void): Promise<Settings> {
  const draft = structuredClone(await getSettings());
  const next = mutator(draft) ?? draft;
  await saveSettings(next);
  return next;
}

/** Ayarlar değiştiğinde çağrılır; aboneliği iptal eden fonksiyon döner. */
export function watchSettings(callback: (settings: Settings) => void): () => void {
  const listener = (changes: Record<string, { newValue?: unknown }>, area: string) => {
    if (area !== 'local' || !(STORAGE_KEYS.settings in changes)) return;
    callback(mergeSettings(createDefaultSettings(), changes[STORAGE_KEYS.settings].newValue));
  };
  browser.storage.onChanged.addListener(listener);
  return () => browser.storage.onChanged.removeListener(listener);
}

/** Alan adı listesindeki desenlerle eşleşme (example.com, *.example.com). */
export function siteMatches(host: string, patterns: string[]): boolean {
  const h = host.toLowerCase().replace(/^www\./, '');
  return patterns.some((raw) => {
    const p = raw.trim().toLowerCase().replace(/^https?:\/\//, '').replace(/\/.*$/, '').replace(/^www\./, '');
    if (!p) return false;
    if (p.startsWith('*.')) {
      const base = p.slice(2);
      return h === base || h.endsWith('.' + base);
    }
    return h === p || h.endsWith('.' + p);
  });
}

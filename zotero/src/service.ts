import { setFetch } from '@/core/http';
import { TranslationService, type CacheStore } from '@/core/service';
import type { TranslateRequest } from '@/core/types';
import { loadSettings } from './prefs';

/** Bellekte tutulan basit LRU önbellek. */
class MemoryCache implements CacheStore {
  private map = new Map<string, string>();
  constructor(private max = 4000) {}
  async getMany(keys: string[]) {
    return keys.map((k) => {
      const v = this.map.get(k);
      if (v !== undefined) {
        this.map.delete(k);
        this.map.set(k, v);
      }
      return v;
    });
  }
  async setMany(entries: [string, string][]) {
    for (const [k, v] of entries) {
      this.map.set(k, v);
      if (this.map.size > this.max) this.map.delete(this.map.keys().next().value!);
    }
  }
}

let service: TranslationService | null = null;

export function getService(fetchImpl?: typeof fetch): TranslationService {
  if (fetchImpl) setFetch((input, init) => fetchImpl(input, init));
  service ??= new TranslationService({ getSettings: loadSettings, cache: new MemoryCache() });
  return service;
}

export async function translateText(text: string, extra: Partial<TranslateRequest> = {}): Promise<string> {
  const s = loadSettings();
  const res = await getService().translate({
    texts: [text],
    from: s.sourceLanguage,
    to: s.targetLanguage,
    engineId: s.engineId,
    expertId: s.expertId,
    rich: false,
    ...extra,
  });
  return res.translations[0] ?? '';
}

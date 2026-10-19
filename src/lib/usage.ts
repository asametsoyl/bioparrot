import { browser } from 'wxt/browser';
import { STORAGE_KEYS } from '@/core/constants';
import type { EngineUsage } from '@/core/engines/types';
import type { UsageRecorder } from '@/core/service';

export interface UsageStat {
  chars: number;
  requests: number;
  input: number;
  output: number;
}

/** Gün → motor → istatistik */
export type UsageData = Record<string, Record<string, UsageStat>>;

const KEEP_DAYS = 60;

function today(): string {
  return new Date().toISOString().slice(0, 10);
}

/** Kullanım istatistiklerini bellekte toplayıp birkaç saniyede bir kaydeder. */
export class StorageUsageRecorder implements UsageRecorder {
  private pending: Record<string, UsageStat> = {};
  private timer: ReturnType<typeof setTimeout> | null = null;

  record(engineId: string, chars: number, usage?: EngineUsage): void {
    const s = (this.pending[engineId] ??= { chars: 0, requests: 0, input: 0, output: 0 });
    s.chars += chars;
    s.requests += 1;
    s.input += usage?.input ?? 0;
    s.output += usage?.output ?? 0;
    if (!this.timer) this.timer = setTimeout(() => void this.flush(), 4000);
  }

  async flush(): Promise<void> {
    this.timer = null;
    const pending = this.pending;
    this.pending = {};
    if (!Object.keys(pending).length) return;
    const data = ((await browser.storage.local.get(STORAGE_KEYS.usage))[STORAGE_KEYS.usage] ?? {}) as UsageData;
    const day = (data[today()] ??= {});
    for (const [id, s] of Object.entries(pending)) {
      const d = (day[id] ??= { chars: 0, requests: 0, input: 0, output: 0 });
      d.chars += s.chars;
      d.requests += s.requests;
      d.input += s.input;
      d.output += s.output;
    }
    const cutoff = new Date(Date.now() - KEEP_DAYS * 86_400_000).toISOString().slice(0, 10);
    for (const key of Object.keys(data)) if (key < cutoff) delete data[key];
    await browser.storage.local.set({ [STORAGE_KEYS.usage]: data });
  }
}

export async function readUsage(): Promise<UsageData> {
  return ((await browser.storage.local.get(STORAGE_KEYS.usage))[STORAGE_KEYS.usage] ?? {}) as UsageData;
}

import { afterEach, describe, expect, it, vi } from 'vitest';
import { createDefaultSettings } from '@/core/defaults';
import { setFetch } from '@/core/http';
import { TranslationService, type CacheStore } from '@/core/service';
import type { Settings } from '@/core/types';

class MapCache implements CacheStore {
  map = new Map<string, string>();
  async getMany(keys: string[]) {
    return keys.map((k) => this.map.get(k));
  }
  async setMany(entries: [string, string][]) {
    for (const [k, v] of entries) this.map.set(k, v);
  }
}

function json(data: unknown, status = 200) {
  return new Response(JSON.stringify(data), { status, headers: { 'Content-Type': 'application/json' } });
}

afterEach(() => setFetch(null));

describe('TranslationService', () => {
  it('translates with Google, caches results and reports detected language', async () => {
    const fetchMock = vi.fn(async (_url: string, init?: RequestInit) => {
      const body = new URLSearchParams(String(init?.body));
      return json(body.getAll('q').map((q) => [`TR(${q})`, 'en']));
    });
    setFetch(fetchMock as any);
    const settings = createDefaultSettings();
    const cache = new MapCache();
    const service = new TranslationService({ getSettings: () => settings, cache });
    const first = await service.translate({ texts: ['Hello', 'World'], from: 'auto', to: 'tr', engineId: 'google', rich: true });
    expect(first.translations).toEqual(['TR(Hello)', 'TR(World)']);
    expect(first.detected).toEqual(['en', 'en']);
    await new Promise((r) => setTimeout(r, 10));
    const second = await service.translate({ texts: ['Hello'], from: 'auto', to: 'tr', engineId: 'google', rich: true });
    expect(second.translations).toEqual(['TR(Hello)']);
    expect(fetchMock).toHaveBeenCalledTimes(1);
  });

  it('protects glossary terms for MT engines', async () => {
    let sent = '';
    setFetch((async (_url: string, init?: RequestInit) => {
      const q = new URLSearchParams(String(init?.body)).getAll('q');
      sent = q[0];
      return json(q.map((x) => [x.replace('Each', 'Her').replace('carries RNA.', 'RNA taşır.'), 'en']));
    }) as any);
    const settings: Settings = { ...createDefaultSettings(), glossary: [{ id: '1', source: 'extracellular vesicle', target: 'hücre dışı vezikül' }] };
    const service = new TranslationService({ getSettings: () => settings });
    const res = await service.translate({ texts: ['Each extracellular vesicle carries RNA.'], from: 'auto', to: 'tr', engineId: 'google', rich: true, noCache: true });
    expect(sent).toContain('translate="no"');
    expect(res.translations[0]).toBe('Her hücre dışı vezikül RNA taşır.');
  });

  it('splits LLM batches and retries segments the model skipped', async () => {
    const calls: string[] = [];
    setFetch((async (_url: string, init?: RequestInit) => {
      const body = JSON.parse(String(init?.body));
      const user: string = body.messages[1].content;
      calls.push(user);
      if (user.includes('<seg id="2">')) return json({ choices: [{ message: { content: '<seg id="1">Bir</seg>' } }] });
      return json({ choices: [{ message: { content: 'İki' } }], usage: { prompt_tokens: 10, completion_tokens: 2 } });
    }) as any);
    const settings = createDefaultSettings();
    settings.engines.openai.apiKey = 'sk-test';
    const service = new TranslationService({ getSettings: () => settings });
    const res = await service.translate({ texts: ['One', 'Two'], from: 'en', to: 'tr', engineId: 'openai', rich: false, noCache: true });
    expect(res.translations).toEqual(['Bir', 'İki']);
    expect(calls).toHaveLength(2);
  });

  it('requires an API key for keyed engines', async () => {
    const service = new TranslationService({ getSettings: () => createDefaultSettings() });
    await expect(service.translate({ texts: ['x'], from: 'auto', to: 'tr', engineId: 'deepl' })).rejects.toThrow(/API/);
  });
});

import { getEngine } from './engines';
import { ENGINE_META, engineRequiresKey, isLlmEngine } from './engines/meta';
import type { EngineUsage } from './engines/types';
import { getExpert } from './experts';
import { hashString } from './hash';
import { makeBatches } from './batch';
import { Limiter, sleep } from './limiter';
import { applyGlossaryMarkup, markupToPlain, relevantGlossary, stripGlossaryMarkup } from './markup';
import { buildExplainPrompt, buildTermExtractionPrompt, parseTermList } from './prompt';
import {
  TranslateError,
  type EngineConfig,
  type Settings,
  type TranslateRequest,
  type TranslateResponse,
  type TranslationContext,
} from './types';

export interface CacheStore {
  getMany(keys: string[]): Promise<(string | undefined)[]>;
  setMany(entries: [string, string][]): Promise<void>;
}

export interface UsageRecorder {
  record(engineId: string, chars: number, usage?: EngineUsage): void;
}

export interface ServiceDeps {
  getSettings: () => Settings | Promise<Settings>;
  cache?: CacheStore;
  usage?: UsageRecorder;
}

const MAX_RETRIES = 2;

export class TranslationService {
  private limiters = new Map<string, Limiter>();

  constructor(private deps: ServiceDeps) {}

  private limiterFor(config: EngineConfig): Limiter {
    let l = this.limiters.get(config.id);
    if (!l) {
      l = new Limiter(config.maxConcurrent, config.maxRps);
      this.limiters.set(config.id, l);
    } else {
      l.configure(config.maxConcurrent, config.maxRps);
    }
    return l;
  }

  resolveEngine(settings: Settings, engineId?: string): EngineConfig {
    const config = settings.engines[engineId || settings.engineId] ?? settings.engines.google;
    if (!config) throw new TranslateError('Çeviri motoru bulunamadı.', 'config');
    if (config.type === 'deeplx' && !config.baseUrl?.trim()) {
      throw new TranslateError(`${config.name}: sunucu adresi girilmemiş. Ayarlar → Çeviri motorları.`, 'config');
    }
    if (engineRequiresKey(config) && !config.apiKey?.trim()) {
      throw new TranslateError(`${config.name}: API anahtarı girilmemiş. Ayarlar → Çeviri motorları.`, 'config');
    }
    return config;
  }

  private cacheKey(config: EngineConfig, req: TranslateRequest, expertKey: string, glossaryKey: string, text: string): string {
    const engineKey = `${config.type}|${config.model ?? ''}|${config.type === 'openai' || config.type === 'deeplx' ? config.baseUrl ?? '' : ''}`;
    return hashString(`${engineKey}|${req.from}|${req.to}|${expertKey}|${glossaryKey}|${text}`) + ':' + text.length;
  }

  async translate(req: TranslateRequest, signal?: AbortSignal): Promise<TranslateResponse> {
    const settings = await this.deps.getSettings();
    const config = this.resolveEngine(settings, req.engineId);
    const meta = ENGINE_META[config.type];
    const llm = meta.kind === 'llm';
    const expert = getExpert(settings, req.expertId || settings.expertId);
    const rich = !!req.rich && meta.richText === 'html';

    // Motor zengin metni desteklemiyorsa etiketleri kaldır.
    const texts = req.texts.map((t) => (req.rich && !rich ? markupToPlain(t) : t));
    const out: (string | undefined)[] = new Array(texts.length).fill(undefined);
    const detected: (string | undefined)[] = new Array(texts.length).fill(undefined);

    const glossary = relevantGlossary(texts, settings.glossary);
    const glossaryKey = glossary.length ? hashString(glossary.map((g) => `${g.source}=${g.target}`).join('|')) : '';
    const expertKey = llm ? `${expert.id}:${hashString(expert.systemPrompt)}:${settings.academic.keepTerms ? 1 : 0}` : '';

    // 1) Önbellek
    const useCache = settings.cache.enabled && !req.noCache && !!this.deps.cache;
    const keys = texts.map((t) => this.cacheKey(config, req, expertKey, glossaryKey, t));
    if (useCache) {
      try {
        const hits = await this.deps.cache!.getMany(keys);
        hits.forEach((v, i) => {
          if (v !== undefined) out[i] = v;
        });
      } catch {
        /* önbellek hatası çeviriyi engellemesin */
      }
    }

    // 2) Önbellekte olmayanları gruplara ayır
    const pending = texts.map((_, i) => i).filter((i) => out[i] === undefined && texts[i].trim() !== '');
    texts.forEach((t, i) => {
      if (t.trim() === '') out[i] = t;
    });
    const batches = makeBatches(pending, texts, config.batchSize, config.batchChars);

    const runBatch = async (indices: number[]) => {
      const inputs = indices.map((i) => texts[i]);
      const result = await this.callEngine(config, {
        texts: llm ? inputs : inputs.map((t) => (rich ? applyGlossaryMarkup(t, glossary) : t)),
        from: req.from,
        to: req.to,
        rich,
        llm,
        expert,
        context: req.context,
        glossary: llm ? relevantGlossary(inputs, glossary) : [],
        academic: settings.academic,
        signal,
      });
      result.translations.forEach((t, j) => {
        if (t !== undefined) out[indices[j]] = llm ? t : stripGlossaryMarkup(t);
        detected[indices[j]] = result.detected?.[j];
      });
      // Model bazı parçaları atladıysa onları tek tek çevir.
      const missing = indices.filter((i) => out[i] === undefined);
      if (missing.length && indices.length > 1) {
        await Promise.all(
          missing.map(async (i) => {
            const single = await this.callEngine(config, {
              texts: [llm ? texts[i] : rich ? applyGlossaryMarkup(texts[i], glossary) : texts[i]],
              from: req.from,
              to: req.to,
              rich,
              llm,
              expert,
              context: req.context,
              glossary: llm ? relevantGlossary([texts[i]], glossary) : [],
              academic: settings.academic,
              signal,
            });
            const t = single.translations[0];
            if (t !== undefined) out[i] = llm ? t : stripGlossaryMarkup(t);
          }),
        );
      }
    };

    await Promise.all(batches.map(runBatch));

    const failed = pending.filter((i) => out[i] === undefined);
    if (failed.length) {
      throw new TranslateError(`${failed.length} parça çevrilemedi (model yanıtı eksik).`, 'parse');
    }

    // 3) Önbelleğe yaz
    if (useCache && pending.length) {
      const entries: [string, string][] = pending.map((i) => [keys[i], out[i]!]);
      this.deps.cache!.setMany(entries).catch(() => {});
    }

    return { translations: out as string[], detected, engineId: config.id, rich };
  }

  private async callEngine(
    config: EngineConfig,
    args: {
      texts: string[];
      from: string;
      to: string;
      rich: boolean;
      llm: boolean;
      expert: ReturnType<typeof getExpert>;
      context?: TranslationContext;
      glossary: Settings['glossary'];
      academic: Settings['academic'];
      signal?: AbortSignal;
    },
  ) {
    const engine = getEngine(config.type);
    const limiter = this.limiterFor(config);
    const chars = args.texts.reduce((n, t) => n + t.length, 0);
    let attempt = 0;
    for (;;) {
      try {
        const result = await limiter.run(() =>
          engine.translate(
            {
              texts: args.texts,
              from: args.from,
              to: args.to,
              rich: args.rich,
              expert: args.expert,
              context: args.context,
              glossary: args.glossary,
              academic: args.academic,
              signal: args.signal,
            },
            config,
          ),
        );
        this.deps.usage?.record(config.id, chars, result.usage);
        return result;
      } catch (err) {
        const e = err instanceof TranslateError ? err : new TranslateError(String((err as Error)?.message ?? err), 'unknown');
        const retryable = e.code === 'quota' || e.code === 'server' || e.code === 'network';
        // Kalıcı kota hataları (402/456) yeniden denenmez.
        if (!retryable || e.status === 402 || e.status === 456 || attempt >= MAX_RETRIES || args.signal?.aborted) throw e;
        attempt++;
        await sleep(e.retryAfterMs ?? 800 * 2 ** attempt + Math.random() * 400, args.signal);
      }
    }
  }

  /** Bir metinden belgeye özel tutarlı terim listesi çıkarır (yalnızca yapay zeka motorları). */
  async extractTerms(text: string, to: string, engineId?: string, expertId?: string): Promise<{ source: string; target: string }[]> {
    const settings = await this.deps.getSettings();
    const config = this.resolveEngine(settings, engineId);
    if (!isLlmEngine(config)) return [];
    const engine = getEngine(config.type);
    if (!engine.complete) return [];
    const expert = getExpert(settings, expertId || settings.academic.expertId);
    const prompt = buildTermExtractionPrompt(text.slice(0, 6000), to, expert);
    const res = await this.limiterFor(config).run(() => engine.complete!(prompt.system, prompt.user, config));
    this.deps.usage?.record(config.id, text.length, res.usage);
    return parseTermList(res.text);
  }

  /** Seçili metni açıklar (yapay zeka motoru gerekir). */
  async explain(text: string, to: string, context?: TranslationContext, engineId?: string, expertId?: string): Promise<string> {
    const settings = await this.deps.getSettings();
    let config = this.resolveEngine(settings, engineId);
    if (!isLlmEngine(config)) {
      const fallback = Object.values(settings.engines).find((e) => isLlmEngine(e) && !!e.apiKey?.trim());
      if (!fallback) throw new TranslateError('Açıklama için ayarlardan bir yapay zeka motoru (OpenAI, Claude, Gemini…) yapılandırın.', 'config');
      config = fallback;
    }
    const engine = getEngine(config.type);
    const expert = getExpert(settings, expertId || settings.expertId);
    const prompt = buildExplainPrompt(text.slice(0, 4000), to, expert, context);
    const res = await this.limiterFor(config).run(() => engine.complete!(prompt.system, prompt.user, config));
    this.deps.usage?.record(config.id, text.length, res.usage);
    return res.text.replace(/<think(?:ing)?>[\s\S]*?<\/think(?:ing)?>/gi, '').trim();
  }
}

import Anthropic from '@anthropic-ai/sdk';
import { getFetch } from '../http';
import { TranslateError, type EngineConfig } from '../types';
import { llmTranslate } from './llm';
import type { Engine, EngineUsage } from './types';

export const DEFAULT_CLAUDE_MODEL = 'claude-opus-5-5';

/** Sunucu tarafı yedek model (fallbacks: "default") destekleyen modeller. */
const FALLBACK_MODELS = /^claude-(fable-5-1|opus-5-5|opus-5|sonnet-5-5)$/;

/** output_config.effort destekleyen modeller (Claude 5 ailesi ve Opus/Sonnet 4.6+). */
function supportsEffort(model: string): boolean {
  return /^claude-(fable|mythos|opus|sonnet|haiku)-5/.test(model) || /^claude-(opus|sonnet)-4-[678]/.test(model);
}

const clients = new Map<string, Anthropic>();

function getClient(config: EngineConfig): Anthropic {
  const apiKey = config.apiKey?.trim();
  if (!apiKey) throw new TranslateError('Claude için API anahtarı girin.', 'config');
  const baseURL = config.baseUrl?.trim() || undefined;
  const cacheKey = `${apiKey}|${baseURL ?? ''}`;
  let client = clients.get(cacheKey);
  if (!client) {
    client = new Anthropic({
      apiKey,
      baseURL,
      // Anahtar kullanıcının kendi tarayıcısında saklanır ve yalnızca Anthropic'e gönderilir.
      dangerouslyAllowBrowser: true,
      maxRetries: 2,
      timeout: 180_000,
      fetch: getFetch() as unknown as typeof fetch,
    });
    clients.set(cacheKey, client);
  }
  return client;
}

function mapError(err: unknown): unknown {
  if (err instanceof TranslateError) return err;
  if (err instanceof Anthropic.AuthenticationError || err instanceof Anthropic.PermissionDeniedError) {
    return new TranslateError(`Claude yetkilendirme hatası: ${err.message}`, 'auth', err.status);
  }
  if (err instanceof Anthropic.RateLimitError) {
    return new TranslateError(`Claude istek sınırı aşıldı: ${err.message}`, 'quota', err.status);
  }
  if (err instanceof Anthropic.NotFoundError) {
    return new TranslateError(`Claude modeli bulunamadı: ${err.message}`, 'config', err.status);
  }
  if (err instanceof Anthropic.BadRequestError) {
    return new TranslateError(`Claude isteği geçersiz: ${err.message}`, 'config', err.status);
  }
  if (err instanceof Anthropic.APIConnectionError) {
    return new TranslateError(`Claude'a bağlanılamadı: ${err.message}`, 'network');
  }
  if (err instanceof Anthropic.InternalServerError) {
    return new TranslateError(`Claude sunucu hatası (${err.status}): ${err.message}`, 'server', err.status);
  }
  if (err instanceof Anthropic.APIError) {
    return new TranslateError(`Claude hatası (${err.status ?? '?'}): ${err.message}`, 'unknown', err.status);
  }
  return err;
}

async function complete(
  config: EngineConfig,
  stableSystem: string,
  contextBlock: string,
  user: string,
  signal?: AbortSignal,
): Promise<{ text: string; usage?: EngineUsage }> {
  const client = getClient(config);
  const model = config.model?.trim() || DEFAULT_CLAUDE_MODEL;
  const official = !config.baseUrl?.trim() || /api\.anthropic\.com/.test(config.baseUrl);

  // Sabit sistem istemi önbelleğe alınır; sayfaya özel bağlam ondan sonra gelir.
  const system: Anthropic.Beta.BetaTextBlockParam[] = [
    { type: 'text', text: stableSystem, cache_control: { type: 'ephemeral' } },
  ];
  if (contextBlock) system.push({ type: 'text', text: contextBlock });

  const params: Anthropic.Beta.MessageCreateParamsNonStreaming = {
    model,
    max_tokens: 16000,
    system,
    messages: [{ role: 'user', content: user }],
  };
  if (config.effort && supportsEffort(model)) params.output_config = { effort: config.effort };
  if (official && FALLBACK_MODELS.test(model)) {
    // Güvenlik sınıflandırıcısı reddederse istek sunucu tarafında önerilen yedek modelle yeniden çalışır.
    params.betas = ['server-side-fallback-2026-07-01'];
    params.fallbacks = 'default';
  }

  let response: Anthropic.Beta.BetaMessage;
  try {
    response = await client.beta.messages.create(params, { signal });
  } catch (err) {
    throw mapError(err);
  }
  if (response.stop_reason === 'refusal') {
    throw new TranslateError('Claude bu metni çevirmeyi reddetti.', 'refused');
  }
  const text = response.content
    .filter((block): block is Anthropic.Beta.BetaTextBlock => block.type === 'text')
    .map((block) => block.text)
    .join('');
  const u = response.usage;
  return {
    text,
    usage: {
      input: (u.input_tokens ?? 0) + (u.cache_read_input_tokens ?? 0) + (u.cache_creation_input_tokens ?? 0),
      output: u.output_tokens,
    },
  };
}

/** Anthropic Claude (resmi SDK). */
export const anthropicEngine: Engine = {
  translate(input, config) {
    return llmTranslate(input, config, (prompt, cfg, signal) =>
      complete(cfg, prompt.stableSystem, prompt.contextBlock, prompt.user, signal),
    );
  },
  complete(system, user, config, signal) {
    return complete(config, system, '', user, signal);
  },
};

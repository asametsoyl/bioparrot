import { TranslateError } from './types';

export type FetchLike = (input: string, init?: RequestInit) => Promise<Response>;

let customFetch: FetchLike | null = null;

/** Zotero gibi ortamlarda farklı bir fetch uygulaması kullanmak için. */
export function setFetch(fn: FetchLike | null): void {
  customFetch = fn;
}

export function getFetch(): FetchLike {
  return customFetch ?? ((input, init) => fetch(input, init));
}

export interface HttpOptions extends Omit<RequestInit, 'body'> {
  body?: BodyInit | null;
  json?: unknown;
  timeoutMs?: number;
}

/** Bir HTTP isteği yapar; hata durumunda anlamlı bir TranslateError fırlatır. */
export async function http(url: string, opts: HttpOptions = {}): Promise<Response> {
  const { json, timeoutMs = 60_000, signal, headers, ...rest } = opts;
  const controller = new AbortController();
  const timer = setTimeout(() => controller.abort(new DOMException('timeout', 'TimeoutError')), timeoutMs);
  const onAbort = () => controller.abort(signal?.reason);
  signal?.addEventListener('abort', onAbort, { once: true });

  const finalHeaders = new Headers(headers);
  let body = rest.body;
  if (json !== undefined) {
    finalHeaders.set('Content-Type', 'application/json');
    body = JSON.stringify(json);
  }

  let res: Response;
  try {
    res = await getFetch()(url, { ...rest, body, headers: finalHeaders, signal: controller.signal });
  } catch (err) {
    const timedOut = controller.signal.reason instanceof DOMException && controller.signal.reason.name === 'TimeoutError';
    throw new TranslateError(
      timedOut ? 'İstek zaman aşımına uğradı.' : `Ağ hatası: ${(err as Error)?.message ?? err}`,
      'network',
    );
  } finally {
    clearTimeout(timer);
    signal?.removeEventListener('abort', onAbort);
  }

  if (!res.ok) throw await errorFromResponse(res);
  return res;
}

export async function httpJson<T = any>(url: string, opts: HttpOptions = {}): Promise<T> {
  const res = await http(url, opts);
  const text = await res.text();
  try {
    return JSON.parse(text) as T;
  } catch {
    throw new TranslateError(`Beklenmeyen yanıt: ${text.slice(0, 200)}`, 'parse', res.status);
  }
}

export async function errorFromResponse(res: Response): Promise<TranslateError> {
  let detail = '';
  try {
    const text = await res.text();
    try {
      const data = JSON.parse(text);
      detail =
        data?.error?.message ??
        data?.error?.code ??
        (typeof data?.error === 'string' ? data.error : '') ??
        '';
      if (!detail) detail = data?.message ?? data?.detail ?? data?.msg ?? '';
      if (typeof detail !== 'string') detail = JSON.stringify(detail);
    } catch {
      detail = text.replace(/<[^>]+>/g, ' ').replace(/\s+/g, ' ').trim().slice(0, 200);
    }
  } catch {
    /* gövde okunamadı */
  }
  return httpError(res.status, detail, res.headers.get('retry-after'));
}

export function httpError(status: number, detail: string, retryAfter?: string | null): TranslateError {
  const suffix = detail ? ` – ${detail}` : '';
  const retryAfterMs = parseRetryAfter(retryAfter);
  if (status === 401 || status === 403) return new TranslateError(`Yetkilendirme hatası (${status}): API anahtarını kontrol edin${suffix}`, 'auth', status);
  if (status === 402) return new TranslateError(`Bakiye/kota yetersiz (402)${suffix}`, 'quota', status);
  if (status === 429) return new TranslateError(`Çok fazla istek veya kota aşıldı (429)${suffix}`, 'quota', status, retryAfterMs);
  if (status === 456) return new TranslateError(`Kota doldu (456)${suffix}`, 'quota', status);
  if (status === 404) return new TranslateError(`Adres veya model bulunamadı (404)${suffix}`, 'config', status);
  if (status === 400 || status === 422) return new TranslateError(`Geçersiz istek (${status})${suffix}`, 'config', status);
  if (status >= 500) return new TranslateError(`Sunucu hatası (${status})${suffix}`, 'server', status, retryAfterMs);
  return new TranslateError(`HTTP ${status}${suffix}`, 'unknown', status);
}

function parseRetryAfter(value: string | null | undefined): number | undefined {
  if (!value) return undefined;
  const sec = Number(value);
  if (Number.isFinite(sec)) return Math.min(sec * 1000, 60_000);
  const date = Date.parse(value);
  if (Number.isFinite(date)) return Math.max(0, Math.min(date - Date.now(), 60_000));
  return undefined;
}

export function parseJsonObject(value: string | undefined, label: string): Record<string, any> {
  if (!value || !value.trim()) return {};
  try {
    const obj = JSON.parse(value);
    if (obj && typeof obj === 'object' && !Array.isArray(obj)) return obj;
  } catch {
    /* aşağıda hata */
  }
  throw new TranslateError(`${label} geçerli bir JSON nesnesi değil.`, 'config');
}

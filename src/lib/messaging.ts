import { browser } from 'wxt/browser';
import {
  TranslateError,
  type DisplayMode,
  type TranslateRequest,
  type TranslateResponse,
  type TranslationContext,
} from '@/core/types';

/** İçerik betiği / açılır pencere → arka plan mesajları. */
export type BackgroundMessage =
  | { type: 'translate'; req: TranslateRequest }
  | { type: 'extractTerms'; text: string; to: string; engineId?: string; expertId?: string }
  | { type: 'explain'; text: string; to: string; context?: TranslationContext; engineId?: string; expertId?: string }
  | { type: 'testEngine'; engineId: string; to: string; text: string }
  | { type: 'openOptions'; section?: string }
  | { type: 'openPdf'; url?: string }
  | { type: 'broadcast'; action: TabAction }
  | { type: 'cacheStats' }
  | { type: 'clearCache' };

export interface BackgroundResults {
  translate: TranslateResponse;
  extractTerms: { source: string; target: string }[];
  explain: string;
  testEngine: { translation: string; ms: number };
  openOptions: void;
  openPdf: void;
  broadcast: void;
  cacheStats: { entries: number };
  clearCache: void;
}

/** Arka plan / açılır pencere → sekmedeki içerik betiklerine gönderilen eylemler. */
export type TabAction =
  | { type: 'toggle' }
  | { type: 'translate' }
  | { type: 'restore' }
  | { type: 'toggleMode' }
  | { type: 'setMode'; mode: DisplayMode }
  | { type: 'translateSelection'; text?: string }
  | { type: 'toggleHover' }
  | { type: 'getState' }
  | { type: 'settingsHint' };

export interface PageState {
  active: boolean;
  mode: DisplayMode;
  pageLang?: string;
  isPaper: boolean;
  isPdf: boolean;
  pending: number;
  expertId?: string;
  url: string;
}

export interface TabActionMessage {
  type: 'bp:tab';
  action: TabAction;
}

type Result<T> = { ok: true; data: T } | { ok: false; error: string; code?: TranslateError['code'] };

export async function sendToBackground<K extends BackgroundMessage['type']>(
  msg: Extract<BackgroundMessage, { type: K }>,
): Promise<BackgroundResults[K]> {
  let res: Result<BackgroundResults[K]> | undefined;
  try {
    res = (await browser.runtime.sendMessage(msg)) as Result<BackgroundResults[K]> | undefined;
  } catch (err) {
    throw new TranslateError(String((err as Error)?.message ?? err), 'network');
  }
  if (!res) throw new TranslateError('Arka plan yanıt vermedi.', 'network');
  if (!res.ok) throw new TranslateError(res.error, res.code ?? 'unknown');
  return res.data;
}

export function okResult<T>(data: T): Result<T> {
  return { ok: true, data };
}

export function errorResult(err: unknown): Result<never> {
  if (err instanceof TranslateError) return { ok: false, error: err.message, code: err.code };
  return { ok: false, error: String((err as Error)?.message ?? err), code: 'unknown' };
}

/** Bir sekmenin tüm çerçevelerine (veya yalnızca ana çerçeveye) eylem gönderir. */
export async function sendTabAction<T = unknown>(tabId: number, action: TabAction, topFrameOnly = false): Promise<T | undefined> {
  const msg: TabActionMessage = { type: 'bp:tab', action };
  try {
    return (await browser.tabs.sendMessage(tabId, msg, topFrameOnly ? { frameId: 0 } : undefined)) as T;
  } catch {
    return undefined;
  }
}

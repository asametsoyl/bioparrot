import type { UiLocale, UiLocalePref } from '../types';
import { en } from './en';
import { tr, type MessageKey } from './tr';

export type { MessageKey };

const DICTS: Record<UiLocale, Record<MessageKey, string>> = { tr, en };

let current: UiLocale = detectLocale();

function detectLocale(): UiLocale {
  try {
    const nav = (globalThis as any).navigator;
    return String(nav?.language ?? '').toLowerCase().startsWith('tr') ? 'tr' : 'en';
  } catch {
    return 'en';
  }
}

export function resolveLocale(pref: UiLocalePref | undefined): UiLocale {
  return pref === 'tr' || pref === 'en' ? pref : detectLocale();
}

export function setLocale(pref: UiLocalePref | undefined): UiLocale {
  current = resolveLocale(pref);
  return current;
}

export function getLocale(): UiLocale {
  return current;
}

export function t(key: MessageKey, params?: Record<string, string | number>): string {
  let s = DICTS[current][key] ?? tr[key] ?? key;
  if (params) s = s.replace(/\{(\w+)\}/g, (m, k: string) => (params[k] !== undefined ? String(params[k]) : m));
  return s;
}

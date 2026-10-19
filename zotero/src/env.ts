/* Zotero eklenti kapsamında tarayıcı API'lerini (fetch, Headers, AbortController…) ana pencereden sağlar. */

declare const Zotero: any;

const NAMES = [
  'window',
  'document',
  'navigator',
  'fetch',
  'Headers',
  'Request',
  'Response',
  'AbortController',
  'AbortSignal',
  'URL',
  'URLSearchParams',
  'TextEncoder',
  'TextDecoder',
  'DOMException',
  'DOMParser',
  'Blob',
  'File',
  'FormData',
  'ReadableStream',
  'setTimeout',
  'clearTimeout',
  'setInterval',
  'clearInterval',
  'queueMicrotask',
  'structuredClone',
  'crypto',
  'performance',
  'atob',
  'btoa',
  'Event',
  'EventTarget',
  'CustomEvent',
] as const;

function mainWindow(): any {
  return Zotero.getMainWindow?.() ?? Services.wm.getMostRecentWindow('navigator:browser');
}

declare const Services: any;

/** Eksik genel nesneleri, ilk erişimde ana pencereden alınacak şekilde tanımlar. */
export function installGlobals(target: any): void {
  for (const name of NAMES) {
    for (const obj of [target, globalThis as any]) {
      if (!obj || name in obj) continue;
      try {
        Object.defineProperty(obj, name, {
          configurable: true,
          get() {
            const win = mainWindow();
            const value = win?.[name];
            return typeof value === 'function' && !/^[A-Z]/.test(name) ? value.bind(win) : value;
          },
        });
      } catch {
        /* salt okunur özellik */
      }
    }
  }
}

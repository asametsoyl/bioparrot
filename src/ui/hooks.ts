import { useCallback, useEffect, useState } from 'react';
import { setLocale } from '@/core/i18n';
import type { Settings } from '@/core/types';
import { getSettings, updateSettings, watchSettings } from '@/lib/settings';

export type UpdateSettings = (mutator: (draft: Settings) => void) => Promise<void>;

/** Ayarları yükler, değişiklikleri izler ve güncelleme fonksiyonu sağlar. */
export function useSettings(): { settings: Settings | null; update: UpdateSettings } {
  const [settings, setSettings] = useState<Settings | null>(null);

  useEffect(() => {
    let alive = true;
    void getSettings().then((s) => alive && setSettings(s));
    const unwatch = watchSettings((s) => setSettings(s));
    return () => {
      alive = false;
      unwatch();
    };
  }, []);

  const update = useCallback<UpdateSettings>(async (mutator) => {
    setSettings((prev) => {
      if (!prev) return prev;
      const draft = structuredClone(prev);
      mutator(draft);
      return draft;
    });
    await updateSettings((d) => {
      mutator(d);
    });
  }, []);

  if (settings) setLocale(settings.uiLanguage);
  return { settings, update };
}

/** Tema tercihini <html data-theme> özniteliğine uygular. */
export function useTheme(theme: Settings['theme'] | undefined): void {
  useEffect(() => {
    document.documentElement.dataset.theme = theme ?? 'system';
  }, [theme]);
}

export { useDebouncedCallback } from './use-debounced';

import { useCallback, useEffect, useRef } from 'react';

/** Bir değeri belirtilen süre sonra işleyen gecikmeli geri çağırma. */
export function useDebouncedCallback<T extends unknown[]>(fn: (...args: T) => void, ms: number): (...args: T) => void {
  const ref = useRef(fn);
  ref.current = fn;
  const timer = useRef<ReturnType<typeof setTimeout> | null>(null);
  useEffect(() => () => {
    if (timer.current) clearTimeout(timer.current);
  }, []);
  return useCallback(
    (...args: T) => {
      if (timer.current) clearTimeout(timer.current);
      timer.current = setTimeout(() => ref.current(...args), ms);
    },
    [ms],
  );
}

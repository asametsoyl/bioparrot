/** Motor başına eşzamanlılık ve saniyedeki istek sınırı uygulayan basit kuyruk. */
export class Limiter {
  private active = 0;
  private queue: (() => void)[] = [];
  private lastStart = 0;
  private timer: ReturnType<typeof setTimeout> | null = null;

  constructor(
    public maxConcurrent: number,
    public maxRps: number,
  ) {}

  configure(maxConcurrent: number, maxRps: number): void {
    this.maxConcurrent = Math.max(1, maxConcurrent || 1);
    this.maxRps = Math.max(0.1, maxRps || 1);
    this.pump();
  }

  run<T>(task: () => Promise<T>): Promise<T> {
    return new Promise<T>((resolve, reject) => {
      this.queue.push(() => {
        this.active++;
        this.lastStart = Date.now();
        task()
          .then(resolve, reject)
          .finally(() => {
            this.active--;
            this.pump();
          });
      });
      this.pump();
    });
  }

  private pump(): void {
    if (this.timer) return;
    while (this.queue.length && this.active < this.maxConcurrent) {
      const interval = 1000 / this.maxRps;
      const wait = this.lastStart + interval - Date.now();
      if (wait > 0) {
        this.timer = setTimeout(() => {
          this.timer = null;
          this.pump();
        }, wait);
        return;
      }
      this.queue.shift()!();
    }
  }
}

export function sleep(ms: number, signal?: AbortSignal): Promise<void> {
  return new Promise((resolve, reject) => {
    const t = setTimeout(resolve, ms);
    signal?.addEventListener(
      'abort',
      () => {
        clearTimeout(t);
        reject(signal.reason);
      },
      { once: true },
    );
  });
}

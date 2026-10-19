import type { CacheStore } from '@/core/service';

const DB_NAME = 'bioparrot-cache';
const STORE = 'translations';

interface Entry {
  v: string;
  ts: number;
}

function promisify<T>(req: IDBRequest<T>): Promise<T> {
  return new Promise((resolve, reject) => {
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error);
  });
}

function txDone(tx: IDBTransaction): Promise<void> {
  return new Promise((resolve, reject) => {
    tx.oncomplete = () => resolve();
    tx.onerror = () => reject(tx.error);
    tx.onabort = () => reject(tx.error);
  });
}

/** Arka planda IndexedDB üzerinde çalışan çeviri önbelleği. */
export class IdbCache implements CacheStore {
  private dbPromise: Promise<IDBDatabase> | null = null;

  private db(): Promise<IDBDatabase> {
    if (!this.dbPromise) {
      this.dbPromise = new Promise((resolve, reject) => {
        const req = indexedDB.open(DB_NAME, 1);
        req.onupgradeneeded = () => {
          const store = req.result.createObjectStore(STORE);
          store.createIndex('ts', 'ts');
        };
        req.onsuccess = () => resolve(req.result);
        req.onerror = () => {
          this.dbPromise = null;
          reject(req.error);
        };
      });
    }
    return this.dbPromise;
  }

  async getMany(keys: string[]): Promise<(string | undefined)[]> {
    if (!keys.length) return [];
    const db = await this.db();
    const store = db.transaction(STORE, 'readonly').objectStore(STORE);
    const values = await Promise.all(keys.map((k) => promisify<Entry | undefined>(store.get(k))));
    return values.map((e) => e?.v);
  }

  async setMany(entries: [string, string][]): Promise<void> {
    if (!entries.length) return;
    const db = await this.db();
    const tx = db.transaction(STORE, 'readwrite');
    const store = tx.objectStore(STORE);
    const ts = Date.now();
    for (const [k, v] of entries) store.put({ v, ts } satisfies Entry, k);
    await txDone(tx);
  }

  async count(): Promise<number> {
    const db = await this.db();
    return promisify(db.transaction(STORE, 'readonly').objectStore(STORE).count());
  }

  async clear(): Promise<void> {
    const db = await this.db();
    const tx = db.transaction(STORE, 'readwrite');
    tx.objectStore(STORE).clear();
    await txDone(tx);
  }

  /** Süresi dolan kayıtları siler. */
  async prune(ttlDays: number): Promise<void> {
    const db = await this.db();
    const tx = db.transaction(STORE, 'readwrite');
    const index = tx.objectStore(STORE).index('ts');
    const range = IDBKeyRange.upperBound(Date.now() - ttlDays * 86_400_000);
    await new Promise<void>((resolve, reject) => {
      const req = index.openCursor(range);
      req.onsuccess = () => {
        const cursor = req.result;
        if (!cursor) return resolve();
        cursor.delete();
        cursor.continue();
      };
      req.onerror = () => reject(req.error);
    });
    await txDone(tx);
  }
}

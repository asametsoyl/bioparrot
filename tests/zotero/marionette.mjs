// Basit Marionette istemcisi (Zotero/Firefox uzaktan kontrol protokolü).
import net from 'node:net';

export class Marionette {
  constructor(port = 2829, host = '127.0.0.1') {
    this.port = port;
    this.host = host;
    this.id = 0;
    this.pending = new Map();
    this.buffer = Buffer.alloc(0);
  }

  async connect(timeoutMs = 60000) {
    const start = Date.now();
    for (;;) {
      try {
        await new Promise((resolve, reject) => {
          this.socket = net.connect(this.port, this.host, resolve);
          this.socket.once('error', reject);
        });
        break;
      } catch (err) {
        if (Date.now() - start > timeoutMs) throw err;
        await new Promise((r) => setTimeout(r, 1000));
      }
    }
    this.socket.on('data', (chunk) => this.onData(chunk));
    this.hello = await new Promise((resolve) => (this.onHello = resolve));
    return this.hello;
  }

  onData(chunk) {
    this.buffer = Buffer.concat([this.buffer, chunk]);
    for (;;) {
      const colon = this.buffer.indexOf(':');
      if (colon < 0) return;
      const len = Number(this.buffer.subarray(0, colon).toString());
      if (this.buffer.length < colon + 1 + len) return;
      const json = JSON.parse(this.buffer.subarray(colon + 1, colon + 1 + len).toString());
      this.buffer = this.buffer.subarray(colon + 1 + len);
      if (!Array.isArray(json)) {
        this.onHello?.(json);
        continue;
      }
      const [, id, error, result] = json;
      const p = this.pending.get(id);
      this.pending.delete(id);
      if (error) p?.reject(new Error(`${error.error}: ${error.message}\n${error.stacktrace ?? ''}`));
      else p?.resolve(result);
    }
  }

  send(name, params = {}) {
    const id = ++this.id;
    const msg = Buffer.from(JSON.stringify([0, id, name, params]));
    this.socket.write(`${msg.length}:`);
    this.socket.write(msg);
    return new Promise((resolve, reject) => this.pending.set(id, { resolve, reject }));
  }

  async exec(fn, args = [], timeoutMs = 120000) {
    const script = `return (${fn.toString()})(...arguments);`;
    const res = await this.send('WebDriver:ExecuteScript', { script, args, scriptTimeout: timeoutMs, newSandbox: true });
    return res?.value;
  }

  close() {
    this.socket?.destroy();
  }
}

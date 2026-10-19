import { t } from '@/core/i18n';
import { sendToBackground } from '@/lib/messaging';
import { ICONS, createUiHost } from './ui-host';

const CSS = `
.bp-toast { position: fixed; right: 18px; bottom: 18px; z-index: 2147483647; max-width: 380px; display: flex; gap: 10px; align-items: flex-start; padding: 12px 12px 12px 14px; border-radius: 12px; background: var(--bp-bg); color: var(--bp-fg); box-shadow: var(--bp-shadow); border: 1px solid var(--bp-line); border-left: 4px solid #dc2626; font-size: 13px; animation: bp-in .18s ease-out; }
@keyframes bp-in { from { transform: translateY(8px); opacity: 0; } }
.bp-toast .msg { flex: 1; }
.bp-toast .msg b { display: block; margin-bottom: 2px; }
.bp-toast .act { margin-top: 8px; display: inline-flex; height: 26px; align-items: center; padding: 0 10px; border-radius: 7px; background: var(--bp-green); color: #fff; font-size: 12px; font-weight: 600; }
.bp-toast .x { color: var(--bp-muted); width: 22px; height: 22px; display: grid; place-items: center; border-radius: 6px; }
.bp-toast .x:hover { background: var(--bp-soft); }
`;

let host: HTMLElement | null = null;
let lastMessage = '';

/** Sağ altta kısa bir hata bildirimi gösterir (aynı mesajı tekrar göstermez). */
export function showErrorToast(message: string, withSettings: boolean): void {
  if (message === lastMessage && host?.isConnected) return;
  lastMessage = message;
  host?.remove();
  const ui = createUiHost('toast', CSS);
  host = ui.host;
  const box = document.createElement('div');
  box.className = 'bp-root bp-toast';
  const msg = document.createElement('div');
  msg.className = 'msg';
  const title = document.createElement('b');
  title.textContent = `BioParrot · ${t('cnt.failed')}`;
  msg.appendChild(title);
  msg.appendChild(document.createTextNode(message));
  if (withSettings) {
    const btn = document.createElement('button');
    btn.className = 'act';
    btn.textContent = t('common.settings');
    btn.addEventListener('click', () => void sendToBackground({ type: 'openOptions', section: 'engines' }).catch(() => {}));
    msg.appendChild(document.createElement('br'));
    msg.appendChild(btn);
  }
  const close = document.createElement('button');
  close.className = 'x';
  close.innerHTML = ICONS.close;
  close.addEventListener('click', () => host?.remove());
  box.append(msg, close);
  ui.root.appendChild(box);
  setTimeout(() => {
    if (host === ui.host) host?.remove();
  }, 9000);
}

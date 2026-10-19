import { browser } from 'wxt/browser';
import { t } from '@/core/i18n';
import { isSameLanguage, normalizeLang } from '@/core/languages';
import type { Settings } from '@/core/types';
import { sendToBackground } from '@/lib/messaging';
import { createUiHost, eventTarget, isEditable } from './ui-host';

const CSS = `
.bp-input-badge { position: fixed; z-index: 2147483647; display: inline-flex; align-items: center; gap: 6px; padding: 4px 10px; border-radius: 999px; font-size: 12px; background: var(--bp-bg); color: var(--bp-fg); box-shadow: var(--bp-shadow); border: 1px solid var(--bp-line); }
.bp-input-badge.err { color: #dc2626; }
`;

/** Metin kutusunda art arda boşluk tuşuyla yazılanı çevirir. */
export class InputTranslator {
  private presses: number[] = [];
  private host: HTMLElement | null = null;
  private busy = false;

  constructor(private settings: () => Settings) {}

  attach(): void {
    document.addEventListener('keydown', this.onKeyDown, true);
  }

  detach(): void {
    document.removeEventListener('keydown', this.onKeyDown, true);
    this.host?.remove();
  }

  private onKeyDown = (e: KeyboardEvent) => {
    const s = this.settings();
    if (s.input.trigger === 'off') return;
    if (e.key !== ' ' || e.ctrlKey || e.metaKey || e.altKey || e.isComposing) {
      this.presses = [];
      return;
    }
    const el = eventTarget(e);
    if (!isEditable(el)) return;
    const now = Date.now();
    this.presses = this.presses.filter((ts) => now - ts < 450);
    this.presses.push(now);
    const needed = s.input.trigger === 'double-space' ? 2 : 3;
    if (this.presses.length < needed) return;
    this.presses = [];
    e.preventDefault();
    e.stopPropagation();
    void this.translate(el);
  };

  private async translate(el: HTMLElement): Promise<void> {
    if (this.busy) return;
    const s = this.settings();
    const text = readText(el).replace(/[\s ]+$/, '');
    if (!text.trim()) return;
    this.busy = true;
    const badge = this.badge(el, t('common.translating'));
    try {
      let to = s.input.targetLanguage;
      if (s.input.smartReverse) {
        try {
          const r = await browser.i18n.detectLanguage(text);
          const lang = normalizeLang(r?.languages?.[0]?.language);
          if (r?.isReliable && isSameLanguage(lang, to)) to = s.targetLanguage;
        } catch {
          /* algılama yok */
        }
      }
      const res = await sendToBackground({
        type: 'translate',
        req: {
          texts: [text],
          from: 'auto',
          to,
          engineId: s.input.engineId || s.engineId,
          expertId: s.input.expertId || undefined,
          rich: false,
        },
      });
      const out = res.translations[0];
      if (out) writeText(el, out);
      badge.remove();
    } catch (err) {
      badge.classList.add('err');
      badge.textContent = t('err.generic', { msg: String((err as Error)?.message ?? err) });
      setTimeout(() => badge.remove(), 4000);
    } finally {
      this.busy = false;
    }
  }

  private badge(el: HTMLElement, text: string): HTMLElement {
    if (!this.host) {
      const { host, root } = createUiHost('input', CSS);
      this.host = host;
      const wrap = document.createElement('div');
      wrap.className = 'bp-root';
      root.appendChild(wrap);
    }
    const wrap = this.host.shadowRoot!.querySelector('.bp-root')!;
    const badge = document.createElement('div');
    badge.className = 'bp-input-badge';
    badge.innerHTML = '<span class="bp-spin"></span>';
    badge.appendChild(document.createTextNode(text));
    const r = el.getBoundingClientRect();
    badge.style.left = `${Math.max(8, Math.min(r.right - 150, window.innerWidth - 220))}px`;
    badge.style.top = `${r.bottom + 6 > window.innerHeight - 40 ? Math.max(8, r.top - 34) : r.bottom + 6}px`;
    wrap.appendChild(badge);
    return badge;
  }
}

function readText(el: HTMLElement): string {
  if (el instanceof HTMLInputElement || el instanceof HTMLTextAreaElement) return el.value;
  return el.innerText;
}

/** Metni geri alınabilir şekilde değiştirir (React/Vue ve zengin editörlerle uyumlu). */
function writeText(el: HTMLElement, text: string): void {
  el.focus();
  if (el instanceof HTMLInputElement || el instanceof HTMLTextAreaElement) {
    el.select();
    if (!document.execCommand('insertText', false, text)) {
      const proto = el instanceof HTMLInputElement ? HTMLInputElement.prototype : HTMLTextAreaElement.prototype;
      Object.getOwnPropertyDescriptor(proto, 'value')?.set?.call(el, text);
      el.dispatchEvent(new Event('input', { bubbles: true }));
    }
    return;
  }
  const sel = window.getSelection();
  const range = document.createRange();
  range.selectNodeContents(el);
  sel?.removeAllRanges();
  sel?.addRange(range);
  if (!document.execCommand('insertText', false, text)) {
    el.textContent = text;
    el.dispatchEvent(new InputEvent('input', { bubbles: true, inputType: 'insertText', data: text }));
  }
}

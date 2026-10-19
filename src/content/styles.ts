import type { TranslationStyle } from '@/core/types';

/** Sayfaya (ve açık gölge köklerine) eklenen çeviri stilleri. */
export const PAGE_CSS = `
font.bp-wrap { display: inline; }
font.bp-wrap.bp-block { display: block !important; margin: 0.3em 0 0 0 !important; }
font.bp-wrap.bp-block.bp-flush { margin-top: 0 !important; }
.bp-text { font-size: calc(1em * var(--bp-font-scale, 1)); color: var(--bp-text-color, inherit); }
.bp-text a { color: inherit; }
.bp-wrap .bp-text code, .bp-wrap .bp-text kbd { font-size: 0.92em; }

.bp-style-underline .bp-text { text-decoration: underline dotted var(--bp-color, #22c55e) 2px !important; text-underline-offset: 0.3em; text-decoration-skip-ink: none; }
.bp-style-dashed .bp-text { text-decoration: underline dashed var(--bp-color, #22c55e) 1.5px !important; text-underline-offset: 0.3em; text-decoration-skip-ink: none; }
.bp-style-wavy .bp-text { text-decoration: underline wavy var(--bp-color, #22c55e) 1px !important; text-underline-offset: 0.3em; text-decoration-skip-ink: none; }
.bp-style-highlight .bp-text { background-color: color-mix(in srgb, var(--bp-color, #22c55e) 20%, transparent); border-radius: 3px; padding: 0 0.15em; -webkit-box-decoration-break: clone; box-decoration-break: clone; }
.bp-style-marker .bp-text { background-image: linear-gradient(transparent 62%, color-mix(in srgb, var(--bp-color, #22c55e) 40%, transparent) 62%); -webkit-box-decoration-break: clone; box-decoration-break: clone; }
.bp-style-blur .bp-text { filter: blur(4px); transition: filter 0.2s ease; }
.bp-style-blur:hover .bp-text { filter: none; }
.bp-style-faded .bp-text { opacity: 0.62; }
.bp-style-italic .bp-text { font-style: italic; }
.bp-style-bold .bp-text { font-weight: 600; }
.bp-style-border.bp-block .bp-text { display: block; border: 1px dashed var(--bp-color, #22c55e); border-radius: 6px; padding: 0.3em 0.55em; }
.bp-style-border.bp-inline .bp-text { border: 1px dashed var(--bp-color, #22c55e); border-radius: 4px; padding: 0 0.3em; }
.bp-style-quote.bp-block .bp-text { display: block; border-left: 3px solid var(--bp-color, #22c55e); padding-left: 0.6em; }
.bp-style-card.bp-block .bp-text { display: block; background: color-mix(in srgb, var(--bp-color, #22c55e) 8%, transparent); border: 1px solid color-mix(in srgb, var(--bp-color, #22c55e) 28%, transparent); border-radius: 8px; padding: 0.45em 0.7em; }
.bp-style-card.bp-inline .bp-text { background: color-mix(in srgb, var(--bp-color, #22c55e) 10%, transparent); border-radius: 4px; padding: 0 0.3em; }

font.bp-loading { display: inline-block !important; width: 0.75em; height: 0.75em; margin: 0 0.35em; vertical-align: -0.08em; border: 2px solid var(--bp-color, #22c55e); border-right-color: transparent; border-radius: 50%; animation: bp-spin 0.7s linear infinite; opacity: 0.85; box-sizing: border-box; }
@keyframes bp-spin { to { transform: rotate(360deg); } }
font.bp-error { display: inline-block; margin: 0 0.35em; color: #dc2626; cursor: pointer; font-size: 0.85em; font-weight: 600; }
font.bp-error:hover { text-decoration: underline; }

.bp-orig-hidden { display: none !important; }
font.bp-wrap.bp-flush { margin-top: 0 !important; }
font.bp-wrap.bp-flush > .bp-sp { display: none; }

.bp-hover-target { outline: 2px dashed color-mix(in srgb, var(--bp-color, #22c55e) 70%, transparent) !important; outline-offset: 3px; border-radius: 2px; }
`;

const STYLE_ATTR = 'data-bp-style';

/** Belgeye veya gölge köke stil ekler (bir kez). */
export function injectPageStyle(root: Document | ShadowRoot): void {
  const host = root instanceof Document ? root.head ?? root.documentElement : root;
  if (!host || host.querySelector?.(`style[${STYLE_ATTR}]`)) return;
  const style = (root instanceof Document ? root : root.ownerDocument).createElement('style');
  style.setAttribute(STYLE_ATTR, '');
  style.textContent = PAGE_CSS;
  host.appendChild(style);
}

/** Renk ve boyut değişkenlerini kök öğeye uygular (gölge köklere miras geçer). */
export function applyStyleVariables(style: TranslationStyle): void {
  const el = document.documentElement;
  el.style.setProperty('--bp-color', style.color || '#22c55e');
  el.style.setProperty('--bp-font-scale', String(style.fontScale || 1));
  if (style.textColor) el.style.setProperty('--bp-text-color', style.textColor);
  else el.style.removeProperty('--bp-text-color');
}

export function clearStyleVariables(): void {
  const el = document.documentElement;
  el.style.removeProperty('--bp-color');
  el.style.removeProperty('--bp-font-scale');
  el.style.removeProperty('--bp-text-color');
}

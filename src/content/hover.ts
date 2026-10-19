import { CLS } from '@/core/constants';
import type { Settings } from '@/core/types';
import type { PageTranslator } from './page-translator';
import { deepElementFromPoint, isEditable } from './ui-host';

/** Tetikleyici tuşa basılıyken farenin altındaki paragrafı çevirir. */
export class HoverTranslator {
  private x = -1;
  private y = -1;
  private highlighted: Element | null = null;
  /** Kısayolla geçici olarak kapatılabilir. */
  enabled = true;

  constructor(
    private translator: PageTranslator,
    private settings: () => Settings,
  ) {}

  attach(): void {
    document.addEventListener('mousemove', this.onMove, { passive: true, capture: true });
    window.addEventListener('keydown', this.onKeyDown, true);
    window.addEventListener('keyup', this.onKeyUp, true);
    window.addEventListener('blur', this.clear);
  }

  detach(): void {
    document.removeEventListener('mousemove', this.onMove, true);
    window.removeEventListener('keydown', this.onKeyDown, true);
    window.removeEventListener('keyup', this.onKeyUp, true);
    window.removeEventListener('blur', this.clear);
    this.clear();
  }

  private onMove = (e: MouseEvent) => {
    this.x = e.clientX;
    this.y = e.clientY;
  };

  private onKeyDown = (e: KeyboardEvent) => {
    const key = this.settings().hover.key;
    if (!this.enabled || key === 'off' || e.key !== key || e.repeat) return;
    if (isEditable(document.activeElement) || this.x < 0) return;
    // Başka tuşlarla birlikte basıldıysa (kısayol) karışma.
    const others = (['Shift', 'Alt', 'Control', 'Meta'] as const).filter((k) => k !== key && e.getModifierState(k));
    if (others.length) return;
    const el = deepElementFromPoint(this.x, this.y);
    if (!el || el.closest('[data-bp-ui]')) return;
    const target = this.translator.hoverTarget(el);
    this.highlight(target);
    this.translator.toggleAt(el);
  };

  private onKeyUp = (e: KeyboardEvent) => {
    if (e.key === this.settings().hover.key) this.clear();
  };

  private highlight(el: Element | null): void {
    this.clear();
    if (!el || el === document.body || el === document.documentElement) return;
    el.classList.add(CLS.hoverTarget);
    this.highlighted = el;
  }

  private clear = () => {
    this.highlighted?.classList.remove(CLS.hoverTarget);
    this.highlighted = null;
  };
}

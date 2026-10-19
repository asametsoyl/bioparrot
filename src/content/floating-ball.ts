import { markSvg } from '@/core/brand';
import { t } from '@/core/i18n';
import { ICONS, createUiHost } from './ui-host';

export interface BallState {
  active: boolean;
  pending: number;
  isPdf: boolean;
}

export interface BallActions {
  toggle(): void;
  toggleMode(): void;
  settings(): void;
  hideSite(): void;
  openPdf(): void;
  savePosition(top: number): void;
}

const CSS = `
.bp-ball-wrap { position: fixed; z-index: 2147483646; display: flex; flex-direction: column; align-items: center; gap: 6px; user-select: none; touch-action: none; }
.bp-ball-wrap.right { right: 0; }
.bp-ball-wrap.left { left: 0; }
.bp-ball { position: relative; width: 40px; height: 40px; border-radius: 22px 0 0 22px; background: var(--bp-bg); box-shadow: var(--bp-shadow); display: grid; place-items: center; transition: width .18s ease, transform .18s ease; cursor: pointer; border: 1px solid var(--bp-line); border-right: 0; }
.left .bp-ball { border-radius: 0 22px 22px 0; border-left: 0; border-right: 1px solid var(--bp-line); }
.bp-ball-wrap:hover .bp-ball { width: 46px; }
.bp-ball svg { display: block; }
.bp-ball.active { background: linear-gradient(135deg, var(--bp-green-2), var(--bp-teal)); border-color: transparent; }
.bp-ball .badge { position: absolute; top: 2px; right: 4px; width: 14px; height: 14px; border-radius: 50%; background: #fff; color: var(--bp-green); display: none; place-items: center; box-shadow: 0 1px 3px rgba(0,0,0,.25); }
.left .bp-ball .badge { right: auto; left: 4px; }
.bp-ball.active .badge { display: grid; }
.bp-ball .ring { position: absolute; inset: 4px; border-radius: 50%; border: 2px solid rgba(255,255,255,.0); border-top-color: var(--bp-beak); animation: bp-rot .8s linear infinite; display: none; }
.bp-ball.busy .ring { display: block; }
.bp-tools { display: flex; flex-direction: column; gap: 4px; padding: 4px; border-radius: 12px; background: var(--bp-bg); box-shadow: var(--bp-shadow); border: 1px solid var(--bp-line); opacity: 0; transform: translateY(-4px) scale(.96); pointer-events: none; transition: opacity .15s ease, transform .15s ease; }
.bp-ball-wrap:hover .bp-tools, .bp-ball-wrap.open .bp-tools { opacity: 1; transform: none; pointer-events: auto; }
.right .bp-tools { margin-right: 6px; }
.left .bp-tools { margin-left: 6px; }
.bp-tool { width: 30px; height: 30px; border-radius: 8px; display: grid; place-items: center; color: var(--bp-muted); position: relative; }
.bp-tool:hover { background: var(--bp-soft); color: var(--bp-green); }
.bp-tool[data-tip]:hover::after { content: attr(data-tip); position: absolute; top: 50%; transform: translateY(-50%); white-space: nowrap; background: #0f172a; color: #fff; font-size: 12px; padding: 4px 8px; border-radius: 6px; pointer-events: none; }
.right .bp-tool[data-tip]:hover::after { right: 38px; }
.left .bp-tool[data-tip]:hover::after { left: 38px; }
`;

/** Sayfanın kenarındaki yüzen BioParrot düğmesi. */
export class FloatingBall {
  private host: HTMLElement | null = null;
  private wrap!: HTMLElement;
  private ball!: HTMLElement;
  private mark!: HTMLElement;
  private dragging = false;

  constructor(
    private actions: BallActions,
    private side: 'left' | 'right',
    private top: number,
  ) {}

  mount(): void {
    if (this.host) return;
    const { host, root } = createUiHost('floating-ball', CSS);
    this.host = host;
    const container = document.createElement('div');
    root.appendChild(container);
    container.innerHTML = `
      <div class="bp-root bp-ball-wrap ${this.side}">
        <div class="bp-ball" role="button" tabindex="0">
          <span class="mark"></span><span class="ring"></span><span class="badge">${ICONS.check}</span>
        </div>
        <div class="bp-tools">
          <button class="bp-tool" data-act="mode">${ICONS.mode}</button>
          <button class="bp-tool" data-act="pdf" hidden>${ICONS.pdf}</button>
          <button class="bp-tool" data-act="settings">${ICONS.settings}</button>
          <button class="bp-tool" data-act="hide">${ICONS.eyeOff}</button>
        </div>
      </div>`;
    this.wrap = root.querySelector('.bp-ball-wrap')!;
    this.ball = root.querySelector('.bp-ball')!;
    this.mark = root.querySelector('.mark')!;
    this.wrap.style.top = `${this.top}%`;
    this.applyLabels();

    root.querySelectorAll<HTMLButtonElement>('.bp-tool').forEach((btn) =>
      btn.addEventListener('click', (e) => {
        e.stopPropagation();
        const act = btn.dataset.act;
        if (act === 'mode') this.actions.toggleMode();
        else if (act === 'settings') this.actions.settings();
        else if (act === 'pdf') this.actions.openPdf();
        else if (act === 'hide') this.actions.hideSite();
      }),
    );
    this.ball.addEventListener('keydown', (e) => {
      if (e.key === 'Enter' || e.key === ' ') {
        e.preventDefault();
        this.actions.toggle();
      }
    });
    this.setupDrag();
    this.update({ active: false, pending: 0, isPdf: false });
  }

  unmount(): void {
    this.host?.remove();
    this.host = null;
  }

  applyLabels(): void {
    if (!this.host) return;
    const root = this.host.shadowRoot!;
    const set = (act: string, tip: string) => {
      const b = root.querySelector<HTMLElement>(`[data-act="${act}"]`);
      if (b) {
        b.dataset.tip = tip;
        b.setAttribute('aria-label', tip);
      }
    };
    set('mode', t('fb.toggleMode'));
    set('pdf', t('fb.pdf'));
    set('settings', t('fb.settings'));
    set('hide', t('fb.hideSite'));
  }

  update(state: BallState): void {
    if (!this.host) return;
    this.ball.classList.toggle('active', state.active);
    this.ball.classList.toggle('busy', state.pending > 0);
    this.mark.innerHTML = markSvg(24, state.active ? '#ffffff' : '#16a34a');
    this.ball.title = state.isPdf ? t('fb.pdf') : state.active ? t('fb.original') : `${t('fb.translate')} (Alt+A)`;
    const pdf = this.host.shadowRoot!.querySelector<HTMLElement>('[data-act="pdf"]');
    if (pdf) pdf.hidden = !state.isPdf;
  }

  setPosition(side: 'left' | 'right', top: number): void {
    this.side = side;
    this.top = top;
    if (!this.host) return;
    this.wrap.classList.remove('left', 'right');
    this.wrap.classList.add(side);
    this.wrap.style.top = `${top}%`;
  }

  private setupDrag(): void {
    let startY = 0;
    let startTop = 0;
    let moved = false;
    this.ball.addEventListener('pointerdown', (e) => {
      if (e.button !== 0) return;
      startY = e.clientY;
      startTop = this.wrap.getBoundingClientRect().top;
      moved = false;
      this.dragging = true;
      this.ball.setPointerCapture(e.pointerId);
    });
    this.ball.addEventListener('pointermove', (e) => {
      if (!this.dragging) return;
      const dy = e.clientY - startY;
      if (Math.abs(dy) > 4) moved = true;
      if (!moved) return;
      const top = Math.min(Math.max(startTop + dy, 8), window.innerHeight - 60);
      this.wrap.style.top = `${top}px`;
    });
    this.ball.addEventListener('pointerup', (e) => {
      if (!this.dragging) return;
      this.dragging = false;
      this.ball.releasePointerCapture(e.pointerId);
      if (moved) {
        const pct = Math.round((this.wrap.getBoundingClientRect().top / window.innerHeight) * 1000) / 10;
        this.top = pct;
        this.wrap.style.top = `${pct}%`;
        this.actions.savePosition(pct);
      } else {
        this.actions.toggle();
      }
    });
  }
}

import { isLlmEngine } from '../engines/meta';
import { escapeHtml, markupToPlain } from '../markup';
import { lookupSectionHeading, matchInitialCase } from '../section-headings';
import { getAllExperts } from '../experts';
import type { Settings, TranslateRequest, TranslateResponse, TranslationContext } from '../types';
import { analyzePage, estimateBodySize, finalizeDocument, guessTitleAndAbstract, type PageLayout } from './layout';
import { readRawPage, type PDFDocumentProxy } from './pdfjs';

export interface PdfJobHost {
  translate(req: TranslateRequest): Promise<TranslateResponse>;
  extractTerms(text: string, to: string, engineId?: string, expertId?: string): Promise<{ source: string; target: string }[]>;
}

export type PageStatus = 'pending' | 'translating' | 'done' | 'error';

export interface JobState {
  phase: 'idle' | 'analyzing' | 'translating' | 'done' | 'error' | 'stopped';
  analyzed: number;
  total: number;
  pagesDone: number;
  pageStatus: PageStatus[];
  error?: string;
  rich: boolean;
}

/** Bir PDF'in sayfa düzenini çıkarıp tüm çevrilebilir blokları çeviren iş. */
export class PdfTranslationJob {
  layouts: PageLayout[] = [];
  translations = new Map<string, string>();
  state: JobState;
  context: TranslationContext = {};
  private stopped = false;
  private expertId: string;

  constructor(
    private doc: PDFDocumentProxy,
    private host: PdfJobHost,
    private settings: Settings,
    private onUpdate: (job: PdfTranslationJob) => void,
  ) {
    this.state = {
      phase: 'idle',
      analyzed: 0,
      total: doc.numPages,
      pagesDone: 0,
      pageStatus: new Array(doc.numPages).fill('pending'),
      rich: true,
    };
    const s = settings;
    // PDF çevirmeni akademik metinler içindir: varsayılan "Genel" ise akademik uzman kullanılır.
    this.expertId = s.academic.autoExpert && s.expertId === 'general' ? s.academic.expertId : s.expertId;
  }

  get expert(): string {
    return this.expertId;
  }

  stop(): void {
    this.stopped = true;
    if (this.state.phase === 'translating' || this.state.phase === 'analyzing') this.update({ phase: 'stopped' });
  }

  private update(patch: Partial<JobState>): void {
    this.state = { ...this.state, ...patch };
    this.onUpdate(this);
  }

  /** Tüm sayfaların düzenini çıkarır (hızlıdır, çeviriden önce yapılır). */
  async analyze(): Promise<void> {
    if (this.layouts.length === this.doc.numPages) return;
    this.update({ phase: 'analyzing', analyzed: 0 });
    const raws = [];
    for (let i = 0; i < this.doc.numPages; i++) {
      if (this.stopped) return;
      raws.push(await readRawPage(this.doc, i));
      this.update({ analyzed: i + 1 });
    }
    // Gövde yazı boyutu tüm belgeden hesaplanır (şekil içi metinleri ayırmak için).
    const prelim = raws.map((r) => analyzePage(r));
    const body = estimateBodySize(prelim.flatMap((p) => p.blocks.flatMap((b) => b.lines)));
    this.layouts = raws.map((r) => analyzePage(r, body));
    finalizeDocument(this.layouts, { skipReferences: this.settings.academic.skipReferences });
    const { title, abstract } = guessTitleAndAbstract(this.layouts);
    this.context = { title, summary: abstract };
    this.onUpdate(this);
  }

  async run(): Promise<void> {
    this.stopped = false;
    try {
      await this.analyze();
      if (this.stopped) return;
      const s = this.settings;
      const engine = s.engines[s.engineId];
      const academic = getAllExperts(s).find((e) => e.id === this.expertId)?.group === 'academic';
      const context: TranslationContext = s.page.contextAware ? { ...this.context } : {};
      this.update({ phase: 'translating' });
      if (isLlmEngine(engine) && academic && s.academic.termConsistency && (this.context.title || this.context.summary)) {
        try {
          const terms = await this.host.extractTerms(`${this.context.title}\n\n${this.context.summary}`, s.targetLanguage, s.engineId, this.expertId);
          if (terms.length) context.terms = terms;
          this.context.terms = terms;
        } catch {
          /* terim listesi isteğe bağlı */
        }
      }

      const queue = this.layouts.map((_, i) => i).filter((i) => this.state.pageStatus[i] !== 'done');
      const workers = Array.from({ length: Math.min(3, queue.length) }, async () => {
        while (queue.length && !this.stopped) {
          const index = queue.shift()!;
          await this.translatePage(index, context);
        }
      });
      await Promise.all(workers);
      if (!this.stopped) {
        const failed = this.state.pageStatus.filter((p) => p === 'error').length;
        this.update({ phase: failed && failed === this.layouts.length ? 'error' : 'done' });
      }
    } catch (err) {
      this.update({ phase: 'error', error: String((err as Error)?.message ?? err) });
    }
  }

  private async translatePage(index: number, context: TranslationContext): Promise<void> {
    const layout = this.layouts[index];
    const s0 = this.settings;
    let blocks = layout.blocks.filter((b) => b.translatable && !this.translations.has(b.id) && b.markup.trim());
    // Standart bölüm başlıkları klasik motorlara gönderilmeden sözlükten çevrilir ("Abstract" → "Özet").
    if (!isLlmEngine(s0.engines[s0.engineId])) {
      blocks = blocks.filter((b) => {
        const hit = b.atoms.length === 0 ? lookupSectionHeading(b.plain, s0.targetLanguage) : undefined;
        if (hit) this.translations.set(b.id, escapeHtml(hit));
        return !hit;
      });
    }
    const status = [...this.state.pageStatus];
    if (!blocks.length) {
      status[index] = 'done';
      this.update({ pageStatus: status, pagesDone: this.state.pagesDone + 1 });
      return;
    }
    status[index] = 'translating';
    this.update({ pageStatus: status });
    const s = this.settings;
    try {
      const res = await this.host.translate({
        texts: blocks.map((b) => b.markup),
        from: s.sourceLanguage,
        to: s.targetLanguage,
        engineId: s.engineId,
        expertId: this.expertId,
        rich: true,
        context,
      });
      const norm = (x: string) => x.replace(/\s+/g, ' ').trim().toLowerCase();
      blocks.forEach((b, i) => {
        const t = res.translations[i];
        if (!t) return;
        const plain = res.rich ? markupToPlain(t) : t;
        // Aynı kalan metinler (isimler, kısaltmalar) için özgün sayfa görünür bırakılır.
        if (norm(plain) === norm(b.plain)) return;
        this.translations.set(b.id, b.kind === 'heading' && !res.rich ? matchInitialCase(b.plain, t, s.targetLanguage) : t);
      });
      const next = [...this.state.pageStatus];
      next[index] = 'done';
      this.update({ pageStatus: next, pagesDone: this.state.pagesDone + 1, rich: res.rich });
    } catch (err) {
      const next = [...this.state.pageStatus];
      next[index] = 'error';
      this.update({ pageStatus: next, error: String((err as Error)?.message ?? err) });
    }
  }
}

import { t } from '@/core/i18n';
import { getService, translateText } from './service';
import { getPref, loadSettings } from './prefs';

declare const Zotero: any;

/** Son çevrilen seçimler (vurguya yorum eklemek için). */
const recent = new Map<string, string>();

export function normalizeSelection(text: string): string {
  return text
    .replace(/([a-zçğıöşü])-\s*\n\s*([a-zçğıöşü])/gi, '$1$2')
    .replace(/\s*\n\s*/g, ' ')
    .replace(/\s{2,}/g, ' ')
    .trim();
}

function remember(source: string, translation: string): void {
  recent.set(source.toLowerCase(), translation);
  if (recent.size > 30) recent.delete(recent.keys().next().value!);
}

export function translationFor(annotationText: string): string | undefined {
  return recent.get(normalizeSelection(annotationText).toLowerCase());
}

function contextFor(reader: any): { title?: string; summary?: string } {
  try {
    const item = Zotero.Items.get(reader.itemID);
    const parent = item?.parentItem ?? item;
    return { title: parent?.getField?.('title') || undefined, summary: parent?.getField?.('abstractNote') || undefined };
  } catch {
    return {};
  }
}

function button(doc: Document, label: string, onClick: () => void): HTMLButtonElement {
  const b = doc.createElement('button');
  b.textContent = label;
  b.style.cssText =
    'font: inherit; font-size: 12px; padding: 3px 9px; border-radius: 6px; border: 1px solid rgba(127,127,127,.35); background: transparent; color: inherit; cursor: pointer;';
  b.addEventListener('click', (e) => {
    e.preventDefault();
    e.stopPropagation();
    onClick();
  });
  return b;
}

/** PDF okuyucunun seçim penceresine BioParrot çevirisini ekler. */
export function onTextSelectionPopup(event: { reader: any; doc: Document; params: any; append: (el: Element) => void }): void {
  const { reader, doc, params, append } = event;
  const raw: string = params?.annotation?.text ?? '';
  const text = normalizeSelection(raw);
  if (!text) return;

  const box = doc.createElement('div');
  box.style.cssText =
    'margin-top: 6px; padding: 8px 10px; max-width: 380px; border-top: 1px solid rgba(127,127,127,.25); font-size: 13px; line-height: 1.5; user-select: text;';
  const head = doc.createElement('div');
  head.style.cssText = 'display:flex; align-items:center; gap:6px; font-weight:600; font-size:12px; opacity:.8; margin-bottom:4px;';
  const icon = doc.createElement('span');
  icon.textContent = '🦜';
  const title = doc.createElement('span');
  const s = loadSettings();
  title.textContent = `${t('zot.popupTitle')} → ${s.targetLanguage}`;
  head.append(icon, title);
  const body = doc.createElement('div');
  body.style.cssText = 'white-space: pre-wrap; max-height: 220px; overflow: auto;';
  const actions = doc.createElement('div');
  actions.style.cssText = 'display:flex; gap:6px; margin-top:6px; flex-wrap: wrap;';
  box.append(head, body, actions);
  append(box);

  let translation = '';
  const run = async () => {
    body.textContent = t('common.translating');
    try {
      translation = await translateText(text, { context: contextFor(reader) });
      body.textContent = translation;
      remember(text, translation);
    } catch (err) {
      body.textContent = t('err.generic', { msg: String((err as Error)?.message ?? err) });
    }
  };

  actions.append(
    button(doc, t('common.copy'), () => {
      if (translation) Zotero.Utilities.Internal.copyTextToClipboard(translation);
    }),
    button(doc, t('sel.explain'), async () => {
      const explain = doc.createElement('div');
      explain.style.cssText = 'margin-top:8px; padding-top:6px; border-top: 1px dashed rgba(127,127,127,.3); white-space: pre-wrap;';
      explain.textContent = '…';
      box.appendChild(explain);
      try {
        const st = loadSettings();
        explain.textContent = await getService().explain(text, st.targetLanguage, contextFor(reader), st.engineId, st.expertId);
      } catch (err) {
        explain.textContent = String((err as Error)?.message ?? err);
      }
    }),
  );

  if (getPref<boolean>('autoTranslateSelection', true)) void run();
  else {
    body.textContent = '';
    actions.prepend(button(doc, t('zot.translateBtn'), () => void run()));
  }
}

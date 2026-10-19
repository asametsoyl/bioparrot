import { ENGINE_META } from '@/core/engines/meta';
import { expertName, getAllExperts } from '@/core/experts';
import { getLocale, t } from '@/core/i18n';
import { LANGUAGES, languageName } from '@/core/languages';
import { getPref, loadSettings, setPref } from './prefs';
import { getService } from './service';

const HTML = 'http://www.w3.org/1999/xhtml';

/** Zotero tercih penceresindeki BioParrot panelini oluşturur. */
export function initPrefsPane(doc: Document): void {
  const root = doc.getElementById('bioparrot-prefs-root');
  if (!root || root.childElementCount) return;
  const ui = getLocale();
  const el = <K extends keyof HTMLElementTagNameMap>(tag: K, props: { style?: string; [key: string]: unknown } = {}, ...children: (Node | string)[]) => {
    const node = doc.createElementNS(HTML, tag) as HTMLElementTagNameMap[K];
    const { style, ...rest } = props;
    Object.assign(node, rest);
    if (style) node.setAttribute('style', style);
    for (const c of children) node.append(c);
    return node;
  };
  const row = (label: string, control: HTMLElement) =>
    el('div', { style: 'display:flex; align-items:center; gap:12px; margin:8px 0;' }, el('label', { style: 'width:190px; flex:none;' }, label), control);
  const settings = loadSettings();

  const select = (pref: string, options: [string, string][], value: string) => {
    const s = el('select', { style: 'min-width:260px;' });
    for (const [v, label] of options) s.append(el('option', { value: v, textContent: label, selected: v === value }));
    s.addEventListener('change', () => setPref(pref, s.value));
    return s;
  };
  const input = (pref: string, type = 'text', placeholder = '') => {
    const i = el('input', { type, value: getPref<string>(pref, ''), placeholder, style: 'min-width:320px;' });
    i.addEventListener('change', () => setPref(pref, i.value.trim()));
    return i;
  };
  const checkbox = (pref: string, label: string) => {
    const c = el('input', { type: 'checkbox', checked: !!getPref<boolean>(pref, false) });
    c.addEventListener('change', () => setPref(pref, c.checked));
    return el('label', { style: 'display:flex; align-items:center; gap:8px; margin:6px 0;' }, c, label);
  };

  const languages: [string, string][] = [...LANGUAGES].sort((a, b) => languageName(a.code, ui).localeCompare(languageName(b.code, ui), ui)).map((l) => [l.code, `${languageName(l.code, ui)} — ${l.native}`]);
  const engines: [string, string][] = Object.values(settings.engines).map((e) => [e.id, `${e.name}${ENGINE_META[e.type].free ? ` (${t('common.free')})` : ''}`]);
  const experts: [string, string][] = getAllExperts(settings).map((e) => [e.id, `${e.icon ?? ''} ${expertName(e, ui)}`]);

  const status = el('span', { style: 'margin-left:10px;' });
  const testBtn = el('button', { textContent: t('zot.prefs.test') });
  testBtn.addEventListener('click', async () => {
    status.textContent = t('common.testing');
    try {
      const s = loadSettings();
      const res = await getService().translate({ texts: [t('eng.testSample')], from: 'auto', to: s.targetLanguage, engineId: s.engineId, noCache: true });
      status.textContent = `✓ ${res.translations[0]}`;
    } catch (err) {
      status.textContent = `✕ ${String((err as Error)?.message ?? err)}`;
    }
  });

  const importArea = el('textarea', { rows: 4, style: 'width:100%; font-family: monospace; font-size: 11px;', placeholder: '{"app":"bioparrot","settings":{…}}' });
  const importStatus = el('span', { style: 'margin-left:10px;' });
  const importBtn = el('button', { textContent: t('zot.prefs.importBtn') });
  importBtn.addEventListener('click', () => {
    try {
      const parsed = JSON.parse(importArea.value);
      const s = parsed?.settings ?? parsed;
      setPref('settingsJson', JSON.stringify(s));
      // Tekil tercihler içe aktarılanı ezmesin diye temizlenir.
      for (const k of ['openaiKey', 'anthropicKey', 'geminiKey', 'deeplKey', 'deepseekKey', 'openrouterKey', 'openaiModel', 'anthropicModel', 'geminiModel']) setPref(k, '');
      if (s?.targetLanguage) setPref('targetLanguage', s.targetLanguage);
      if (s?.engineId) setPref('engineId', s.engineId);
      if (s?.academic?.expertId) setPref('expertId', s.academic.expertId);
      importStatus.textContent = `✓ ${t('zot.prefs.imported')}`;
    } catch (err) {
      importStatus.textContent = `✕ ${String((err as Error)?.message ?? err)}`;
    }
  });

  root.append(
    el('h2', { textContent: 'BioParrot 🦜', style: 'margin: 0 0 4px;' }),
    el('p', { textContent: t('zot.prefs.intro'), style: 'margin: 0 0 12px; opacity: .8;' }),
    row(t('zot.prefs.target'), select('targetLanguage', languages, settings.targetLanguage)),
    row(t('zot.prefs.engine'), select('engineId', engines, settings.engineId)),
    row(t('zot.prefs.expert'), select('expertId', experts, getPref<string>('expertId', 'academic'))),
    el('div', { style: 'margin: 4px 0 12px;' }, testBtn, status),
    checkbox('autoTranslateSelection', t('zot.prefs.auto')),
    checkbox('annotationComment', t('zot.prefs.annotation')),
    el('h3', { textContent: t('zot.prefs.keys'), style: 'margin: 16px 0 4px;' }),
    el('p', { textContent: t('zot.prefs.keysHint'), style: 'margin: 0 0 6px; opacity: .75; font-size: 12px;' }),
    row('Claude (Anthropic)', input('anthropicKey', 'password', 'sk-ant-…')),
    row(`Claude ${t('zot.prefs.model')}`, input('anthropicModel', 'text', 'claude-opus-5-5 / claude-haiku-5-5')),
    row('OpenAI', input('openaiKey', 'password', 'sk-…')),
    row(`OpenAI ${t('zot.prefs.model')}`, input('openaiModel', 'text', 'gpt-4.1-mini')),
    row('Gemini', input('geminiKey', 'password')),
    row(`Gemini ${t('zot.prefs.model')}`, input('geminiModel', 'text', 'gemini-flash-latest')),
    row('DeepL', input('deeplKey', 'password', '…:fx')),
    row('DeepSeek', input('deepseekKey', 'password')),
    row('OpenRouter', input('openrouterKey', 'password')),
    el('h3', { textContent: t('zot.prefs.import'), style: 'margin: 16px 0 4px;' }),
    el('p', { textContent: t('zot.prefs.importHint'), style: 'margin: 0 0 6px; opacity: .75; font-size: 12px;' }),
    importArea,
    el('div', { style: 'margin-top:6px;' }, importBtn, importStatus),
  );
}

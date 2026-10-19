import { useEffect, useMemo, useState, type ReactNode } from 'react';
import { LOGO_SVG, markSvg } from '@/core/brand';
import { ENGINE_META, isEngineReady } from '@/core/engines/meta';
import { expertName, getAllExperts } from '@/core/experts';
import { getLocale, t } from '@/core/i18n';
import { LANGUAGES, languageName } from '@/core/languages';
import type { Settings } from '@/core/types';
import { useDebouncedCallback } from './use-debounced';

const LOGO_URL = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(LOGO_SVG)}`;

export function Logo({ size = 26, withText = true, sub }: { size?: number; withText?: boolean; sub?: string }) {
  return (
    <div className="logo">
      <img src={LOGO_URL} width={size} height={size} alt="" />
      {withText && (
        <div className="logo-text">
          <strong>BioParrot</strong>
          {sub && <span>{sub}</span>}
        </div>
      )}
    </div>
  );
}

export function Mark({ size = 18, color }: { size?: number; color?: string }) {
  return <span className="mark" dangerouslySetInnerHTML={{ __html: markSvg(size, color) }} />;
}

export function Switch({ checked, onChange, label }: { checked: boolean; onChange: (v: boolean) => void; label?: string }) {
  return (
    <label className="switch" aria-label={label}>
      <input type="checkbox" checked={checked} onChange={(e) => onChange(e.target.checked)} />
      <span className="track" />
    </label>
  );
}

export function Segmented<T extends string>({
  value,
  options,
  onChange,
}: {
  value: T;
  options: { value: T; label: string }[];
  onChange: (v: T) => void;
}) {
  return (
    <div className="seg" role="radiogroup">
      {options.map((o) => (
        <button key={o.value} type="button" role="radio" aria-checked={value === o.value} className={value === o.value ? 'on' : ''} onClick={() => onChange(o.value)}>
          {o.label}
        </button>
      ))}
    </div>
  );
}

export function Select<T extends string>({
  value,
  onChange,
  children,
  className = '',
  ariaLabel,
}: {
  value: T;
  onChange: (v: T) => void;
  children: ReactNode;
  className?: string;
  ariaLabel?: string;
}) {
  return (
    <select className={`select ${className}`} value={value} aria-label={ariaLabel} onChange={(e) => onChange(e.target.value as T)}>
      {children}
    </select>
  );
}

const POPULAR = ['en', 'tr', 'de', 'fr', 'es', 'it', 'pt', 'ru', 'zh-CN', 'zh-TW', 'ja', 'ko', 'ar'];

export function LanguageOptions({ includeAuto = false }: { includeAuto?: boolean }) {
  const ui = getLocale();
  const sorted = useMemo(
    () => [...LANGUAGES].sort((a, b) => languageName(a.code, ui).localeCompare(languageName(b.code, ui), ui)),
    [ui],
  );
  const label = (code: string) => {
    const lang = LANGUAGES.find((l) => l.code === code)!;
    const name = languageName(code, ui);
    return lang.native && lang.native !== name ? `${name} — ${lang.native}` : name;
  };
  return (
    <>
      {includeAuto && <option value="auto">{t('lang.auto')}</option>}
      <optgroup label="★">
        {POPULAR.map((c) => (
          <option key={`p-${c}`} value={c}>
            {label(c)}
          </option>
        ))}
      </optgroup>
      <optgroup label="A–Z">
        {sorted
          .filter((l) => !POPULAR.includes(l.code))
          .map((l) => (
            <option key={l.code} value={l.code}>
              {label(l.code)}
            </option>
          ))}
      </optgroup>
    </>
  );
}

export function LanguageSelect({
  value,
  onChange,
  includeAuto,
  ariaLabel,
}: {
  value: string;
  onChange: (v: string) => void;
  includeAuto?: boolean;
  ariaLabel?: string;
}) {
  return (
    <Select value={value} onChange={onChange} ariaLabel={ariaLabel}>
      <LanguageOptions includeAuto={includeAuto} />
    </Select>
  );
}

export function EngineOptions({ settings }: { settings: Settings }) {
  const engines = Object.values(settings.engines);
  const groups: { label: string; items: typeof engines }[] = [
    { label: t('eng.groupFree'), items: engines.filter((e) => ENGINE_META[e.type].free) },
    { label: t('eng.groupMt'), items: engines.filter((e) => ENGINE_META[e.type].kind === 'mt' && !ENGINE_META[e.type].free) },
    { label: t('eng.groupLlm'), items: engines.filter((e) => ENGINE_META[e.type].kind === 'llm') },
  ];
  return (
    <>
      {groups
        .filter((g) => g.items.length)
        .map((g) => (
          <optgroup key={g.label} label={g.label}>
            {g.items.map((e) => (
              <option key={e.id} value={e.id}>
                {e.name}
                {isEngineReady(e) ? '' : ` · ${t('common.needsKey')}`}
              </option>
            ))}
          </optgroup>
        ))}
    </>
  );
}

export function ExpertOptions({ settings }: { settings: Settings }) {
  const ui = getLocale();
  const experts = getAllExperts(settings);
  const groups = [
    { label: t('exp.academic'), items: experts.filter((e) => e.builtin && e.group === 'academic') },
    { label: t('exp.general'), items: experts.filter((e) => e.builtin && e.group !== 'academic') },
    { label: t('exp.custom'), items: experts.filter((e) => !e.builtin) },
  ];
  return (
    <>
      {groups
        .filter((g) => g.items.length)
        .map((g) => (
          <optgroup key={g.label} label={g.label}>
            {g.items.map((e) => (
              <option key={e.id} value={e.id}>
                {e.icon ? `${e.icon} ` : ''}
                {expertName(e, ui)}
              </option>
            ))}
          </optgroup>
        ))}
    </>
  );
}

/** Yazarken değil, kısa bir duraklamadan veya odaktan çıkınca kaydeden metin kutusu. */
export function TextInput({
  value,
  onCommit,
  placeholder,
  type = 'text',
  className = '',
  monospace,
  ariaLabel,
}: {
  value: string;
  onCommit: (v: string) => void;
  placeholder?: string;
  type?: 'text' | 'password' | 'url' | 'number';
  className?: string;
  monospace?: boolean;
  ariaLabel?: string;
}) {
  const [local, setLocal] = useState(value);
  useEffect(() => setLocal(value), [value]);
  const commit = useDebouncedCallback((v: string) => {
    if (v !== value) onCommit(v);
  }, 600);
  return (
    <input
      className={`input ${className}`}
      style={monospace ? { fontFamily: 'ui-monospace, SFMono-Regular, Menlo, monospace', fontSize: 12.5 } : undefined}
      type={type}
      value={local}
      placeholder={placeholder}
      aria-label={ariaLabel}
      spellCheck={false}
      autoComplete="off"
      onChange={(e) => {
        setLocal(e.target.value);
        commit(e.target.value);
      }}
      onBlur={() => local !== value && onCommit(local)}
      onKeyDown={(e) => {
        if (e.key === 'Enter') (e.target as HTMLInputElement).blur();
      }}
    />
  );
}

export function SecretInput({ value, onCommit, placeholder }: { value: string; onCommit: (v: string) => void; placeholder?: string }) {
  const [show, setShow] = useState(false);
  return (
    <div className="input-group">
      <TextInput value={value} onCommit={onCommit} placeholder={placeholder} type={show ? 'text' : 'password'} monospace />
      <button type="button" className="btn btn-icon" onClick={() => setShow((s) => !s)} aria-label="show">
        {show ? '🙈' : '👁'}
      </button>
    </div>
  );
}

export function TextArea({ value, onCommit, rows = 8, placeholder }: { value: string; onCommit: (v: string) => void; rows?: number; placeholder?: string }) {
  const [local, setLocal] = useState(value);
  useEffect(() => setLocal(value), [value]);
  const commit = useDebouncedCallback((v: string) => v !== value && onCommit(v), 800);
  return (
    <textarea
      className="textarea"
      rows={rows}
      value={local}
      placeholder={placeholder}
      spellCheck={false}
      onChange={(e) => {
        setLocal(e.target.value);
        commit(e.target.value);
      }}
      onBlur={() => local !== value && onCommit(local)}
    />
  );
}

/** Ayar satırı: solda başlık ve açıklama, sağda kontrol. */
export function Row({ title, desc, children, vertical }: { title: ReactNode; desc?: ReactNode; children?: ReactNode; vertical?: boolean }) {
  return (
    <div className={`row ${vertical ? 'row-v' : ''}`}>
      <div className="row-text">
        <div className="row-title">{title}</div>
        {desc && <div className="row-desc">{desc}</div>}
      </div>
      {children !== undefined && <div className="row-control">{children}</div>}
    </div>
  );
}

export function Card({ title, desc, children, actions }: { title?: ReactNode; desc?: ReactNode; children: ReactNode; actions?: ReactNode }) {
  return (
    <section className="card">
      {(title || actions) && (
        <header className="card-head">
          <div>
            {title && <h3>{title}</h3>}
            {desc && <p className="card-desc">{desc}</p>}
          </div>
          {actions && <div className="card-actions">{actions}</div>}
        </header>
      )}
      <div className="card-body">{children}</div>
    </section>
  );
}

/** Alan adı / dil gibi dize listelerini düzenler. */
export function TagList({
  items,
  onChange,
  placeholder,
  render,
}: {
  items: string[];
  onChange: (items: string[]) => void;
  placeholder?: string;
  render?: (item: string) => ReactNode;
}) {
  const [draft, setDraft] = useState('');
  const add = () => {
    const v = draft.trim();
    if (!v || items.includes(v)) return setDraft('');
    onChange([...items, v]);
    setDraft('');
  };
  return (
    <div className="taglist">
      <div className="tags">
        {items.length === 0 && <span className="faint small">{t('rul.empty')}</span>}
        {items.map((item) => (
          <span key={item} className="tag">
            {render ? render(item) : item}
            <button type="button" aria-label={t('common.delete')} onClick={() => onChange(items.filter((x) => x !== item))}>
              ×
            </button>
          </span>
        ))}
      </div>
      <div className="input-group">
        <input className="input" value={draft} placeholder={placeholder} onChange={(e) => setDraft(e.target.value)} onKeyDown={(e) => e.key === 'Enter' && add()} />
        <button type="button" className="btn" onClick={add}>
          {t('common.add')}
        </button>
      </div>
    </div>
  );
}

export function LanguageTagList({ items, onChange }: { items: string[]; onChange: (items: string[]) => void }) {
  const [draft, setDraft] = useState('en');
  const ui = getLocale();
  return (
    <div className="taglist">
      <div className="tags">
        {items.length === 0 && <span className="faint small">{t('rul.empty')}</span>}
        {items.map((code) => (
          <span key={code} className="tag">
            {languageName(code, ui)}
            <button type="button" aria-label={t('common.delete')} onClick={() => onChange(items.filter((x) => x !== code))}>
              ×
            </button>
          </span>
        ))}
      </div>
      <div className="input-group">
        <LanguageSelect value={draft} onChange={setDraft} />
        <button type="button" className="btn" onClick={() => !items.includes(draft) && onChange([...items, draft])}>
          {t('common.add')}
        </button>
      </div>
    </div>
  );
}

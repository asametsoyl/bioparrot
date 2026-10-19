import { useMemo, useRef, useState } from 'react';
import { t } from '@/core/i18n';
import type { GlossaryEntry } from '@/core/types';
import { Card, Switch } from '@/ui/components';
import type { SectionProps } from '../App';

/** Tırnaklı alanları destekleyen basit CSV ayrıştırıcı (virgül, noktalı virgül veya sekme). */
export function parseCsv(text: string): string[][] {
  const firstLine = text.split(/\r?\n/, 1)[0] ?? '';
  const sep = firstLine.includes('\t') ? '\t' : firstLine.split(';').length > firstLine.split(',').length ? ';' : ',';
  const rows: string[][] = [];
  let row: string[] = [];
  let field = '';
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (quoted) {
      if (c === '"' && text[i + 1] === '"') {
        field += '"';
        i++;
      } else if (c === '"') quoted = false;
      else field += c;
    } else if (c === '"') quoted = true;
    else if (c === sep) {
      row.push(field);
      field = '';
    } else if (c === '\n' || c === '\r') {
      if (c === '\r' && text[i + 1] === '\n') i++;
      row.push(field);
      rows.push(row);
      row = [];
      field = '';
    } else field += c;
  }
  if (field || row.length) {
    row.push(field);
    rows.push(row);
  }
  return rows.filter((r) => r.some((x) => x.trim()));
}

function toCsv(entries: GlossaryEntry[]): string {
  const esc = (s: string) => (/[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s);
  return ['source,target,case_sensitive', ...entries.map((e) => `${esc(e.source)},${esc(e.target)},${e.caseSensitive ? 1 : 0}`)].join('\n');
}

export function GlossarySection({ settings, update }: SectionProps) {
  const [source, setSource] = useState('');
  const [target, setTarget] = useState('');
  const [query, setQuery] = useState('');
  const [message, setMessage] = useState('');
  const fileRef = useRef<HTMLInputElement>(null);

  const filtered = useMemo(() => {
    const q = query.trim().toLowerCase();
    return q ? settings.glossary.filter((e) => e.source.toLowerCase().includes(q) || e.target.toLowerCase().includes(q)) : settings.glossary;
  }, [settings.glossary, query]);

  const add = () => {
    if (!source.trim() || !target.trim()) return;
    update((d) => {
      d.glossary = d.glossary.filter((e) => e.source.toLowerCase() !== source.trim().toLowerCase());
      d.glossary.unshift({ id: crypto.randomUUID(), source: source.trim(), target: target.trim() });
    });
    setSource('');
    setTarget('');
  };

  const importCsv = async (file: File) => {
    const rows = parseCsv(await file.text());
    const entries = rows
      .filter((r, i) => !(i === 0 && /^source$/i.test(r[0]?.trim() ?? '')))
      .filter((r) => r[0]?.trim() && r[1]?.trim())
      .map((r) => ({ id: crypto.randomUUID(), source: r[0].trim(), target: r[1].trim(), caseSensitive: /^(1|true|yes|evet)$/i.test(r[2]?.trim() ?? '') }));
    await update((d) => {
      const map = new Map(d.glossary.map((e) => [e.source.toLowerCase(), e]));
      for (const e of entries) map.set(e.source.toLowerCase(), e);
      d.glossary = Array.from(map.values());
    });
    setMessage(t('glo.imported', { n: entries.length }));
  };

  const exportCsv = () => {
    const blob = new Blob(['﻿' + toCsv(settings.glossary)], { type: 'text/csv;charset=utf-8' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'bioparrot-glossary.csv';
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  };

  return (
    <>
      <div className="note">{t('glo.intro')}</div>
      <Card
        title={t('glo.count', { n: settings.glossary.length })}
        actions={
          <>
            <button className="btn btn-sm" onClick={() => fileRef.current?.click()}>
              {t('glo.import')}
            </button>
            <button className="btn btn-sm" onClick={exportCsv} disabled={!settings.glossary.length}>
              {t('glo.export')}
            </button>
            <input
              ref={fileRef}
              type="file"
              accept=".csv,.tsv,.txt"
              hidden
              onChange={(e) => {
                const f = e.target.files?.[0];
                if (f) void importCsv(f);
                e.target.value = '';
              }}
            />
          </>
        }
      >
        {message && <div className="note green" style={{ marginTop: 8 }}>{message}</div>}
        <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr auto', gap: 8, padding: '10px 0' }}>
          <input className="input" placeholder={`${t('glo.source')} (ör. extracellular vesicle)`} value={source} onChange={(e) => setSource(e.target.value)} onKeyDown={(e) => e.key === 'Enter' && add()} />
          <input className="input" placeholder={`${t('glo.target')} (ör. hücre dışı vezikül)`} value={target} onChange={(e) => setTarget(e.target.value)} onKeyDown={(e) => e.key === 'Enter' && add()} />
          <button className="btn btn-primary" onClick={add}>
            {t('common.add')}
          </button>
        </div>
        {settings.glossary.length > 8 && (
          <input className="input" placeholder={t('common.search')} value={query} onChange={(e) => setQuery(e.target.value)} style={{ marginBottom: 8 }} />
        )}
        {settings.glossary.length === 0 ? (
          <p className="small faint" style={{ padding: '6px 0 10px' }}>
            {t('glo.empty')}
          </p>
        ) : (
          <table className="table">
            <thead>
              <tr>
                <th>{t('glo.source')}</th>
                <th>{t('glo.target')}</th>
                <th style={{ width: 120 }}>{t('glo.caseSensitive')}</th>
                <th style={{ width: 60 }} />
              </tr>
            </thead>
            <tbody>
              {filtered.slice(0, 500).map((e) => (
                <tr key={e.id}>
                  <td>
                    <input
                      className="input"
                      defaultValue={e.source}
                      onBlur={(ev) => ev.target.value.trim() && ev.target.value !== e.source && update((d) => void (d.glossary.find((x) => x.id === e.id)!.source = ev.target.value.trim()))}
                    />
                  </td>
                  <td>
                    <input
                      className="input"
                      defaultValue={e.target}
                      onBlur={(ev) => ev.target.value.trim() && ev.target.value !== e.target && update((d) => void (d.glossary.find((x) => x.id === e.id)!.target = ev.target.value.trim()))}
                    />
                  </td>
                  <td>
                    <Switch checked={!!e.caseSensitive} onChange={(v) => update((d) => void (d.glossary.find((x) => x.id === e.id)!.caseSensitive = v))} />
                  </td>
                  <td>
                    <button className="btn btn-sm btn-ghost btn-danger" onClick={() => update((d) => void (d.glossary = d.glossary.filter((x) => x.id !== e.id)))}>
                      ✕
                    </button>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        )}
      </Card>
    </>
  );
}

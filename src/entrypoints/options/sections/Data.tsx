import { useEffect, useMemo, useRef, useState } from 'react';
import { browser } from 'wxt/browser';
import { STORAGE_KEYS } from '@/core/constants';
import { createDefaultSettings, mergeSettings } from '@/core/defaults';
import { t } from '@/core/i18n';
import { sendToBackground } from '@/lib/messaging';
import { saveSettings } from '@/lib/settings';
import { readUsage, type UsageData } from '@/lib/usage';
import { Card, Row, Switch } from '@/ui/components';
import type { SectionProps } from '../App';

const nf = new Intl.NumberFormat();

export function DataSection({ settings, update }: SectionProps) {
  const [entries, setEntries] = useState<number | null>(null);
  const [usage, setUsage] = useState<UsageData>({});
  const [withKeys, setWithKeys] = useState(false);
  const [message, setMessage] = useState('');
  const fileRef = useRef<HTMLInputElement>(null);

  const refresh = () => {
    void sendToBackground({ type: 'cacheStats' })
      .then((r) => setEntries(r.entries))
      .catch(() => setEntries(null));
    void readUsage().then(setUsage);
  };
  useEffect(refresh, []);

  const totals = useMemo(() => {
    const cutoff = new Date(Date.now() - 30 * 86_400_000).toISOString().slice(0, 10);
    const out: Record<string, { chars: number; requests: number; input: number; output: number }> = {};
    for (const [day, engines] of Object.entries(usage)) {
      if (day < cutoff) continue;
      for (const [id, s] of Object.entries(engines)) {
        const o = (out[id] ??= { chars: 0, requests: 0, input: 0, output: 0 });
        o.chars += s.chars;
        o.requests += s.requests;
        o.input += s.input;
        o.output += s.output;
      }
    }
    return out;
  }, [usage]);

  const exportSettings = () => {
    const data = structuredClone(settings);
    if (!withKeys) for (const e of Object.values(data.engines)) e.apiKey = '';
    const blob = new Blob([JSON.stringify({ app: 'bioparrot', exportedAt: new Date().toISOString(), settings: data }, null, 2)], { type: 'application/json' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = `bioparrot-settings-${new Date().toISOString().slice(0, 10)}.json`;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  };

  const importSettings = async (file: File) => {
    try {
      const json = JSON.parse(await file.text());
      const incoming = json?.settings ?? json;
      const merged = mergeSettings(createDefaultSettings(), incoming);
      // Dosyada boş bırakılan anahtarlar mevcut anahtarları silmesin.
      for (const [id, e] of Object.entries(merged.engines)) {
        if (!e.apiKey && settings.engines[id]?.apiKey) e.apiKey = settings.engines[id].apiKey;
      }
      await saveSettings(merged);
      setMessage(t('dat.imported'));
    } catch (err) {
      setMessage(t('dat.importFailed', { msg: String((err as Error)?.message ?? err) }));
    }
  };

  return (
    <>
      {message && <div className="note green">{message}</div>}
      <Card title={t('dat.cache')} desc={t('dat.cache.desc')}>
        <Row title={t('dat.cacheEnabled')}>
          <Switch checked={settings.cache.enabled} onChange={(v) => update((d) => void (d.cache.enabled = v))} />
        </Row>
        <Row title={t('dat.cacheTtl')}>
          <input
            className="input"
            type="number"
            min={1}
            max={365}
            style={{ width: 100 }}
            value={settings.cache.ttlDays}
            onChange={(e) => e.target.value && update((d) => void (d.cache.ttlDays = Math.max(1, Number(e.target.value))))}
          />
        </Row>
        <Row title={entries === null ? '…' : t('dat.cacheEntries', { n: nf.format(entries) })}>
          <button
            className="btn"
            onClick={async () => {
              await sendToBackground({ type: 'clearCache' });
              setMessage(t('dat.cleared'));
              refresh();
            }}
          >
            {t('dat.clearCache')}
          </button>
        </Row>
      </Card>

      <Card title={t('dat.usage')}>
        {Object.keys(totals).length === 0 ? (
          <p className="small faint" style={{ padding: '8px 0 10px' }}>
            {t('dat.usageEmpty')}
          </p>
        ) : (
          <table className="table" style={{ margin: '6px 0 8px' }}>
            <thead>
              <tr>
                <th>{t('engine.label')}</th>
                <th className="num">{t('dat.requests')}</th>
                <th className="num">{t('dat.chars')}</th>
                <th className="num">{t('dat.tokens')}</th>
              </tr>
            </thead>
            <tbody>
              {Object.entries(totals).map(([id, s]) => (
                <tr key={id}>
                  <td>{settings.engines[id]?.name ?? id}</td>
                  <td className="num">{nf.format(s.requests)}</td>
                  <td className="num">{nf.format(s.chars)}</td>
                  <td className="num">{s.input || s.output ? `${nf.format(s.input)} / ${nf.format(s.output)}` : '—'}</td>
                </tr>
              ))}
            </tbody>
          </table>
        )}
      </Card>

      <Card title={t('dat.backup')}>
        <Row title={t('dat.export')}>
          <div className="input-group">
            <label className="small muted" style={{ display: 'flex', alignItems: 'center', gap: 6 }}>
              <input type="checkbox" checked={withKeys} onChange={(e) => setWithKeys(e.target.checked)} />
              {t('dat.exportWithKeys')}
            </label>
            <button className="btn" onClick={exportSettings}>
              {t('dat.export')}
            </button>
          </div>
        </Row>
        <Row title={t('dat.import')}>
          <>
            <button className="btn" onClick={() => fileRef.current?.click()}>
              {t('dat.import')}
            </button>
            <input
              ref={fileRef}
              type="file"
              accept="application/json,.json"
              hidden
              onChange={(e) => {
                const f = e.target.files?.[0];
                if (f) void importSettings(f);
                e.target.value = '';
              }}
            />
          </>
        </Row>
        <Row title={t('dat.reset')}>
          <button
            className="btn btn-danger"
            onClick={async () => {
              if (!confirm(t('dat.resetConfirm'))) return;
              await browser.storage.local.remove(STORAGE_KEYS.settings);
              location.reload();
            }}
          >
            {t('common.reset')}
          </button>
        </Row>
      </Card>
    </>
  );
}

import { useState } from 'react';
import { ENGINE_META, OPENAI_PRESETS, getPreset, isEngineReady } from '@/core/engines/meta';
import { t } from '@/core/i18n';
import type { EngineConfig } from '@/core/types';
import { sendToBackground } from '@/lib/messaging';
import { Card, Row, SecretInput, Select, TextArea, TextInput } from '@/ui/components';
import type { SectionProps } from '../App';

function NumberField({ label, value, onCommit, step = 1 }: { label: string; value: number; onCommit: (v: number) => void; step?: number }) {
  return (
    <label>
      {label}
      <input className="input" type="number" min={step < 1 ? 0.1 : 1} step={step} value={value} onChange={(e) => e.target.value && onCommit(Number(e.target.value))} />
    </label>
  );
}

function EngineCard({
  config,
  isDefault,
  open,
  onToggle,
  onChange,
  onMakeDefault,
  onRemove,
  target,
}: {
  config: EngineConfig;
  isDefault: boolean;
  open: boolean;
  onToggle: () => void;
  onChange: (patch: Partial<EngineConfig>) => void;
  onMakeDefault: () => void;
  onRemove?: () => void;
  target: string;
}) {
  const meta = ENGINE_META[config.type];
  const preset = config.type === 'openai' ? getPreset(config.preset) : undefined;
  const ready = isEngineReady(config);
  const [test, setTest] = useState<{ state: 'idle' | 'busy' | 'ok' | 'err'; text?: string }>({ state: 'idle' });
  const models = preset?.models ?? meta.models ?? [];
  const keyUrl = preset?.keyUrl ?? meta.keyUrl;

  const runTest = async () => {
    setTest({ state: 'busy' });
    try {
      const res = await sendToBackground({ type: 'testEngine', engineId: config.id, to: target, text: t('eng.testSample') });
      setTest({ state: 'ok', text: `${res.translation}  ·  ${res.ms} ms` });
    } catch (err) {
      setTest({ state: 'err', text: String((err as Error)?.message ?? err) });
    }
  };

  return (
    <div className={`engine ${isDefault ? 'is-default' : ''}`}>
      <div className="engine-head" onClick={onToggle}>
        <span className={`chev ${open ? 'open' : ''}`}>▶</span>
        <div className="engine-name">
          {config.name}
          {meta.free && <span className="badge badge-green">{t('common.free')}</span>}
          {meta.kind === 'llm' && <span className="badge badge-orange">AI</span>}
          {ready ? <span className="badge">{t('common.ready')}</span> : <span className="badge badge-gray">{t('common.notConfigured')}</span>}
        </div>
        {isDefault ? (
          <span className="badge badge-green">✓ {t('eng.isDefault')}</span>
        ) : (
          <button
            className="btn btn-sm"
            disabled={!ready}
            onClick={(e) => {
              e.stopPropagation();
              onMakeDefault();
            }}
          >
            {t('eng.setDefault')}
          </button>
        )}
      </div>
      {open && (
        <div className="engine-body">
          {config.type === 'anthropic' && <div className="note green" style={{ marginTop: 10 }}>{t('eng.claudeHint')}</div>}
          {preset?.note && <div className="note" style={{ marginTop: 10 }}>{preset.note}</div>}
          {config.type === 'openai' && (
            <>
              <Row title={t('eng.preset')}>
                <Select
                  value={config.preset ?? 'custom'}
                  onChange={(id) => {
                    const p = getPreset(id)!;
                    onChange({ preset: id, baseUrl: p.baseUrl, model: p.model, ...(config.builtin ? {} : { name: p.label }) });
                  }}
                >
                  {OPENAI_PRESETS.map((p) => (
                    <option key={p.id} value={p.id}>
                      {p.label}
                    </option>
                  ))}
                </Select>
              </Row>
              {!config.builtin && (
                <Row title={t('eng.name')}>
                  <TextInput value={config.name} onCommit={(v) => onChange({ name: v || 'OpenAI' })} />
                </Row>
              )}
            </>
          )}
          {meta.fields.includes('apiKey') && (
            <Row
              title={t('eng.apiKey')}
              desc={
                keyUrl ? (
                  <a href={keyUrl} target="_blank" rel="noreferrer">
                    {t('eng.getKey')} ↗
                  </a>
                ) : undefined
              }
            >
              <SecretInput value={config.apiKey ?? ''} onCommit={(v) => onChange({ apiKey: v.trim() })} placeholder={config.type === 'deeplx' ? '(isteğe bağlı)' : 'sk-…'} />
            </Row>
          )}
          {meta.fields.includes('region') && (
            <Row title={t('eng.region')} desc={t('eng.region.desc')}>
              <TextInput value={config.region ?? ''} onCommit={(v) => onChange({ region: v.trim() })} placeholder="westeurope" />
            </Row>
          )}
          {meta.fields.includes('baseUrl') && (
            <Row title={t('eng.baseUrl')}>
              <TextInput
                value={config.baseUrl ?? ''}
                onCommit={(v) => onChange({ baseUrl: v.trim() })}
                monospace
                placeholder={preset?.baseUrl ?? (config.type === 'deeplx' ? 'http://localhost:1188/translate' : config.type === 'microsoft' ? 'https://api.cognitive.microsofttranslator.com' : '')}
              />
            </Row>
          )}
          {meta.fields.includes('model') && (
            <Row title={t('eng.model')}>
              <>
                <input
                  className="input"
                  list={`models-${config.id}`}
                  defaultValue={config.model ?? ''}
                  key={config.model}
                  spellCheck={false}
                  onBlur={(e) => e.target.value !== config.model && onChange({ model: e.target.value.trim() })}
                  onKeyDown={(e) => e.key === 'Enter' && (e.target as HTMLInputElement).blur()}
                />
                <datalist id={`models-${config.id}`}>
                  {models.map((m) => (
                    <option key={m} value={m} />
                  ))}
                </datalist>
              </>
            </Row>
          )}
          {meta.fields.includes('effort') && (
            <Row title={t('eng.effort')} desc={t('eng.effort.desc')}>
              <Select value={config.effort ?? ''} onChange={(v) => onChange({ effort: v as EngineConfig['effort'] })}>
                <option value="">{t('eng.effort.default')}</option>
                <option value="low">{t('eng.effort.low')}</option>
                <option value="medium">{t('eng.effort.medium')}</option>
                <option value="high">{t('eng.effort.high')}</option>
              </Select>
            </Row>
          )}
          {meta.fields.includes('temperature') && (
            <Row title={t('eng.temperature')} desc={t('eng.temperature.desc')}>
              <TextInput
                value={config.temperature === undefined ? '' : String(config.temperature)}
                onCommit={(v) => onChange({ temperature: v.trim() === '' ? undefined : Math.min(2, Math.max(0, Number(v))) })}
                placeholder="0.2"
              />
            </Row>
          )}

          <details className="adv">
            <summary>{t('eng.advanced')}</summary>
            <p className="small muted" style={{ margin: '4px 0 8px' }}>
              {t('eng.limits')}
            </p>
            <div className="grid-4">
              <NumberField label={t('eng.maxConcurrent')} value={config.maxConcurrent} onCommit={(v) => onChange({ maxConcurrent: v })} />
              <NumberField label={t('eng.maxRps')} value={config.maxRps} step={0.5} onCommit={(v) => onChange({ maxRps: v })} />
              <NumberField label={t('eng.batchSize')} value={config.batchSize} onCommit={(v) => onChange({ batchSize: v })} />
              <NumberField label={t('eng.batchChars')} value={config.batchChars} step={100} onCommit={(v) => onChange({ batchChars: v })} />
            </div>
            {meta.fields.includes('extraBody') && (
              <Row title={t('eng.extraBody')} vertical>
                <TextArea rows={3} value={config.extraBody ?? ''} onCommit={(v) => onChange({ extraBody: v })} placeholder='{"reasoning_effort": "minimal"}' />
              </Row>
            )}
            {meta.fields.includes('customHeaders') && (
              <Row title={t('eng.customHeaders')} vertical>
                <TextArea rows={2} value={config.customHeaders ?? ''} onCommit={(v) => onChange({ customHeaders: v })} placeholder='{"X-Title": "BioParrot"}' />
              </Row>
            )}
          </details>

          <div style={{ display: 'flex', gap: 8, marginTop: 12 }}>
            <button className="btn" onClick={runTest} disabled={test.state === 'busy'}>
              {test.state === 'busy' ? (
                <>
                  <span className="spinner" /> {t('common.testing')}
                </>
              ) : (
                t('common.test')
              )}
            </button>
            {onRemove && (
              <button className="btn btn-danger" onClick={onRemove}>
                {t('eng.remove')}
              </button>
            )}
          </div>
          {(test.state === 'ok' || test.state === 'err') && (
            <div className={`test-result ${test.state}`}>
              {test.state === 'ok' ? '✓ ' : '✕ '}
              {test.text}
            </div>
          )}
        </div>
      )}
    </div>
  );
}

export function EnginesSection({ settings, update }: SectionProps) {
  const [open, setOpen] = useState<string | null>(null);
  const engines = Object.values(settings.engines);
  const order = (e: EngineConfig) => (ENGINE_META[e.type].free ? 0 : ENGINE_META[e.type].kind === 'llm' ? 1 : 2);
  const sorted = [...engines].sort((a, b) => order(a) - order(b));

  const change = (id: string, patch: Partial<EngineConfig>) =>
    update((d) => {
      d.engines[id] = { ...d.engines[id], ...patch };
    });

  const addCustom = () => {
    const id = `custom-${Math.random().toString(36).slice(2, 8)}`;
    const p = getPreset('custom')!;
    update((d) => {
      d.engines[id] = {
        id,
        type: 'openai',
        name: p.label,
        preset: 'custom',
        baseUrl: '',
        model: '',
        apiKey: '',
        maxConcurrent: 4,
        maxRps: 4,
        batchSize: 12,
        batchChars: 2400,
      };
    });
    setOpen(id);
  };

  return (
    <>
      <div className="note">{t('eng.intro')}</div>
      <div className="note green">🎓 {t('eng.llmHint')}</div>
      {sorted.map((config) => (
        <EngineCard
          key={config.id}
          config={config}
          target={settings.targetLanguage}
          isDefault={settings.engineId === config.id}
          open={open === config.id}
          onToggle={() => setOpen(open === config.id ? null : config.id)}
          onChange={(patch) => change(config.id, patch)}
          onMakeDefault={() => update((d) => void (d.engineId = config.id))}
          onRemove={
            config.builtin
              ? undefined
              : () =>
                  update((d) => {
                    delete d.engines[config.id];
                    if (d.engineId === config.id) d.engineId = 'google';
                  })
          }
        />
      ))}
      <Card>
        <Row title={t('eng.addCustom')} desc="DeepSeek, OpenRouter, Groq, Mistral, Ollama, LM Studio, Azure OpenAI…">
          <button className="btn btn-primary" onClick={addCustom}>
            + {t('common.add')}
          </button>
        </Row>
      </Card>
    </>
  );
}

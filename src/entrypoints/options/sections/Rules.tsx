import { t } from '@/core/i18n';
import { Card, LanguageTagList, Row, Switch, TagList, TextInput } from '@/ui/components';
import type { SectionProps } from '../App';

export function RulesSection({ settings, update }: SectionProps) {
  const rules = settings.rules;
  return (
    <>
      <Card title={t('rul.alwaysSites')}>
        <div style={{ padding: '10px 0' }}>
          <TagList items={rules.alwaysSites} placeholder={t('rul.sitePlaceholder')} onChange={(items) => update((d) => void (d.rules.alwaysSites = items))} />
        </div>
      </Card>
      <Card title={t('rul.neverSites')}>
        <div style={{ padding: '10px 0' }}>
          <TagList items={rules.neverSites} placeholder={t('rul.sitePlaceholder')} onChange={(items) => update((d) => void (d.rules.neverSites = items))} />
        </div>
      </Card>
      <Card title={t('rul.alwaysLangs')} desc={t('rul.alwaysLangs.desc')}>
        <div style={{ padding: '10px 0' }}>
          <LanguageTagList items={rules.alwaysLanguages} onChange={(items) => update((d) => void (d.rules.alwaysLanguages = items))} />
        </div>
      </Card>
      <Card title={t('rul.neverLangs')} desc={t('rul.neverLangs.desc')}>
        <div style={{ padding: '10px 0' }}>
          <LanguageTagList items={rules.neverLanguages} onChange={(items) => update((d) => void (d.rules.neverLanguages = items))} />
        </div>
      </Card>
      <Card
        title={t('rul.custom')}
        desc={t('rul.custom.desc')}
        actions={
          <button
            className="btn btn-sm btn-primary"
            onClick={() =>
              update((d) => {
                d.rules.custom.push({ id: crypto.randomUUID(), pattern: '', include: '', exclude: '', enabled: true });
              })
            }
          >
            + {t('rul.addRule')}
          </button>
        }
      >
        {rules.custom.length === 0 && <p className="small faint" style={{ padding: '8px 0' }}>{t('rul.empty')}</p>}
        {rules.custom.map((rule) => {
          const patch = (p: Partial<typeof rule>) =>
            update((d) => {
              d.rules.custom = d.rules.custom.map((r) => (r.id === rule.id ? { ...r, ...p } : r));
            });
          return (
            <div key={rule.id} style={{ borderBottom: '1px solid var(--border)', padding: '10px 0' }}>
              <Row title={t('rul.pattern')}>
                <div className="input-group">
                  <Switch checked={rule.enabled} onChange={(v) => patch({ enabled: v })} />
                  <TextInput value={rule.pattern} onCommit={(v) => patch({ pattern: v.trim() })} placeholder={t('rul.sitePlaceholder')} />
                  <button className="btn btn-ghost btn-danger" onClick={() => update((d) => void (d.rules.custom = d.rules.custom.filter((r) => r.id !== rule.id)))}>
                    ✕
                  </button>
                </div>
              </Row>
              <Row title={t('rul.include')}>
                <TextInput value={rule.include ?? ''} onCommit={(v) => patch({ include: v })} placeholder="article, .post-content" monospace />
              </Row>
              <Row title={t('rul.exclude')}>
                <TextInput value={rule.exclude ?? ''} onCommit={(v) => patch({ exclude: v })} placeholder=".comments, .sidebar" monospace />
              </Row>
            </div>
          );
        })}
      </Card>
    </>
  );
}

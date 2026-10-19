import { isLlmEngine } from '@/core/engines/meta';
import { expertName, getAllExperts } from '@/core/experts';
import { getLocale, t } from '@/core/i18n';
import { Card, Row, Select, Switch } from '@/ui/components';
import type { SectionProps } from '../App';

export function AcademicSection({ settings, update, navigate }: SectionProps) {
  const ui = getLocale();
  const academicExperts = getAllExperts(settings).filter((e) => e.group === 'academic' || !e.builtin);
  const llm = isLlmEngine(settings.engines[settings.engineId]);
  return (
    <>
      <div className="note green">{t('acad.intro')}</div>
      {!llm && (
        <div className="note">
          ✨ {t('eng.llmHint')}{' '}
          <a href="#engines" onClick={() => navigate('engines')}>
            {t('opt.nav.engines')} →
          </a>
        </div>
      )}
      <Card>
        <Row title={t('acad.autoExpert')} desc={t('acad.autoExpert.desc')}>
          <Switch checked={settings.academic.autoExpert} onChange={(v) => update((d) => void (d.academic.autoExpert = v))} />
        </Row>
        <Row title={t('acad.expert')} desc={t('acad.expert.desc')}>
          <Select value={settings.academic.expertId} onChange={(v) => update((d) => void (d.academic.expertId = v))}>
            {academicExperts.map((e) => (
              <option key={e.id} value={e.id}>
                {e.icon ? `${e.icon} ` : ''}
                {expertName(e, ui)}
              </option>
            ))}
          </Select>
        </Row>
        <Row title={t('acad.termConsistency')} desc={t('acad.termConsistency.desc')}>
          <Switch checked={settings.academic.termConsistency} onChange={(v) => update((d) => void (d.academic.termConsistency = v))} />
        </Row>
        <Row title={t('acad.keepTerms')} desc={t('acad.keepTerms.desc')}>
          <Switch checked={settings.academic.keepTerms} onChange={(v) => update((d) => void (d.academic.keepTerms = v))} />
        </Row>
        <Row title={t('acad.skipReferences')} desc={t('acad.skipReferences.desc')}>
          <Switch checked={settings.academic.skipReferences} onChange={(v) => update((d) => void (d.academic.skipReferences = v))} />
        </Row>
      </Card>
      <Card title={t('acad.protected')}>
        <p className="small muted" style={{ padding: '6px 0 8px' }}>
          {t('acad.protected.desc')}
        </p>
      </Card>
      <Card title={t('acad.tips')}>
        <ul className="small muted" style={{ margin: '6px 0 10px', paddingLeft: 18, lineHeight: 1.8 }}>
          <li>
            {t('acad.tip.glossary')}{' '}
            <a href="#glossary" onClick={() => navigate('glossary')}>
              {t('opt.nav.glossary')} →
            </a>
          </li>
          <li>{t('acad.tip.pdf')}</li>
          <li>
            {t('acad.tip.writing')}{' '}
            <a href="#interaction" onClick={() => navigate('interaction')}>
              {t('opt.nav.interaction')} →
            </a>
          </li>
        </ul>
      </Card>
    </>
  );
}

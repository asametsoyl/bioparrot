import { t } from '@/core/i18n';
import { Card, EngineOptions, ExpertOptions, LanguageSelect, Row, Segmented, Select, Switch } from '@/ui/components';
import type { SectionProps } from '../App';

export function GeneralSection({ settings, update, navigate }: SectionProps) {
  return (
    <>
      <Card title={t('gen.languages')}>
        <Row title={t('lang.target')} desc={t('gen.target.desc')}>
          <LanguageSelect value={settings.targetLanguage} onChange={(v) => update((d) => void (d.targetLanguage = v))} />
        </Row>
        <Row title={t('lang.source')} desc={t('gen.source.desc')}>
          <LanguageSelect includeAuto value={settings.sourceLanguage} onChange={(v) => update((d) => void (d.sourceLanguage = v))} />
        </Row>
      </Card>

      <Card title={t('gen.services')}>
        <Row
          title={t('engine.label')}
          desc={
            <>
              {t('gen.engine.desc')}{' '}
              <a href="#engines" onClick={() => navigate('engines')}>
                {t('opt.nav.engines')} →
              </a>
            </>
          }
        >
          <Select value={settings.engineId} onChange={(v) => update((d) => void (d.engineId = v))}>
            <EngineOptions settings={settings} />
          </Select>
        </Row>
        <Row title={t('expert.label')} desc={t('gen.expert.desc')}>
          <Select value={settings.expertId} onChange={(v) => update((d) => void (d.expertId = v))}>
            <ExpertOptions settings={settings} />
          </Select>
        </Row>
        <Row title={t('mode.label')}>
          <Segmented
            value={settings.displayMode}
            onChange={(v) => update((d) => void (d.displayMode = v))}
            options={[
              { value: 'bilingual', label: t('mode.bilingual') },
              { value: 'translation', label: t('mode.translation') },
            ]}
          />
        </Row>
      </Card>

      <Card title={t('gen.page')}>
        <Row title={t('gen.autoAll')} desc={t('gen.autoAll.desc')}>
          <Switch checked={settings.page.autoTranslateAll} onChange={(v) => update((d) => void (d.page.autoTranslateAll = v))} />
        </Row>
        <Row title={t('gen.scope')} desc={t('gen.scope.desc')}>
          <Segmented
            value={settings.page.scope}
            onChange={(v) => update((d) => void (d.page.scope = v))}
            options={[
              { value: 'main', label: t('gen.scope.main') },
              { value: 'full', label: t('gen.scope.full') },
            ]}
          />
        </Row>
        <Row title={t('gen.lazy')} desc={t('gen.lazy.desc')}>
          <Switch checked={settings.page.lazy} onChange={(v) => update((d) => void (d.page.lazy = v))} />
        </Row>
        <Row title={t('gen.contextAware')} desc={t('gen.contextAware.desc')}>
          <Switch checked={settings.page.contextAware} onChange={(v) => update((d) => void (d.page.contextAware = v))} />
        </Row>
        <Row title={t('gen.inlineShort')} desc={t('gen.inlineShort.desc')}>
          <Switch checked={settings.page.inlineShort} onChange={(v) => update((d) => void (d.page.inlineShort = v))} />
        </Row>
        <Row title={t('gen.skipCode')}>
          <Switch checked={settings.page.skipCode} onChange={(v) => update((d) => void (d.page.skipCode = v))} />
        </Row>
        <Row title={t('gen.iframes')}>
          <Switch checked={settings.page.translateIframes} onChange={(v) => update((d) => void (d.page.translateIframes = v))} />
        </Row>
      </Card>

      <Card title={t('gen.interface')}>
        <Row title={t('gen.uiLanguage')}>
          <Select value={settings.uiLanguage} onChange={(v) => update((d) => void (d.uiLanguage = v))}>
            <option value="auto">{t('common.auto')}</option>
            <option value="tr">Türkçe</option>
            <option value="en">English</option>
          </Select>
        </Row>
        <Row title={t('gen.theme')}>
          <Segmented
            value={settings.theme}
            onChange={(v) => update((d) => void (d.theme = v))}
            options={[
              { value: 'system', label: t('common.system') },
              { value: 'light', label: t('common.light') },
              { value: 'dark', label: t('common.dark') },
            ]}
          />
        </Row>
      </Card>
    </>
  );
}

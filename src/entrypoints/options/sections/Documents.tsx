import { browser } from 'wxt/browser';
import { t } from '@/core/i18n';
import { Card, Row, Switch } from '@/ui/components';
import type { SectionProps } from '../App';

export function DocumentsSection({ settings, update }: SectionProps) {
  return (
    <>
      <div className="note green">{t('doc.intro')}</div>
      <Card>
        <Row title={t('pop.pdf')} desc={t('acad.tip.pdf')}>
          <button className="btn btn-primary" onClick={() => void browser.tabs.create({ url: browser.runtime.getURL('/pdf.html') })}>
            {t('doc.openPdf')}
          </button>
        </Row>
        <Row title={t('doc.autoOpen')} desc={t('doc.autoOpen.desc')}>
          <Switch checked={settings.pdf.autoOpen} onChange={(v) => update((d) => void (d.pdf.autoOpen = v))} />
        </Row>
      </Card>
    </>
  );
}

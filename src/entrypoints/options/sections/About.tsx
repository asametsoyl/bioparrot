import { browser } from 'wxt/browser';
import { t } from '@/core/i18n';
import { Card, Logo } from '@/ui/components';
import type { SectionProps } from '../App';

export const REPO_URL = 'https://github.com/asametsoyl/bioparrot';

export function AboutSection(_props: SectionProps) {
  const version = browser.runtime.getManifest().version;
  return (
    <Card>
      <div style={{ display: 'flex', flexDirection: 'column', gap: 14, padding: '14px 0' }}>
        <Logo size={56} sub={t('abt.version', { v: version })} />
        <p>{t('abt.desc')}</p>
        <p className="small muted">{t('abt.privacy')}</p>
        <p className="small">
          <a href={REPO_URL} target="_blank" rel="noreferrer">
            {t('abt.source')} ↗
          </a>{' '}
          · MIT
        </p>
      </div>
    </Card>
  );
}

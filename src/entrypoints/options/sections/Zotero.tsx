import { useState } from 'react';
import { t } from '@/core/i18n';
import { Card } from '@/ui/components';
import type { SectionProps } from '../App';

export function ZoteroSection({ settings }: SectionProps) {
  const [copied, setCopied] = useState(false);
  const copy = async () => {
    await navigator.clipboard.writeText(JSON.stringify({ app: 'bioparrot', exportedAt: new Date().toISOString(), settings }));
    setCopied(true);
  };
  return (
    <>
      <div className="note green">🦜 {t('zot.intro')}</div>
      <Card title={t('zot.featuresTitle')}>
        <ul className="small muted" style={{ margin: '6px 0 10px', paddingLeft: 18, lineHeight: 1.8 }}>
          <li>{t('zot.f1')}</li>
          <li>{t('zot.f2')}</li>
          <li>{t('zot.f3')}</li>
          <li>{t('zot.f4')}</li>
        </ul>
      </Card>
      <Card title={t('zot.installTitle')}>
        <ol className="small" style={{ margin: '6px 0 12px', paddingLeft: 18, lineHeight: 1.9 }}>
          <li>{t('zot.step1')}</li>
          <li>{t('zot.step2')}</li>
          <li>{t('zot.step3')}</li>
        </ol>
        <div style={{ display: 'flex', alignItems: 'center', gap: 10, paddingBottom: 12 }}>
          <button className="btn btn-primary" onClick={() => void copy()}>
            {t('zot.copySettings')}
          </button>
          {copied && <span className="small ok-text">✓ {t('zot.copied')}</span>}
        </div>
      </Card>
    </>
  );
}

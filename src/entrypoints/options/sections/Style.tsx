import { useEffect } from 'react';
import { t, type MessageKey } from '@/core/i18n';
import type { StyleId } from '@/core/types';
import { PAGE_CSS } from '@/content/styles';
import { Card, Row, Segmented } from '@/ui/components';
import type { SectionProps } from '../App';

const STYLES: StyleId[] = ['underline', 'dashed', 'wavy', 'highlight', 'marker', 'none', 'italic', 'bold', 'faded', 'blur', 'border', 'quote', 'card'];
const COLORS = ['#22c55e', '#0ea5e9', '#8b5cf6', '#f97316', '#ef4444', '#eab308', '#64748b'];
/** Sayfaya eklenen gerçek öğe <font> olduğu için önizlemede de aynısı kullanılır. */
const Font = 'font' as unknown as 'span';

function Sample({ styleId, inline, text }: { styleId: StyleId; inline?: boolean; text: string }) {
  return (
    <Font className={`bp-wrap ${inline ? 'bp-inline' : 'bp-block'} bp-style-${styleId}`}>
      {inline && <Font className="bp-sp"> </Font>}
      <Font className="bp-text">{text}</Font>
    </Font>
  );
}

export function StyleSection({ settings, update }: SectionProps) {
  const style = settings.style;

  useEffect(() => {
    let el = document.getElementById('bp-preview-css');
    if (!el) {
      el = document.createElement('style');
      el.id = 'bp-preview-css';
      el.textContent = PAGE_CSS;
      document.head.appendChild(el);
    }
  }, []);

  const vars = {
    '--bp-color': style.color,
    '--bp-font-scale': String(style.fontScale),
    ...(style.textColor ? { '--bp-text-color': style.textColor } : {}),
  } as React.CSSProperties;

  return (
    <div style={vars}>
      <Card title={t('sty.preview')}>
        <div className="preview" style={{ margin: '8px 0 12px' }}>
          {t('sty.previewOriginal')}
          <Sample styleId={style.id} text={t('sty.previewTranslation')} />
        </div>
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

      <Card title={t('sty.style')}>
        <div className="style-grid" style={{ padding: '8px 0 10px' }}>
          {STYLES.map((id) => (
            <button key={id} className={`style-opt ${style.id === id ? 'on' : ''}`} onClick={() => update((d) => void (d.style.id = id))}>
              <span className="name">{t(`style.${id}` as MessageKey)}</span>
              <span className="sample">
                Cell membrane
                <Sample styleId={id} text="Hücre zarı" />
              </span>
            </button>
          ))}
        </div>
      </Card>

      <Card>
        <Row title={t('sty.color')}>
          <div className="swatches">
            {COLORS.map((c) => (
              <span key={c} className={`swatch ${style.color === c ? 'on' : ''}`} style={{ background: c }} onClick={() => update((d) => void (d.style.color = c))} />
            ))}
            <input type="color" value={style.color} onChange={(e) => update((d) => void (d.style.color = e.target.value))} />
          </div>
        </Row>
        <Row title={t('sty.textColor')}>
          <div className="swatches">
            <button className={`btn btn-sm ${!style.textColor ? 'btn-primary' : ''}`} onClick={() => update((d) => void (d.style.textColor = ''))}>
              {t('sty.inherit')}
            </button>
            <input type="color" value={style.textColor || '#334155'} onChange={(e) => update((d) => void (d.style.textColor = e.target.value))} />
          </div>
        </Row>
        <Row title={t('sty.fontScale')} desc={`${Math.round(style.fontScale * 100)}%`}>
          <input type="range" min={0.75} max={1.3} step={0.05} value={style.fontScale} onChange={(e) => update((d) => void (d.style.fontScale = Number(e.target.value)))} />
        </Row>
      </Card>
    </div>
  );
}

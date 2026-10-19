import { useEffect, useState } from 'react';
import { browser } from 'wxt/browser';
import { t } from '@/core/i18n';
import type { HoverKey, InputTrigger, SelectionDisplay, SelectionMode } from '@/core/types';
import { Card, EngineOptions, ExpertOptions, LanguageSelect, Row, Segmented, Select, Switch, TagList } from '@/ui/components';
import type { SectionProps } from '../App';

const isMac = navigator.platform.toLowerCase().includes('mac');

export function InteractionSection({ settings, update }: SectionProps) {
  const [commands, setCommands] = useState<{ name?: string; description?: string; shortcut?: string }[]>([]);
  useEffect(() => {
    void browser.commands.getAll().then(setCommands).catch(() => {});
  }, []);

  const keyLabel = (k: HoverKey) =>
    k === 'off' ? t('common.off') : k === 'Alt' ? (isMac ? 'Option (⌥)' : 'Alt') : k === 'Meta' ? (isMac ? 'Command (⌘)' : 'Win') : k === 'Control' ? 'Ctrl' : 'Shift';

  const openShortcuts = () => {
    const firefox = navigator.userAgent.includes('Firefox');
    void browser.tabs.create({ url: firefox ? 'about:addons' : 'chrome://extensions/shortcuts' });
  };

  return (
    <>
      <Card title={t('int.hover')} desc={t('int.hover.desc')}>
        <Row title={t('int.hoverKey')}>
          <Select value={settings.hover.key} onChange={(v: HoverKey) => update((d) => void (d.hover.key = v))}>
            {(['Shift', 'Alt', 'Control', 'Meta', 'off'] as HoverKey[]).map((k) => (
              <option key={k} value={k}>
                {keyLabel(k)}
              </option>
            ))}
          </Select>
        </Row>
      </Card>

      <Card title={t('int.selection')} desc={t('int.selection.desc')}>
        <Row title={t('int.selectionDisplay')}>
          <Segmented<SelectionDisplay>
            value={settings.selection.display}
            onChange={(v) => update((d) => void (d.selection.display = v))}
            options={[
              { value: 'inline', label: t('int.selection.inline') },
              { value: 'card', label: t('int.selection.card') },
            ]}
          />
        </Row>
        <Row title={t('int.selectionMode')}>
          <Segmented<SelectionMode>
            value={settings.selection.mode}
            onChange={(v) => update((d) => void (d.selection.mode = v))}
            options={[
              { value: 'icon', label: t('int.selection.icon') },
              { value: 'direct', label: t('int.selection.direct') },
              { value: 'off', label: t('common.off') },
            ]}
          />
        </Row>
        <Row title={t('int.selectionInInputs')}>
          <Switch checked={settings.selection.inInputs} onChange={(v) => update((d) => void (d.selection.inInputs = v))} />
        </Row>
        <Row title={t('engine.label')}>
          <Select value={settings.selection.engineId} onChange={(v) => update((d) => void (d.selection.engineId = v))}>
            <option value="">{t('int.sameAsDefault')}</option>
            <EngineOptions settings={settings} />
          </Select>
        </Row>
      </Card>

      <Card title={t('int.input')} desc={t('int.input.desc')}>
        <Row title={t('int.inputTrigger')}>
          <Segmented<InputTrigger>
            value={settings.input.trigger}
            onChange={(v) => update((d) => void (d.input.trigger = v))}
            options={[
              { value: 'triple-space', label: t('int.input.triple') },
              { value: 'double-space', label: t('int.input.double') },
              { value: 'off', label: t('common.off') },
            ]}
          />
        </Row>
        <Row title={t('int.inputTarget')}>
          <LanguageSelect value={settings.input.targetLanguage} onChange={(v) => update((d) => void (d.input.targetLanguage = v))} />
        </Row>
        <Row title={t('int.smartReverse')}>
          <Switch checked={settings.input.smartReverse} onChange={(v) => update((d) => void (d.input.smartReverse = v))} />
        </Row>
        <Row title={t('int.inputEngine')}>
          <Select value={settings.input.engineId} onChange={(v) => update((d) => void (d.input.engineId = v))}>
            <option value="">{t('int.sameAsDefault')}</option>
            <EngineOptions settings={settings} />
          </Select>
        </Row>
        <Row title={t('int.inputExpert')} desc={t('int.inputExpert.desc')}>
          <Select value={settings.input.expertId} onChange={(v) => update((d) => void (d.input.expertId = v))}>
            <option value="">{t('int.sameAsDefault')}</option>
            <ExpertOptions settings={settings} />
          </Select>
        </Row>
      </Card>

      <Card title={t('int.floating')}>
        <Row title={t('int.floatingEnabled')}>
          <Switch checked={settings.floatingBall.enabled} onChange={(v) => update((d) => void (d.floatingBall.enabled = v))} />
        </Row>
        <Row title={t('int.floatingSide')}>
          <Segmented
            value={settings.floatingBall.side}
            onChange={(v) => update((d) => void (d.floatingBall.side = v))}
            options={[
              { value: 'left', label: t('common.left') },
              { value: 'right', label: t('common.right') },
            ]}
          />
        </Row>
        <Row title={t('int.hiddenSites')} vertical>
          <TagList items={settings.floatingBall.hiddenSites} placeholder={t('rul.sitePlaceholder')} onChange={(items) => update((d) => void (d.floatingBall.hiddenSites = items))} />
        </Row>
      </Card>

      <Card
        title={t('int.shortcuts')}
        desc={t('int.shortcuts.desc')}
        actions={
          <button className="btn btn-sm" onClick={openShortcuts}>
            {t('int.shortcutsEdit')}
          </button>
        }
      >
        <div className="kv" style={{ padding: '10px 0' }}>
          {commands
            .filter((c) => c.name && !c.name.startsWith('_'))
            .map((c) => (
              <div key={c.name} style={{ display: 'contents' }}>
                <span className="small">{c.description}</span>
                {c.shortcut ? <kbd>{c.shortcut}</kbd> : <span className="small faint">{t('int.notSet')}</span>}
              </div>
            ))}
        </div>
      </Card>
    </>
  );
}

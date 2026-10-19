import { useState } from 'react';
import { BUILTIN_EXPERTS, expertDescription, expertName, getAllExperts } from '@/core/experts';
import { getLocale, t } from '@/core/i18n';
import type { Expert } from '@/core/types';
import { Card, TextArea, TextInput } from '@/ui/components';
import type { SectionProps } from '../App';

export function ExpertsSection({ settings, update }: SectionProps) {
  const ui = getLocale();
  const [editing, setEditing] = useState<string | null>(null);
  const experts = getAllExperts(settings);
  const groups: { title: string; items: Expert[] }[] = [
    { title: t('exp.academic'), items: experts.filter((e) => e.builtin && e.group === 'academic') },
    { title: t('exp.general'), items: experts.filter((e) => e.builtin && e.group !== 'academic') },
    { title: t('exp.custom'), items: experts.filter((e) => !e.builtin) },
  ];

  const saveBuiltin = (id: string, patch: Partial<Pick<Expert, 'systemPrompt' | 'temperature'>>) =>
    update((d) => {
      const original = BUILTIN_EXPERTS.find((e) => e.id === id)!;
      const next = { ...(d.expertOverrides[id] ?? {}), ...patch };
      if (next.systemPrompt === original.systemPrompt) delete next.systemPrompt;
      if (next.temperature === original.temperature) delete next.temperature;
      if (Object.keys(next).length) d.expertOverrides[id] = next;
      else delete d.expertOverrides[id];
    });

  const saveCustom = (id: string, patch: Partial<Expert>) =>
    update((d) => {
      d.experts = d.experts.map((e) => (e.id === id ? { ...e, ...patch } : e));
    });

  const addExpert = () => {
    const id = `expert-${Math.random().toString(36).slice(2, 8)}`;
    const base = BUILTIN_EXPERTS.find((e) => e.id === 'academic')!;
    update((d) => {
      d.experts.push({
        id,
        name: t('exp.new'),
        description: '',
        icon: '🧪',
        group: 'academic',
        temperature: 0.2,
        systemPrompt: base.systemPrompt,
      });
    });
    setEditing(id);
  };

  return (
    <>
      {groups.map((g) =>
        g.items.length || g.title === t('exp.custom') ? (
          <Card
            key={g.title}
            title={g.title}
            actions={
              g.title === t('exp.custom') ? (
                <button className="btn btn-sm btn-primary" onClick={addExpert}>
                  + {t('exp.new')}
                </button>
              ) : undefined
            }
          >
            {g.items.length === 0 && <p className="small faint" style={{ padding: '8px 0' }}>{t('rul.empty')}</p>}
            {g.items.map((e) => {
              const modified = e.builtin && !!settings.expertOverrides[e.id];
              const isDefault = settings.expertId === e.id;
              const open = editing === e.id;
              return (
                <div key={e.id} className="expert">
                  <div className="expert-ico">{e.icon ?? '🧠'}</div>
                  <div className="expert-main">
                    <div className="expert-title">
                      {expertName(e, ui)}
                      {isDefault && <span className="badge badge-green">✓ {t('common.default')}</span>}
                      {settings.academic.expertId === e.id && <span className="badge">🎓</span>}
                      {modified && <span className="badge badge-orange">{t('exp.modified')}</span>}
                    </div>
                    <div className="small muted">{expertDescription(e, ui)}</div>
                    {open && (
                      <div className="expert-editor">
                        {!e.builtin && (
                          <div style={{ display: 'grid', gridTemplateColumns: '80px 1fr', gap: 8 }}>
                            <TextInput value={e.icon ?? ''} onCommit={(v) => saveCustom(e.id, { icon: v })} ariaLabel="icon" />
                            <TextInput value={e.name} onCommit={(v) => saveCustom(e.id, { name: v })} ariaLabel={t('exp.name')} />
                          </div>
                        )}
                        {!e.builtin && (
                          <TextInput value={e.description ?? ''} onCommit={(v) => saveCustom(e.id, { description: v })} placeholder={t('exp.description')} />
                        )}
                        <label className="small muted">{t('exp.prompt')}</label>
                        <TextArea
                          rows={12}
                          value={e.systemPrompt}
                          onCommit={(v) => (e.builtin ? saveBuiltin(e.id, { systemPrompt: v }) : saveCustom(e.id, { systemPrompt: v }))}
                        />
                        <span className="small faint">{t('exp.prompt.desc')}</span>
                        <div className="input-group" style={{ maxWidth: 260 }}>
                          <span className="small muted">{t('exp.temperature')}</span>
                          <TextInput
                            value={e.temperature === undefined ? '' : String(e.temperature)}
                            onCommit={(v) => {
                              const temperature = v.trim() === '' ? undefined : Math.min(2, Math.max(0, Number(v)));
                              if (e.builtin) saveBuiltin(e.id, { temperature });
                              else saveCustom(e.id, { temperature });
                            }}
                          />
                        </div>
                      </div>
                    )}
                  </div>
                  <div className="expert-actions">
                    {!isDefault && (
                      <button className="btn btn-sm" onClick={() => update((d) => void (d.expertId = e.id))}>
                        {t('exp.makeDefault')}
                      </button>
                    )}
                    <button className="btn btn-sm btn-ghost" onClick={() => setEditing(open ? null : e.id)}>
                      {open ? t('common.close') : t('common.edit')}
                    </button>
                    {modified && (
                      <button
                        className="btn btn-sm btn-ghost"
                        onClick={() =>
                          update((d) => {
                            delete d.expertOverrides[e.id];
                          })
                        }
                      >
                        {t('exp.restore')}
                      </button>
                    )}
                    {!e.builtin && (
                      <button
                        className="btn btn-sm btn-ghost btn-danger"
                        onClick={() =>
                          update((d) => {
                            d.experts = d.experts.filter((x) => x.id !== e.id);
                            if (d.expertId === e.id) d.expertId = 'general';
                            if (d.academic.expertId === e.id) d.academic.expertId = 'academic';
                          })
                        }
                      >
                        {t('common.delete')}
                      </button>
                    )}
                  </div>
                </div>
              );
            })}
          </Card>
        ) : null,
      )}
    </>
  );
}

import { useEffect, useState, type ComponentType } from 'react';
import { t, type MessageKey } from '@/core/i18n';
import type { Settings } from '@/core/types';
import { Logo } from '@/ui/components';
import { useSettings, useTheme, type UpdateSettings } from '@/ui/hooks';
import { AboutSection } from './sections/About';
import { AcademicSection } from './sections/Academic';
import { DataSection } from './sections/Data';
import { DocumentsSection } from './sections/Documents';
import { EnginesSection } from './sections/Engines';
import { ExpertsSection } from './sections/Experts';
import { GeneralSection } from './sections/General';
import { GlossarySection } from './sections/Glossary';
import { InteractionSection } from './sections/Interaction';
import { RulesSection } from './sections/Rules';
import { StyleSection } from './sections/Style';
import { ZoteroSection } from './sections/Zotero';

export interface SectionProps {
  settings: Settings;
  update: UpdateSettings;
  navigate: (id: string) => void;
}

const SECTIONS: { id: string; icon: string; label: MessageKey; Component: ComponentType<SectionProps> }[] = [
  { id: 'general', icon: '⚙️', label: 'opt.nav.general', Component: GeneralSection },
  { id: 'engines', icon: '🔌', label: 'opt.nav.engines', Component: EnginesSection },
  { id: 'academic', icon: '🎓', label: 'opt.nav.academic', Component: AcademicSection },
  { id: 'experts', icon: '🧠', label: 'opt.nav.experts', Component: ExpertsSection },
  { id: 'glossary', icon: '📚', label: 'opt.nav.glossary', Component: GlossarySection },
  { id: 'style', icon: '🎨', label: 'opt.nav.style', Component: StyleSection },
  { id: 'rules', icon: '🌐', label: 'opt.nav.rules', Component: RulesSection },
  { id: 'interaction', icon: '🖱️', label: 'opt.nav.interaction', Component: InteractionSection },
  { id: 'documents', icon: '📄', label: 'opt.nav.documents', Component: DocumentsSection },
  { id: 'zotero', icon: '🦜', label: 'opt.nav.zotero', Component: ZoteroSection },
  { id: 'data', icon: '💾', label: 'opt.nav.data', Component: DataSection },
  { id: 'about', icon: 'ℹ️', label: 'opt.nav.about', Component: AboutSection },
];

function readHash(): { section: string; welcome: boolean } {
  const raw = location.hash.replace(/^#/, '');
  if (raw === 'welcome') return { section: 'general', welcome: true };
  return { section: SECTIONS.some((s) => s.id === raw) ? raw : 'general', welcome: false };
}

export default function App() {
  const { settings, update } = useSettings();
  useTheme(settings?.theme);
  const [{ section, welcome }, setRoute] = useState(readHash);
  const [showWelcome, setShowWelcome] = useState(welcome);

  useEffect(() => {
    const onHash = () => setRoute(readHash());
    window.addEventListener('hashchange', onHash);
    return () => window.removeEventListener('hashchange', onHash);
  }, []);

  useEffect(() => {
    document.title = `BioParrot – ${t(SECTIONS.find((s) => s.id === section)!.label)}`;
    window.scrollTo(0, 0);
  }, [section, settings?.uiLanguage]);

  if (!settings) return null;
  const navigate = (id: string) => {
    location.hash = id;
  };
  const Current = SECTIONS.find((s) => s.id === section)!.Component;

  return (
    <div className="layout">
      <aside className="sidebar">
        <div className="sidebar-head">
          <Logo size={34} sub={t('app.tagline')} />
        </div>
        <nav>
          {SECTIONS.map((s) => (
            <a key={s.id} href={`#${s.id}`} className={`nav-item ${s.id === section ? 'active' : ''}`}>
              <span className="nav-ico">{s.icon}</span>
              <span>{t(s.label)}</span>
            </a>
          ))}
        </nav>
      </aside>
      <main className="content">
        {showWelcome && section === 'general' && (
          <div className="welcome">
            <div>
              <h2>{t('wel.title')}</h2>
              <p className="muted">{t('wel.lead')}</p>
              <ol>
                <li>{t('wel.step1')}</li>
                <li>
                  {t('wel.step2')}{' '}
                  <a href="#engines" onClick={() => setShowWelcome(false)}>
                    → {t('opt.nav.engines')}
                  </a>
                </li>
                <li>{t('wel.step3')}</li>
              </ol>
            </div>
            <button className="btn btn-primary" onClick={() => setShowWelcome(false)}>
              {t('wel.dismiss')}
            </button>
          </div>
        )}
        <h1 className="page-title">{t(SECTIONS.find((s) => s.id === section)!.label)}</h1>
        <Current settings={settings} update={update} navigate={navigate} />
      </main>
    </div>
  );
}

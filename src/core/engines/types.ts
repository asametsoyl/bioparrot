import type { EngineConfig, EngineType, Expert, GlossaryEntry, Settings, TranslationContext } from '../types';

export interface EngineInput {
  /** Motora gönderilecek metinler (zengin biçimde ise etiketli). */
  texts: string[];
  from: string;
  to: string;
  rich: boolean;
  expert: Expert;
  context?: TranslationContext;
  /** Bu metinlerde geçen sözlük girdileri. */
  glossary: GlossaryEntry[];
  academic?: Pick<Settings['academic'], 'keepTerms'>;
  signal?: AbortSignal;
}

export interface EngineUsage {
  input?: number;
  output?: number;
}

export interface EngineOutput {
  /** Eksik çeviriler undefined olabilir (servis bunları tek tek yeniden dener). */
  translations: (string | undefined)[];
  detected?: (string | undefined)[];
  usage?: EngineUsage;
}

export interface Engine {
  translate(input: EngineInput, config: EngineConfig): Promise<EngineOutput>;
  /** Serbest metin üretimi (terim çıkarma, açıklama). Yalnızca yapay zeka motorları. */
  complete?(system: string, user: string, config: EngineConfig, signal?: AbortSignal): Promise<{ text: string; usage?: EngineUsage }>;
}

export type EngineField = 'apiKey' | 'baseUrl' | 'model' | 'region' | 'temperature' | 'effort' | 'extraBody' | 'customHeaders';

export interface EngineMeta {
  type: EngineType;
  kind: 'mt' | 'llm';
  label: string;
  /** API anahtarı gerektirmez. */
  free?: boolean;
  /** API anahtarı zorunlu mu? */
  requiresKey: boolean;
  /** Zengin metin desteği: html = etiketli biçimi korur, none = düz metin. */
  richText: 'html' | 'none';
  fields: EngineField[];
  keyUrl?: string;
  docsUrl?: string;
  models?: string[];
}

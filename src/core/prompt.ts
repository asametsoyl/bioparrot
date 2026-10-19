import { englishName } from './languages';
import { formatSegments } from './segments';
import type { Expert, GlossaryEntry, Settings, TranslationContext } from './types';

export interface PromptInput {
  texts: string[];
  from: string;
  to: string;
  rich: boolean;
  expert: Expert;
  context?: TranslationContext;
  glossary: GlossaryEntry[];
  academic?: Pick<Settings['academic'], 'keepTerms'>;
}

export interface BuiltPrompt {
  /** İstekler arasında sabit kalan kısım (önbelleğe alınabilir). */
  stableSystem: string;
  /** Sayfaya/belgeye özel bağlam ve terimler. */
  contextBlock: string;
  /** stableSystem + contextBlock. */
  system: string;
  user: string;
}

function fill(template: string, vars: Record<string, string>): string {
  return template.replace(/\{\{\s*(\w+)\s*\}\}/g, (m, k: string) => vars[k] ?? m);
}

function clip(s: string | undefined, max: number): string {
  if (!s) return '';
  const t = s.replace(/\s+/g, ' ').trim();
  return t.length > max ? t.slice(0, max) + '…' : t;
}

export function buildPrompt(input: PromptInput): BuiltPrompt {
  const to = englishName(input.to);
  const from = input.from === 'auto' ? 'the source language' : englishName(input.from);
  const multi = input.texts.length > 1;
  const vars = {
    to,
    from,
    title: clip(input.context?.title, 300),
    summary: clip(input.context?.summary, 1500),
  };

  const parts: string[] = [fill(input.expert.systemPrompt, vars).trim()];

  if (input.expert.group === 'academic' && input.academic?.keepTerms) {
    parts.push(
      fill(
        `Terminology: when a key technical term has no well-established {{to}} equivalent, write the {{to}} rendering followed by the original term in parentheses the first time it appears in a segment, e.g. "<{{to}} term> (original term)". Do not do this for terms that already have a standard {{to}} translation.`,
        vars,
      ),
    );
  }

  const rules: string[] = [
    `Translate from ${from} into ${to}.`,
    'Output only the translation: no explanations, notes, alternatives, quotation marks or preambles.',
    `If a piece is already in ${to}, or is a name, number, code or URL that must not be translated, return it unchanged.`,
  ];
  if (multi) {
    rules.push(
      'The input consists of numbered segments in the form <seg id="n">…</seg>. Return every segment with the same id, in the same order and in the same form, containing only its translation. Never merge, split, skip or reorder segments, and never move text from one segment to another.',
    );
  }
  if (input.rich) {
    rules.push(
      'The text may contain HTML-like tags with an i attribute, e.g. <a i="0">…</a>, <b i="1">…</b>, <i i="2">…</i>, <span i="3">…</span>, <code class="notranslate" translate="no" i="4">…</code> and placeholders such as <img i="5"> (citation marks, symbols, formulas); <br> is a line break. Keep every tag and placeholder with exactly the same attributes, placed around or next to the words that correspond to the original ones, so formatting, links and citations stay correct. Copy the content of <code> tags verbatim. Keep HTML entities such as &amp; &lt; &gt; as they are.',
    );
  }
  parts.push('## Rules\n' + rules.map((r) => `- ${r}`).join('\n'));

  const stableSystem = parts.join('\n\n');

  const ctx: string[] = [];
  if (vars.title || vars.summary) {
    const lines = ['## Context (for understanding only; do not translate it or mention it)'];
    if (vars.title) lines.push(`Title: ${vars.title}`);
    if (vars.summary) lines.push(`Summary: ${vars.summary}`);
    ctx.push(lines.join('\n'));
  }
  const terms = [
    ...input.glossary.map((g) => ({ source: g.source.trim(), target: g.target.trim() })),
    ...(input.context?.terms ?? []),
  ];
  if (terms.length) {
    const seen = new Set<string>();
    const lines = terms
      .filter((t) => {
        const k = t.source.toLowerCase();
        if (!t.source || !t.target || seen.has(k)) return false;
        seen.add(k);
        return true;
      })
      .slice(0, 80)
      .map((t) => `- ${t.source} → ${t.target}`);
    ctx.push(`## Glossary (always use these translations for these terms)\n${lines.join('\n')}`);
  }
  const contextBlock = ctx.join('\n\n');

  const user = multi ? formatSegments(input.texts) : input.texts[0];

  return {
    stableSystem,
    contextBlock,
    system: contextBlock ? `${stableSystem}\n\n${contextBlock}` : stableSystem,
    user,
  };
}

/** Belge/sayfa için tutarlı terim listesi çıkarma istemi (yalnızca yapay zeka motorları). */
export function buildTermExtractionPrompt(text: string, to: string, expert: Expert): { system: string; user: string } {
  const target = englishName(to);
  return {
    system: `${fill(expert.systemPrompt, { to: target, from: 'the source language', title: '', summary: '' })}

Your task now is terminology extraction, not translation. From the text, select up to 40 domain-specific technical terms, multi-word expressions and recurring key concepts that a translator must render consistently. For each, give the best established ${target} translation used in ${target} scientific literature (keep gene/protein symbols, species names, formulas and acronyms unchanged and omit them). Respond with only a JSON array of objects of the form {"source": "...", "target": "..."}.`,
    user: text,
  };
}

export function parseTermList(output: string): { source: string; target: string }[] {
  const s = output.replace(/<think(?:ing)?>[\s\S]*?<\/think(?:ing)?>/gi, '');
  const start = s.indexOf('[');
  const end = s.lastIndexOf(']');
  if (start < 0 || end <= start) return [];
  try {
    const arr = JSON.parse(s.slice(start, end + 1));
    if (!Array.isArray(arr)) return [];
    return arr
      .filter((x) => x && typeof x.source === 'string' && typeof x.target === 'string')
      .map((x) => ({ source: x.source.trim(), target: x.target.trim() }))
      .filter((x) => x.source && x.target && x.source.toLowerCase() !== x.target.toLowerCase())
      .slice(0, 60);
  } catch {
    return [];
  }
}

/** Seçili metni açıklama istemi (seçim penceresi ve Zotero için). */
export function buildExplainPrompt(text: string, to: string, expert: Expert, context?: TranslationContext): { system: string; user: string } {
  const target = englishName(to);
  const ctx = context?.title ? `\nThe passage comes from: "${clip(context.title, 300)}".` : '';
  return {
    system: `${fill(expert.systemPrompt, { to: target, from: 'the source language', title: '', summary: '' })}

Your task now is to help a reader understand a passage. Write in ${target}. First give a faithful translation, then explain the meaning in plain terms: define the key technical terms, abbreviations and methods, and point out what the passage claims and how strong the claim is. Be concise (at most about 200 words), use short paragraphs or bullet points, and do not invent facts beyond the passage.${ctx}`,
    user: text,
  };
}

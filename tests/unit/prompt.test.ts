import { describe, expect, it } from 'vitest';
import { BUILTIN_EXPERTS } from '@/core/experts';
import { buildPrompt, parseTermList } from '@/core/prompt';

const expert = (id: string) => BUILTIN_EXPERTS.find((e) => e.id === id)!;

describe('prompt builder', () => {
  it('builds a multi-segment academic prompt with context, glossary and terms', () => {
    const p = buildPrompt({
      texts: ['Gene expression was measured.', 'Cells were <a i="0">imaged</a>.'],
      from: 'auto',
      to: 'tr',
      rich: true,
      expert: expert('biomedical'),
      context: { title: 'TP53 in organoids', summary: 'We show…', terms: [{ source: 'organoid', target: 'organoid' }] },
      glossary: [{ id: '1', source: 'gene expression', target: 'gen ifadesi' }],
      academic: { keepTerms: true },
    });
    expect(p.stableSystem).toContain('Turkish');
    expect(p.stableSystem).toContain('<seg id="n">');
    expect(p.stableSystem).toContain('i attribute');
    expect(p.stableSystem).toContain('original term in parentheses');
    expect(p.contextBlock).toContain('Title: TP53 in organoids');
    expect(p.contextBlock).toContain('gene expression → gen ifadesi');
    expect(p.user).toBe('<seg id="1">Gene expression was measured.</seg>\n<seg id="2">Cells were <a i="0">imaged</a>.</seg>');
  });

  it('sends a single text without segment wrappers', () => {
    const p = buildPrompt({ texts: ['Hello'], from: 'en', to: 'de', rich: false, expert: expert('general'), glossary: [] });
    expect(p.user).toBe('Hello');
    expect(p.system).not.toContain('<seg');
    expect(p.contextBlock).toBe('');
  });

  it('parses term lists from model output', () => {
    expect(parseTermList('Here: [{"source":"extracellular vesicle","target":"hücre dışı vezikül"},{"source":"TP53","target":"TP53"}]')).toEqual([
      { source: 'extracellular vesicle', target: 'hücre dışı vezikül' },
    ]);
    expect(parseTermList('no json')).toEqual([]);
  });
});

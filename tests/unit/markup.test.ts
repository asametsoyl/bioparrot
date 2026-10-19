import { describe, expect, it } from 'vitest';
import { applyGlossaryMarkup, escapeHtml, markupToPlain, plainToMarkup, relevantGlossary, stripGlossaryMarkup, tidyMarkup, tokenizeMarkup, unescapeHtml } from '@/core/markup';

describe('markup', () => {
  it('escapes and unescapes HTML entities', () => {
    expect(escapeHtml('a < b & c > d')).toBe('a &lt; b &amp; c &gt; d');
    expect(unescapeHtml('a &lt; b &amp; c &#39;x&#39; &#x2014; &nbsp;')).toBe("a < b & c 'x' —  ");
  });

  it('converts plain text with line breaks', () => {
    expect(plainToMarkup('x < 1\ny')).toBe('x &lt; 1<br>y');
    expect(markupToPlain('x &lt; 1<br>y')).toBe('x < 1\ny');
  });

  it('tokenizes tags with index attributes and void placeholders', () => {
    const tokens = tokenizeMarkup('Click <a i="0">here</a> to see<img i="1"> <code class="notranslate" translate="no" i="2">fetch()</code>');
    expect(tokens.filter((t) => t.type === 'open').map((t) => (t.type === 'open' ? [t.tag, t.index, t.selfClosing] : null))).toEqual([
      ['a', '0', false],
      ['img', '1', true],
      ['code', '2', false],
    ]);
    expect(markupToPlain('Click <a i="0">here</a>')).toBe('Click here');
  });

  it('applies glossary terms as protected spans and strips them afterwards', () => {
    const glossary = [{ id: '1', source: 'extracellular vesicle', target: 'hücre dışı vezikül' }];
    const marked = applyGlossaryMarkup('Each extracellular vesicle carries RNA.', glossary);
    expect(marked).toContain('<span class="notranslate" translate="no" i="g0">hücre dışı vezikül</span>');
    expect(stripGlossaryMarkup('Her <span class="notranslate" translate="no" i="g0">hücre dışı vezikül</span> RNA taşır.')).toBe('Her hücre dışı vezikül RNA taşır.');
  });

  it('matches glossary terms on word boundaries only', () => {
    const glossary = [{ id: '1', source: 'cell', target: 'hücre' }];
    expect(relevantGlossary(['The cellular response'], glossary)).toHaveLength(0);
    expect(relevantGlossary(['The cell response'], glossary)).toHaveLength(1);
  });

  it('tidies spaces that MT engines add around tags', () => {
    expect(tidyMarkup('( <i i="2">p</i> &lt; 0,05; bkz. <a i="3">Şekil 1</a> ).')).toBe('(<i i="2">p</i> &lt; 0,05; bkz. <a i="3">Şekil 1</a>).');
    expect(tidyMarkup('ölçüldü<img i="0"> .')).toBe('ölçüldü<img i="0">.');
  });
});

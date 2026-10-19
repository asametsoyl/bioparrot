import { beforeAll, describe, expect, it } from 'vitest';
import { Segmenter, shouldTranslate } from '@/content/segmenter';

// happy-dom tarayıcı varsayılan display değerlerini hesaplamaz; Chrome'un varsayılanlarını taklit et.
const INLINE = new Set(['A', 'B', 'I', 'EM', 'STRONG', 'SPAN', 'CODE', 'SUP', 'SUB', 'SMALL', 'MARK', 'ABBR', 'CITE', 'Q', 'TIME', 'LABEL', 'FONT', 'IMG', 'BR', 'KBD', 'VAR', 'SAMP']);
beforeAll(() => {
  window.getComputedStyle = ((el: Element) => ({
    display: el.tagName === 'LI' ? 'list-item' : INLINE.has(el.tagName) ? 'inline' : 'block',
  })) as unknown as typeof window.getComputedStyle;
});

const make = () => new Segmenter({ exclude: 'pre, .notranslate', skipCode: true, skipReferences: true });

describe('DOM segmenter', () => {
  it('serializes links, emphasis, code and citations into engine markup', () => {
    document.body.innerHTML = '<p id="p">See <a href="#f1">Figure 1</a> and <b>bold</b> text with <code>fetch()</code> data<sup>2</sup>.</p>';
    const seg = make();
    const { units } = seg.runsOf(document.getElementById('p')!);
    expect(units).toHaveLength(1);
    expect(units[0].whole).toBe(true);
    const ser = seg.serialize(units[0].nodes);
    expect(ser.markup).toBe('See <a i="0">Figure 1</a> and <b i="1">bold</b> text with <code class="notranslate" translate="no" i="2">fetch()</code> data<img i="3">.');
    expect(ser.elements.map((e) => e.tagName)).toEqual(['A', 'B', 'CODE', 'SUP']);
  });

  it('splits mixed content into separate runs around block children', () => {
    document.body.innerHTML = '<div id="d">Intro text <span>here</span><p>Nested paragraph</p>Trailing text</div>';
    const { units, blocks } = make().runsOf(document.getElementById('d')!);
    expect(units.map((u) => u.nodes.map((n) => n.textContent).join(''))).toEqual(['Intro text here', 'Trailing text']);
    expect(blocks.map((b) => (b as Element).tagName)).toEqual(['P']);
  });

  it('skips references after a references heading', () => {
    document.body.innerHTML = '<section id="s"><h2>Results</h2><p>Body</p><h2>References</h2><ol><li>Ref 1</li></ol></section>';
    const { blocks } = make().runsOf(document.getElementById('s')!);
    expect(blocks.map((b) => (b as Element).tagName)).toEqual(['H2', 'P', 'H2']);
  });

  it('filters strings that should not be translated', () => {
    expect(shouldTranslate('https://example.com/a', 'https://example.com/a')).toBe(false);
    expect(shouldTranslate('camelCaseIdentifier', 'camelCaseIdentifier')).toBe(false);
    expect(shouldTranslate('DNA', 'DNA')).toBe(false);
    expect(shouldTranslate('12.5 %', '12.5 %')).toBe(false);
    expect(shouldTranslate('Cells divide', 'Cells divide')).toBe(true);
  });
});

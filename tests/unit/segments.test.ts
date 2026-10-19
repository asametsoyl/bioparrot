import { describe, expect, it } from 'vitest';
import { cleanModelOutput, formatSegments, parseSegments } from '@/core/segments';

describe('LLM segments', () => {
  it('formats and parses numbered segments', () => {
    const input = formatSegments(['Hello', 'World']);
    expect(input).toBe('<seg id="1">Hello</seg>\n<seg id="2">World</seg>');
    expect(parseSegments('<seg id="1">Merhaba</seg>\n<seg id="2">Dünya</seg>', 2)).toEqual(['Merhaba', 'Dünya']);
  });

  it('returns undefined for missing segments', () => {
    expect(parseSegments('<seg id="2">Dünya</seg>', 3)).toEqual([undefined, 'Dünya', undefined]);
  });

  it('accepts untagged output for a single segment and strips wrappers', () => {
    expect(parseSegments('"Merhaba dünya"', 1)).toEqual(['Merhaba dünya']);
    expect(parseSegments('<think>reasoning…</think>\nMerhaba', 1)).toEqual(['Merhaba']);
    expect(cleanModelOutput('```\n<seg id="1">A</seg>\n```')).toBe('<seg id="1">A</seg>');
  });

  it('tolerates unclosed segments', () => {
    expect(parseSegments('<seg id="1">Bir <seg id="2">İki</seg>', 2)).toEqual(['Bir', 'İki']);
  });
});

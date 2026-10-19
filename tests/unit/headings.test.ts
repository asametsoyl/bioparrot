import { describe, expect, it } from 'vitest';
import { lookupSectionHeading, matchInitialCase } from '@/core/section-headings';

describe('section headings', () => {
  it('translates standard academic headings with numbering and case', () => {
    expect(lookupSectionHeading('Abstract', 'tr')).toBe('Özet');
    expect(lookupSectionHeading('2.1 Materials and Methods', 'tr')).toBe('2.1 Gereç ve Yöntemler');
    expect(lookupSectionHeading('IV. RESULTS', 'tr')).toBe('IV. BULGULAR');
    expect(lookupSectionHeading('Introduction', 'de')).toBeUndefined();
    expect(lookupSectionHeading('An introduction to cells', 'tr')).toBeUndefined();
  });

  it('matches the initial capital letter (Turkish dotted İ)', () => {
    expect(matchInitialCase('Introduction', 'giriş', 'tr')).toBe('Giriş');
    expect(matchInitialCase('Introduction', 'ilk adımlar', 'tr')).toBe('İlk adımlar');
    expect(matchInitialCase('lowercase', 'küçük', 'tr')).toBe('küçük');
  });
});

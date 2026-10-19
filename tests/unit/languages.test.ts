import { describe, expect, it } from 'vitest';
import { fromEngineCode, isSameLanguage, normalizeLang, toDeeplSource, toDeeplTarget, toGoogleCode, toMicrosoftCode } from '@/core/languages';

describe('languages', () => {
  it('normalizes codes from different sources', () => {
    expect(normalizeLang('en-US')).toBe('en');
    expect(normalizeLang('zh-Hant')).toBe('zh-TW');
    expect(normalizeLang('iw')).toBe('he');
    expect(normalizeLang('pt-BR')).toBe('pt');
    expect(normalizeLang('tr_TR')).toBe('tr');
  });

  it('compares languages but keeps Chinese scripts apart', () => {
    expect(isSameLanguage('en-GB', 'en')).toBe(true);
    expect(isSameLanguage('zh-CN', 'zh-TW')).toBe(false);
    expect(isSameLanguage('tr', 'az')).toBe(false);
  });

  it('maps to engine specific codes', () => {
    expect(toGoogleCode('he')).toBe('iw');
    expect(toMicrosoftCode('zh-CN')).toBe('zh-Hans');
    expect(toDeeplTarget('en')).toBe('EN-US');
    expect(toDeeplSource('pt')).toBe('PT');
    expect(toDeeplSource('auto')).toBeUndefined();
    expect(fromEngineCode('ku', 'microsoft')).toBe('ckb');
  });
});

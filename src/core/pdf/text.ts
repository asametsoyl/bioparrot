import { tokenizeMarkup } from '../markup';
import type { Block } from './layout';

export interface RichRun {
  text: string;
  /** Formül/kod gibi olduğu gibi korunan parça. */
  atom?: boolean;
  sup?: boolean;
  bold?: boolean;
  italic?: boolean;
}

/**
 * Çeviri metnini (zengin biçim) blok yer tutucularıyla birleştirip düz parçalara ayırır.
 * Hem HTML önizlemede hem PDF dışa aktarımda kullanılır.
 */
export function translationRuns(translation: string, rich: boolean, block: Block): RichRun[] {
  if (!rich) return [{ text: translation }];
  const runs: RichRun[] = [];
  const styles: string[] = [];
  let sink = 0;
  const used = new Set<number>();
  for (const token of tokenizeMarkup(translation)) {
    if (token.type === 'text') {
      if (!sink && token.text) runs.push({ text: token.text, bold: styles.includes('b'), italic: styles.includes('i') });
    } else if (token.type === 'br') {
      runs.push({ text: '\n' });
    } else if (token.type === 'open') {
      const idx = token.index && /^\d+$/.test(token.index) ? Number(token.index) : -1;
      const atom = idx >= 0 ? block.atoms[idx] : undefined;
      if ((token.tag === 'code' || token.tag === 'img') && atom) {
        if (!used.has(idx) || token.tag === 'code') runs.push({ text: atom.text, atom: true, sup: atom.sup });
        used.add(idx);
        if (!token.selfClosing) sink++;
      } else if (!token.selfClosing) {
        styles.push(token.tag);
      }
    } else if (token.type === 'close') {
      if (token.tag === 'code' && sink) sink--;
      else {
        const i = styles.lastIndexOf(token.tag);
        if (i >= 0) styles.splice(i, 1);
      }
    }
  }
  return runs;
}

export function runsToPlain(runs: RichRun[]): string {
  return runs
    .map((r) => (r.sup ? `${r.text}` : r.text))
    .join('')
    .replace(/[ \t]+/g, ' ')
    .trim();
}

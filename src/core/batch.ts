/** İndeksleri motorun parça ve karakter sınırlarına göre gruplar. */
export function makeBatches(indices: number[], texts: string[], maxItems: number, maxChars: number): number[][] {
  const batches: number[][] = [];
  let current: number[] = [];
  let chars = 0;
  const limitItems = Math.max(1, maxItems || 1);
  const limitChars = Math.max(200, maxChars || 2000);
  for (const i of indices) {
    const len = texts[i].length;
    if (current.length && (current.length >= limitItems || chars + len > limitChars)) {
      batches.push(current);
      current = [];
      chars = 0;
    }
    current.push(i);
    chars += len;
  }
  if (current.length) batches.push(current);
  return batches;
}

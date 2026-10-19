/** Yapay zeka modellerine birden çok parçayı tek istekte göndermek için biçim ve ayrıştırıcı. */

export function formatSegments(texts: string[]): string {
  return texts.map((t, i) => `<seg id="${i + 1}">${t}</seg>`).join('\n');
}

const SEG_RE = /<seg\s+id\s*=\s*["']?(\d+)["']?\s*>([\s\S]*?)(?:<\/seg>|(?=<seg\s+id\s*=)|$)/gi;

/** Modelin cevabından düşünme bloklarını, kod çitlerini ve gereksiz sarmalayıcıları temizler. */
export function cleanModelOutput(output: string): string {
  let s = output.replace(/<think(?:ing)?>[\s\S]*?<\/think(?:ing)?>/gi, '').trim();
  const fence = /^```[\w-]*\n([\s\S]*?)\n?```$/.exec(s);
  if (fence) s = fence[1].trim();
  return s;
}

/**
 * Parçalı çıktıyı ayrıştırır. Eksik parçalar undefined döner (çağıran bunları ayrıca çevirir).
 * Tek parçalı isteklerde etiketsiz çıktıyı da kabul eder.
 */
export function parseSegments(output: string, count: number): (string | undefined)[] {
  const s = cleanModelOutput(output);
  const result: (string | undefined)[] = new Array(count).fill(undefined);
  let found = 0;
  SEG_RE.lastIndex = 0;
  let m: RegExpExecArray | null;
  while ((m = SEG_RE.exec(s))) {
    const id = Number(m[1]);
    if (id >= 1 && id <= count && result[id - 1] === undefined) {
      result[id - 1] = m[2].trim();
      found++;
    }
  }
  if (found === 0 && count === 1) {
    result[0] = stripQuotes(s);
  }
  return result;
}

function stripQuotes(s: string): string {
  const t = s.trim();
  if (t.length > 1 && ((t.startsWith('"') && t.endsWith('"')) || (t.startsWith('“') && t.endsWith('”')))) {
    const inner = t.slice(1, -1);
    if (!inner.includes('"')) return inner;
  }
  return t;
}

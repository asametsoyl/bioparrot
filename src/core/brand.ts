/** BioParrot marka renkleri ve logo (papağan başı + iki dilli satırlar). */
export const BRAND = {
  green: '#16a34a',
  greenLight: '#22c55e',
  teal: '#0f766e',
  beak: '#fb923c',
  ink: '#0b3b33',
} as const;

/**
 * Papağan işareti (arka plansız, 128'lik görünüm kutusuna göre).
 * `eye` göz akı, `pupil` göz bebeği rengi.
 */
export function parrotMark(color = '#ffffff', eye: string = BRAND.ink, pupil: string = BRAND.ink): string {
  return `<path d="M57 27c-5-9-2-17 6-20-1 7 2 12 7 16-2-7 1-13 8-15-2 8 0 13 5 17-9-3-18-2-26 2z" fill="${color}" opacity=".92"/>
<path d="M40 104C29 74 35 37 66 29c17-4 30 6 32 19l-13 15c-4 13-8 27-12 41z" fill="${color}"/>
<path d="M93 43c12 0 21 8 19 22-1 6-6 10-11 10 3-6 2-12-4-14l-12 2z" fill="${BRAND.beak}"/>
<circle cx="73" cy="49" r="7.5" fill="${eye}"/><circle cx="74.5" cy="49" r="3.8" fill="${pupil}"/><circle cx="76" cy="47" r="1.7" fill="#fff"/>`;
}

/** Uygulama simgesi (yuvarlatılmış kare, yeşil degrade). */
export const LOGO_SVG = `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 128 128">
<defs><linearGradient id="bpg" x1="0" y1="0" x2="1" y2="1"><stop offset="0" stop-color="${BRAND.greenLight}"/><stop offset="1" stop-color="${BRAND.teal}"/></linearGradient></defs>
<rect width="128" height="128" rx="28" fill="url(#bpg)"/>
${parrotMark('#ffffff')}
<rect x="80" y="84" width="32" height="8" rx="4" fill="#fff" opacity=".95"/>
<rect x="80" y="98" width="22" height="8" rx="4" fill="#fff" opacity=".6"/>
</svg>`;

/** Arayüzde kullanılan küçük işaret (renkli papağan, şeffaf arka plan). */
export function markSvg(size = 24, color: string = BRAND.green): string {
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="24 2 96 106" width="${size}" height="${size}" aria-hidden="true">${parrotMark(color, '#ffffff', BRAND.ink)}</svg>`;
}

/*
 * md2docx.js — manuscript_draft.md / manuscript_draft_ko.md -> reading-copy .docx
 *
 *   node md2docx.js en   -> Manuscript_EN.docx
 *   node md2docx.js ko   -> Manuscript_KO.docx
 *
 * Not a submission format (that is oe_submission/main.tex) — this is the
 *版 people actually read and comment on: embedded figures, real sub/superscripts,
 * styled headings, hanging-indent references.
 */
const fs = require('fs');
const {
  Document, Packer, Paragraph, TextRun, HeadingLevel, AlignmentType, ImageRun,
  BorderStyle, convertInchesToTwip, PageOrientation, TabStopType, ExternalHyperlink,
} = require('docx');

const LANG = (process.argv[2] || 'en').toLowerCase();
const KO = LANG === 'ko';
const SRC = KO ? 'manuscript_draft_ko.md' : 'manuscript_draft.md';
const OUT = KO ? 'Manuscript_KO.docx' : 'Manuscript_EN.docx';

const BODY_FONT = KO ? 'Malgun Gothic' : 'Cambria';
const HEAD_FONT = KO ? 'Malgun Gothic' : 'Calibri';
const ACCENT = '1F3864';       // deep blue for headings
const CAPTION_GREY = '444444';

const FIGS = {
  '1': 'fig1_platform.png',
  '2': 'fig2_achievable_region.png',
  '3': 'fig3_selectivity_map.png',
  '4': 'fig4_recycling_routes.png',
  '5': 'fig5_families.png',
};

// ---------- PNG size from IHDR ----------
function pngSize(f) {
  const b = fs.readFileSync(f);
  return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) };
}

// ---------- LaTeX-ish math -> runs (real sub/superscripts) ----------
const SYM = {
  '\\times': '×', '\\approx': '≈', '\\ge': '≥', '\\geq': '≥',
  '\\le': '≤', '\\leq': '≤', '\\pm': '±', '\\sim': '~',
  '\\lesssim': '≲', '\\infty': '∞', '\\rightarrow': '→',
  '\\eta': 'η', '\\theta': 'θ', '\\lambda': 'λ', '\\mu': 'μ',
  '\\phi': 'φ', '\\gamma': 'γ', '\\Delta': 'Δ', '\\sigma': 'σ',
  '\\in': '∈', '\\cdot': '·', '\\mid': '|', '\\circ': '°',
  '\\max': 'max', '\\min': 'min', '\\log': 'log', '\\sin': 'sin', '\\cos': 'cos',
  '\\quad': '  ', '\\,': ' ', '\\;': ' ', '\\!': '', '\\ ': ' ',
  '\\{': '{', '\\}': '}', '\\%': '%', '\\&': '&',
};

// read a {...} group starting at s[i] === '{', return [content, indexAfter]
function takeGroup(s, i) {
  let depth = 0, j = i;
  for (; j < s.length; j++) {
    if (s[j] === '{') depth++;
    else if (s[j] === '}') { depth--; if (!depth) return [s.slice(i + 1, j), j + 1]; }
  }
  return [s.slice(i + 1), s.length];
}

// \frac{..}{..} -> (..)/(..), brace-matched so nested macros survive
function expandFrac(s) {
  let out = '', i = 0;
  while (i < s.length) {
    if (s.startsWith('\\frac', i)) {
      let k = i + 5;
      while (s[k] === ' ') k++;
      if (s[k] === '{') {
        const [num, a2] = takeGroup(s, k);
        let m = a2; while (s[m] === ' ') m++;
        if (s[m] === '{') {
          const [den, b2] = takeGroup(s, m);
          out += '(' + expandFrac(num) + ')/(' + expandFrac(den) + ')';
          i = b2; continue;
        }
      }
    }
    out += s[i]; i++;
  }
  return out;
}

function mathRuns(src, base = {}) {
  // strip wrappers that carry no visual meaning, innermost first
  let s = src;
  for (let n = 0; n < 4; n++) {
    s = s.replace(/\\(?:mathrm|mathbf|text|textrm|operatorname)\{([^{}]*)\}/g, '$1');
  }
  s = s.replace(/\\left|\\right/g, '');
  s = expandFrac(s);
  for (const [k, v] of Object.entries(SYM)) s = s.split(k).join(v);
  s = s.replace(/\\[a-zA-Z]+/g, '');            // drop any remaining macro
  const runs = [];
  let buf = '';
  const flush = (opts = {}) => {
    if (buf) runs.push(new TextRun({ text: buf, italics: true, font: BODY_FONT, ...base, ...opts }));
    buf = '';
  };
  for (let i = 0; i < s.length; i++) {
    const c = s[i];
    if ((c === '_' || c === '^') && i + 1 < s.length) {
      const sup = c === '^';
      let j = i + 1, txt;
      if (s[j] === '{') {
        let depth = 1; let k = j + 1;
        while (k < s.length && depth) { if (s[k] === '{') depth++; if (s[k] === '}') depth--; k++; }
        txt = s.slice(j + 1, k - 1); i = k - 1;
      } else { txt = s[j]; i = j; }
      flush();
      runs.push(new TextRun({
        text: txt.replace(/[{}]/g, ''), italics: true, font: BODY_FONT,
        subScript: !sup, superScript: sup, ...base,
      }));
    } else if (c === '{' || c === '}') {
      continue;
    } else buf += c;
  }
  flush();
  return runs;
}

// ---------- inline markdown -> runs ----------
function inlineRuns(text, base = {}) {
  const runs = [];
  // tokenizer over **bold** *italic* `code` $math$
  const re = /(\*\*[^*]+\*\*)|(\*[^*\n]+\*)|(`[^`]+`)|(\$[^$]+\$)/g;
  let last = 0, m;
  const plain = (t, extra = {}) => {
    if (!t) return;
    runs.push(new TextRun({ text: t, font: BODY_FONT, ...base, ...extra }));
  };
  while ((m = re.exec(text)) !== null) {
    plain(text.slice(last, m.index));
    const tok = m[0];
    if (tok.startsWith('**')) plain(tok.slice(2, -2), { bold: true });
    else if (tok.startsWith('`')) plain(tok.slice(1, -1), { font: 'Consolas', size: 18 });
    else if (tok.startsWith('$')) runs.push(...mathRuns(tok.slice(1, -1), base));
    else plain(tok.slice(1, -1), { italics: true });
    last = m.index + tok.length;
  }
  plain(text.slice(last));
  return runs;
}

const P = (text, opts = {}) => new Paragraph({
  children: inlineRuns(text, opts.runOpts || {}),
  spacing: { after: 140, line: 276 },
  alignment: AlignmentType.JUSTIFIED,
  ...opts.para,
});

// ---------- build ----------
const src = fs.readFileSync(SRC, 'utf8');
const title = src.match(/^# (.+)/m)[1];
const kids = [];

// title block
kids.push(new Paragraph({
  children: [new TextRun({ text: title, bold: true, size: 34, font: HEAD_FONT, color: ACCENT })],
  alignment: AlignmentType.CENTER, spacing: { after: 160 },
}));
kids.push(new Paragraph({
  children: [new TextRun({
    text: KO ? '저자 — [TBD]   |   소속 — [TBD]' : 'Authors — [TBD]   |   Affiliation — [TBD]',
    size: 20, color: '666666', font: HEAD_FONT,
  })],
  alignment: AlignmentType.CENTER, spacing: { after: 60 },
  border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: 'BFBFBF', space: 8 } },
}));
kids.push(new Paragraph({ text: '', spacing: { after: 160 } }));

const lines = src.split('\n');
let i = 0;
let inRefs = false;

// find where the body starts (after the title line)
while (i < lines.length && !lines[i].startsWith('## ')) i++;

const flushPara = (buf) => {
  const t = buf.join(' ').trim();
  if (!t) return;

  // figure caption block
  const fm = t.match(/^\*\*Fig\.\s*(\d)\s*\|/);
  if (fm && FIGS[fm[1]]) {
    const f = FIGS[fm[1]];
    const { w, h } = pngSize(f);
    const W = 600, H = Math.round((h / w) * W);
    kids.push(new Paragraph({
      children: [new ImageRun({ type: 'png', data: fs.readFileSync(f), transformation: { width: W, height: H } })],
      alignment: AlignmentType.CENTER, spacing: { before: 160, after: 80 },
    }));
    const cap = t.replace(/\s*\*\(`[^`]+`\)\*\s*$/, '');
    kids.push(new Paragraph({
      children: inlineRuns(cap, { size: 17, color: CAPTION_GREY }),
      alignment: AlignmentType.JUSTIFIED, spacing: { after: 200 },
      indent: { left: convertInchesToTwip(0.25), right: convertInchesToTwip(0.25) },
    }));
    return;
  }

  // reference entry: hanging indent, number in its own leading run
  if (inRefs && /^\d+\.\s/.test(t)) {
    const n = t.match(/^(\d+)\.\s+(.*)$/);
    kids.push(new Paragraph({
      children: [
        new TextRun({ text: n[1] + '.', bold: true, size: 18, font: BODY_FONT }),
        new TextRun({ text: '\t', size: 18, font: BODY_FONT }),
        ...inlineRuns(n[2], { size: 18 }),
      ],
      spacing: { after: 60 }, alignment: AlignmentType.JUSTIFIED,
      indent: { left: convertInchesToTwip(0.4), hanging: convertInchesToTwip(0.4) },
      tabStops: [{ type: TabStopType.LEFT, position: convertInchesToTwip(0.4) }],
    }));
    return;
  }

  kids.push(P(t));
};

let buf = [];
for (; i < lines.length; i++) {
  const L = lines[i];
  if (L.trim() === '---') { flushPara(buf); buf = []; continue; }

  if (L.startsWith('## ')) {
    flushPara(buf); buf = [];
    const name = L.slice(3).trim();
    inRefs = /References|참고문헌/.test(name);
    kids.push(new Paragraph({
      children: [new TextRun({ text: name, bold: true, size: 26, font: HEAD_FONT, color: ACCENT })],
      spacing: { before: 320, after: 140 },
      heading: HeadingLevel.HEADING_1,
      border: { bottom: { style: BorderStyle.SINGLE, size: 4, color: 'D0D7E5', space: 6 } },
    }));
    continue;
  }
  if (L.startsWith('### ')) {
    flushPara(buf); buf = [];
    kids.push(new Paragraph({
      children: [new TextRun({ text: L.slice(4).trim(), bold: true, size: 22, font: HEAD_FONT, color: '2E5496' })],
      spacing: { before: 240, after: 120 }, heading: HeadingLevel.HEADING_2,
    }));
    continue;
  }
  if (L.trim().startsWith('$$')) {           // display equation
    flushPara(buf); buf = [];
    const eq = [];
    i++;
    while (i < lines.length && !lines[i].trim().startsWith('$$')) { eq.push(lines[i]); i++; }
    kids.push(new Paragraph({
      children: mathRuns(eq.join(' ').trim(), { size: 22 }),
      alignment: AlignmentType.CENTER, spacing: { before: 120, after: 160 },
    }));
    continue;
  }
  if (L.trim() === '') { flushPara(buf); buf = []; continue; }
  buf.push(L.trim());
}
flushPara(buf);

const doc = new Document({
  creator: 'md2docx.js',
  title,
  styles: {
    default: {
      document: { run: { font: BODY_FONT, size: KO ? 20 : 21 } },
      heading1: { run: { font: HEAD_FONT, color: ACCENT } },
      heading2: { run: { font: HEAD_FONT, color: '2E5496' } },
    },
  },
  sections: [{
    properties: {
      page: {
        margin: {
          top: convertInchesToTwip(1), bottom: convertInchesToTwip(1),
          left: convertInchesToTwip(1), right: convertInchesToTwip(1),
        },
      },
    },
    children: kids,
  }],
});

Packer.toBuffer(doc).then((b) => {
  fs.writeFileSync(OUT, b);
  console.log(`${OUT} written (${kids.length} blocks, ${(b.length / 1024).toFixed(0)} KB)`);
});

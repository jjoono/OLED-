"""
md2optica.py -- manuscript_draft.md / supplementary.md -> Optica (Optics Express) LaTeX

Produces
  oe_submission/main.tex        opticajournal.cls manuscript (class file NOT
                                included -- see oe_submission/README.md)
  oe_submission/supplement.tex  standalone article-class supplement

Deterministic, re-runnable: edit the Markdown, re-run this, get updated .tex.
Conversion rules are conservative -- anything ambiguous is left as plain text
so it fails visibly in the PDF rather than silently.
"""
import re
from pathlib import Path

OUT = Path('oe_submission')
OUT.mkdir(exist_ok=True)

FIGFILES = {
    '1': 'figures/fig1_platform.pdf',
    '2': 'figures/fig2_achievable_region.pdf',
    '3': 'figures/fig3_selectivity_map.pdf',
    '4': 'figures/fig4_recycling_routes.pdf',
    '5': 'figures/fig5_families.pdf',
    'S1': 'figures/figS1_patch_dependence.pdf',
    'S2': 'figures/figS2_warmstart_control.pdf',
    'S3': 'figures/figS3_cost_calibration.pdf',
    'S4': 'figures/figS4_convergence.pdf',
}


# ---------------- character / inline conversions ----------------
def esc_text(t):
    """Escape/convert a Markdown text block into LaTeX body text."""
    t = t.replace('\\%', '\x00PCT\x00')          # already-escaped percents
    t = t.replace('%', '\\%').replace('\x00PCT\x00', '\\%')
    t = t.replace('&', '\\&')
    # unicode -> mode-safe macros (\ensuremath works in text and math alike)
    t = (t.replace('°', '\\ensuremath{^\\circ}')
          .replace('×', '\\ensuremath{\\times}')
          .replace('±', '\\ensuremath{\\pm}')
          .replace('μ', '\\ensuremath{\\mu}')
          .replace('−', '-'))
    t = re.sub(r'(?<![\\\w])~', '\\\\ensuremath{\\\\sim}', t)   # bare tilde
    # `code` -> \texttt{} with escaped underscores
    t = re.sub(r'`([^`]+)`', lambda m: '\\texttt{' + m.group(1).replace('_', '\\_') + '}', t)
    # bold / italics (bold first)
    t = re.sub(r'\*\*([^*]+)\*\*', r'\\textbf{\1}', t)
    t = re.sub(r'(?<!\*)\*([^*\n]+)\*(?!\*)', r'\\emph{\1}', t)
    # [n], [n,m], [n-m] citations -> \cite (avoid [TBD] etc.)
    def cite(m):
        keys = []
        for part in re.split(r'[,\s]+', m.group(1).strip()):
            if '–' in part or '-' in part:
                a, b = re.split(r'[–-]', part)
                keys += [f'r{k}' for k in range(int(a), int(b) + 1)]
            elif part:
                keys.append(f'r{int(part)}')
        return '\\cite{' + ','.join(keys) + '}'
    t = re.sub(r'\[(\d+(?:\s*[,–-]\s*\d+)*)\]', cite, t)
    # figure references (not supplementary S-figures, not other papers' figures)
    t = t.replace('Fig. 3 and Table 2 of ref.', '\x01EXTFIG\x01')   # ref-17's own figure
    t = re.sub(r'(?:Figure|Fig\.)\s?(\d)', r'Fig.~\\ref{fig:\1}', t)
    t = t.replace('\x01EXTFIG\x01', 'Fig.~3 and Table~2 of ref.')
    # daggers / arrows (supplement footnotes)
    t = (t.replace('†', '\\ensuremath{\\dagger}')
          .replace('‡', '\\ensuremath{\\ddagger}')
          .replace('→', '\\ensuremath{\\rightarrow}'))
    return t


def make_figure(num, caption):
    wide = '*' if num in ('1', '2', '5') else ''
    width = '\\textwidth' if wide else '\\linewidth'
    return (f'\\begin{{figure{wide}}}[t]\n\\centering\n'
            f'\\includegraphics[width={width}]{{{FIGFILES[num]}}}\n'
            f'\\caption{{{caption}}}\n\\label{{fig:{num}}}\n\\end{{figure{wide}}}\n')


# ---------------- main manuscript ----------------
src = Path('manuscript_draft.md').read_text(encoding='utf-8')

title = re.match(r'# (.+)', src).group(1)
abstract = src.split('## Abstract', 1)[1].split('---', 1)[0].strip()
body = src.split('## 1. Introduction', 1)[1]
body = '## 1. Introduction' + body
refs_md = body.split('## References', 1)[1]
body = body.split('---\n\n## Funding', 1)[0]
backm = src.split('## Funding', 1)[1].split('## References', 1)[0]

out = []
out.append(r'''% !TEX program = pdflatex
% ============================================================
% Optics Express manuscript -- generated from manuscript_draft.md by md2optica.py
% Requires opticajournal.cls (not distributable; see README.md)
% ============================================================
\documentclass[9pt,twocolumn,twoside]{opticajournal}
\journal{opticajournal}   % Optics Express uses the opticajournal layout
\usepackage{lineno}
\linenumbers

\begin{document}

\title{''' + esc_text(title) + r'''}

% ---- authors: fill in ----
\author{Author One\authormark{1}, Author Two\authormark{1} and Author Three\authormark{1,*}}
\address{\authormark{1}[Affiliation, Address]\\}
\email{\authormark{*}[corresponding@email]}
''')

out.append('\\begin{abstract*}\n' + esc_text(abstract) + '\n\\end{abstract*}\n')

# ---- body: walk line groups ----
paras = body.split('\n\n')
for p in paras:
    p = p.strip()
    if not p or p == '---':
        continue
    if p.startswith('## '):                              # \section
        name = re.sub(r'^## \d+\.\s*', '', p[3:])
        out.append('\\section{' + esc_text(name) + '}')
    elif p.startswith('### '):                           # \subsection
        name = re.sub(r'^### \d+\.\d+\s*', '', p[4:])
        out.append('\\subsection{' + esc_text(name) + '}')
    elif p.startswith('$$'):                             # display equation
        eq = p.strip('$ \n')
        out.append('\\begin{equation}\n' + eq + '\n\\end{equation}')
    elif re.match(r'\*\*Fig\. (\d) \|', p):              # figure caption block
        num = re.match(r'\*\*Fig\. (\d) \|', p).group(1)
        cap = re.sub(r'^\*\*Fig\. \d \|\s*', '**', p)
        cap = re.sub(r'\s*\*\(`[^`]+`\)\*\s*$', '', cap)     # strip file marker
        out.append(make_figure(num, esc_text(cap)))
    else:
        out.append(esc_text(p))

# ---- back matter ----
bm = {}
cur = 'Funding'; bm[cur] = []
for line in ('## Funding\n' + backm).split('\n'):
    m = re.match(r'## (.+)', line)
    if m:
        cur = m.group(1); bm[cur] = []
    elif line.strip() and line.strip() != '---':
        bm[cur].append(line)
out.append('\\begin{backmatter}')
KOREAN_TBD = {'[TBD — 과제/기관 정보 입력 필요]': '[Funding agencies and grant numbers]',
              '[TBD]': '[Acknowledgments]'}
for sec in ['Funding', 'Acknowledgments', 'Disclosures', 'Data availability']:
    raw = '\n'.join(bm.get(sec, [])).strip()
    for k, v in KOREAN_TBD.items():
        raw = raw.replace(k, v)
    body_bm = esc_text(raw)
    out.append(f'\\bmsection{{{sec}}}\n{body_bm}\n')
out.append('\\bmsection{Supplemental document}\nSee Supplement 1 for supporting content.\n')
out.append('\\end{backmatter}')

# ---- references: Markdown "Last, F. M.; ..." -> Optica "F. M. Last, ..." ----
def optica_authors(raw):
    """'Brütting, W.; Frischeisen, J.; et al.' -> 'W. Brütting, J. Frischeisen, et al.'"""
    parts = [a.strip().rstrip('.') + ('.' if a.strip().endswith('.') else '')
             for a in raw.split(';') if a.strip()]
    names, etal = [], False
    for a in parts:
        a = a.strip()
        if a.lower().startswith('et al'):
            etal = True; continue
        if ',' in a:                       # "Last(, compound), Initials"
            last, inits = a.rsplit(',', 1)
            names.append(inits.strip() + ' ' + last.strip())
        else:
            names.append(a)
    if etal:
        return names[0] + ' et al.' if names else 'et al.'
    if len(names) == 1:
        return names[0]
    if len(names) == 2:
        return names[0] + ' and ' + names[1]
    return ', '.join(names[:-1]) + ', and ' + names[-1]

def optica_item(it):
    it = ' '.join(it.split())
    m = re.match(r'(.+?)\*\*(.+?)\*\*\s*(.*)', it)
    if not m:
        return it                                          # 형식 불명 -> 원문 유지
    authors, title, tail = m.group(1).strip(), m.group(2).strip(), m.group(3).strip()
    title = title.rstrip('.').rstrip(',')
    # tail: *Journal* **year**, *vol*, pages/artno .  [ URL/DOI ]
    jm = re.match(r'\*([^*]+)\*\s+\*\*(\d{4})\*\*,\s*\*(\d+)\*,\s*([\w\u2013\u2014-]+)\.?\s*(.*)', tail)
    if jm:
        journal, year, vol, pages, rest = jm.groups()
        ref = (optica_authors(authors) + ', ``' + title + ',\'\' '
               + journal.strip() + ' \\textbf{' + vol + '}, ' + pages + ' (' + year + ').')
    else:
        # 저널 패턴이 안 잡히면 최소 변환 (저자 순서만 교정)
        rest = re.sub(r'\*\*(\d{4})\*\*', r'(\1)', tail)
        rest = re.sub(r'\*([^*]+)\*', r'\1', rest)
        ref = optica_authors(authors) + ', ``' + title + ',\'\' ' + rest
    rest_url = re.search(r'https?://\S+', rest if jm else '')
    if rest_url:
        ref += ' \\url{' + rest_url.group().rstrip('.') + '}'
    ref = ref.replace('&', '\\&').replace('%', '\\%')
    ref = re.sub(r'https?://\S+(?<!\})', lambda m: m.group() if '\\url' in ref[:m.start()][-30:] else m.group(), ref)
    return ref

items = re.findall(r'^\d+\.\s+(.+?)(?=\n\n\d+\.|\n*$)', refs_md.strip(), re.S | re.M)
out.append('\\begin{thebibliography}{%d}' % len(items))
for i, it in enumerate(items, 1):
    out.append(f'\\bibitem{{r{i}}} {optica_item(it)}\n')
out.append('\\end{thebibliography}\n\n\\end{document}')

(Path(OUT) / 'main.tex').write_text('\n\n'.join(out), encoding='utf-8')
print(f'main.tex written ({len(items)} references)')


# ---------------- supplement ----------------
ssrc = Path('supplementary.md').read_text(encoding='utf-8')

def md_table(block):
    rows = [r.strip().strip('|').split('|') for r in block.split('\n') if r.strip().startswith('|')]
    rows = [[c.strip() for c in r] for r in rows]
    rows = [r for r in rows if not all(set(c) <= set('-: ') for c in r)]  # drop rule row
    ncol = max(len(r) for r in rows)
    colspec = 'l' * ncol
    lines = ['\\begin{table}[htbp]', '\\centering', '\\footnotesize',
             '\\begin{tabular}{' + colspec + '}', '\\hline']
    for i, r in enumerate(rows):
        r = r + [''] * (ncol - len(r))
        lines.append(' & '.join(esc_text(c) for c in r) + ' \\\\')
        if i == 0:
            lines.append('\\hline')
    lines += ['\\hline', '\\end{tabular}', '\\end{table}']
    return '\n'.join(lines)

sout = [r'''% Supplement 1 -- generated from supplementary.md by md2optica.py
\documentclass[10pt]{article}
\usepackage[margin=2.5cm]{geometry}
\usepackage{graphicx,amsmath,url}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\title{Supplement 1:\\''' + esc_text(title) + r'''}
\author{}
\date{}
\begin{document}
\maketitle
''']

blocks = ssrc.split('\n\n')
i = 0
while i < len(blocks):
    b = blocks[i].strip()
    if not b or b == '---':
        i += 1; continue
    if b.startswith('# '):
        i += 1; continue                                  # doc title line
    if b.startswith('## '):
        name = b[3:]
        sout.append('\\section*{' + esc_text(name) + '}')
        # attach the matching supplementary figure right after its table section
        m = re.search(r'\(Fig\. (S\d)\)', name)
        if m:
            sn = m.group(1)
            sout.append(f'\\begin{{figure}}[htbp]\n\\centering\n'
                        f'\\includegraphics[width=\\linewidth]{{{FIGFILES[sn]}}}\n'
                        f'\\caption{{Fig. {sn}. Graphical form of this section\'s data.}}\n'
                        f'\\end{{figure}}')
    elif b.startswith('|'):
        sout.append(md_table(b))
    elif b.startswith('```'):
        sout.append('\\begin{verbatim}\n' + b.strip('`\n') + '\n\\end{verbatim}')
    elif re.match(r'\d+\.\s', b):                          # numbered list
        items_ = re.split(r'\n(?=\d+\.\s)', b)
        sout.append('\\begin{enumerate}')
        for it in items_:
            sout.append('\\item ' + esc_text(re.sub(r'^\d+\.\s*', '', ' '.join(it.split()))))
        sout.append('\\end{enumerate}')
    else:
        sout.append(esc_text(b))
    i += 1

sout.append('\\end{document}')
stext = '\n\n'.join(sout)
# 보충 문서에는 본문 그림 label 이 없으므로 \ref 를 일반 텍스트로 되돌린다
stext = re.sub(r'Fig\.~\\ref\{fig:(\d)\}', r'Fig. \1', stext)
(Path(OUT) / 'supplement.tex').write_text(stext, encoding='utf-8')
print('supplement.tex written')

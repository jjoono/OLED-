"""
md2optica.py -- manuscript_draft.md / supplementary.md -> Optica (Optics Express) LaTeX

Produces
  oe_submission/main.tex        optica-article.cls manuscript (Optica universal template; class
                                file NOT included -- see oe_submission/README.md)
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
    '1': 'figures/fig1_mixing_concept.pdf',
    '2': 'figures/fig2_complete_mixing.pdf',
    '3': 'figures/fig3_reduced_model.pdf',
    '4': 'figures/fig4_structure_map.pdf',
    '5': 'figures/fig5_experiment.pdf',
    'S1': 'figures/figS11_device_inputs.pdf',
    'S2': 'figures/figS_model_comparison.pdf',
    'S3': 'figures/figS2_composition_window.pdf',
    'S4': 'figures/fig5_bound_window.pdf',
    'S5': 'figures/fig1_platform.pdf',
    'S6': 'figures/figS5_optimizer_convergence.pdf',
    'S7': 'figures/figS1_restart_control.pdf',
    'S8': 'figures/fig3_efficiency_composition.pdf',
    'S9': 'figures/fig6_hemisphere_benchmark.pdf',
    'S10': 'figures/figS_finite_size.pdf',
}
# source of every figure PDF (copied into oe_submission/figures/ on each run)
import shutil
for _f in FIGFILES.values():
    _n = Path(_f).name
    _src = Path('figures_oe') / _n if (Path('figures_oe') / _n).exists() else Path(_n)
    (OUT / 'figures').mkdir(exist_ok=True)
    if _src.exists():
        shutil.copy(_src, OUT / 'figures' / _n)
    else:
        print('missing figure', _src)


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
          .replace('λ', '\\ensuremath{\\lambda}')
          .replace('Δ', '\\ensuremath{\\Delta}')
          .replace('≥', '\\ensuremath{\\geq}')
          .replace('≤', '\\ensuremath{\\leq}')
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
    t = re.sub(r'\bFigure\s?(\d)', r'Figure~\\ref{fig:\1}', t)
    t = re.sub(r'\bFig\.\s?(\d)', r'Fig.~\\ref{fig:\1}', t)
    t = re.sub(r'(?<![\w|])Table (\d)(?![\d|])', r'Table~\\ref{tab:\1}', t)
    t = re.sub(r'\b(Eq\.|Equation)\s?\((\d+)\)', r'\1~(\\ref{eq:\2})', t)
    t = re.sub(r'"([^"]+)"', r"``\1''", t)                     # straight quotes -> TeX quotes
    # daggers / arrows (supplement footnotes)
    t = (t.replace('†', '\\ensuremath{\\dagger}')
          .replace('‡', '\\ensuremath{\\ddagger}')
          .replace('→', '\\ensuremath{\\rightarrow}'))
    return t


def make_figure(num, caption):
    return (f'\\begin{{figure}}[htbp]\n\\centering\n'
            f'\\includegraphics[width=\\linewidth]{{{FIGFILES[num]}}}\n'
            f'\\caption{{{caption}}}\n\\label{{fig:{num}}}\n\\end{{figure}}\n')


def md_table(block, caption=None, label=None, star=False):
    rows = [r.strip().strip('|').split('|') for r in block.split('\n') if r.strip().startswith('|')]
    rows = [[c.strip() for c in r] for r in rows]
    rows = [r for r in rows if not all(set(c) <= set('-: ') for c in r)]  # drop rule row
    ncol = max(len(r) for r in rows)
    colspec = 'l' * ncol
    env = 'table*' if star else 'table'
    lines = ['\\begin{' + env + '}[htbp]', '\\centering', '\\footnotesize']
    if caption:
        lines.append('\\caption{' + esc_text(caption) + '}')
        if label:
            lines.append('\\label{' + label + '}')
    if star:
        colspec = ''.join(['p{%.3f\\linewidth}' % (0.92 / ncol - 0.02)] * ncol)
    wrap = (not star) and ncol >= 4
    if not star and not wrap and max(len(c) for r in rows for c in r) > 45:
        w = 0.80 / (ncol - 1) if ncol > 1 else 0.9
        colspec = 'p{0.16\\linewidth}' + ''.join('p{%.2f\\linewidth}' % w for _ in range(ncol - 1))
    lines += [('\\resizebox{\\linewidth}{!}{' if wrap else '') + '\\begin{tabular}{' + colspec + '}', '\\hline']
    for i, r in enumerate(rows):
        r = r + [''] * (ncol - len(r))
        lines.append(' & '.join(esc_text(c) for c in r) + ' \\\\')
        if i == 0:
            lines.append('\\hline')
    lines += ['\\hline', '\\end{tabular}' + ('}' if wrap else ''), '\\end{' + env + '}']
    return '\n'.join(lines)



def md_list(block):
    items_ = re.split(r'\n(?=- )', block.strip())
    out_ = ['\\begin{itemize}']
    for it in items_:
        out_.append('\\item ' + esc_text(' '.join(it[2:].split())))
    out_.append('\\end{itemize}')
    return '\n'.join(out_)

# ---------------- main manuscript ----------------
src = Path('manuscript_draft.md').read_text(encoding='utf-8')

title = re.match(r'# (.+)', src).group(1)
abstract = src.split('## Abstract', 1)[1].split('---', 1)[0].strip()
body = '## 1. Introduction' + src.split('## 1. Introduction', 1)[1]
refs_md = body.split('## References', 1)[1]
body = body.split('## Funding', 1)[0]
backm = src.split('## Funding', 1)[1].split('## References', 1)[0]

CAP = r'\*\*(Fig\.|Table) (\d+)(?: \||\.\*\*)\s*'          # '**Fig. 1.** text' (also the older '**Fig. 1 | text**')


def caption_text(p):
    """Caption without its label and without the trailing file marker."""
    p = re.sub(r'\s*\*\(`[^`]+`\)\*\s*$', '', p.strip())
    old = re.match(r'\*\*(?:Fig\.|Table) \d+ \|\s*(.+?)\*\*\s*(.*)$', p, re.S)
    if old:
        return (old.group(1) + ' ' + old.group(2)).strip()
    return re.sub('^' + CAP, '', p).strip()


HEAD = r"""% !TEX program = pdflatex
% ============================================================
% Optics Express manuscript -- generated from manuscript_draft.md by md2optica.py
% Requires the Optica universal template (optica-article.cls, opticajnl.bst; see README.md)
% ============================================================
\documentclass{optica-article}
\journal{opticajournal}
\articletype{Research Article}

\begin{document}

\title{TITLE}

% ---- authors: fill in ----
\author{Author One,\authormark{1} Author Two,\authormark{1} and Author Three\authormark{1,*}}
\address{\authormark{1}[Affiliation, Address]}
\email{\authormark{*}[corresponding@email]}
"""
out = [HEAD.replace('TITLE', esc_text(title))]
out.append('\\begin{abstract*}\n' + esc_text(abstract) + '\n\\end{abstract*}\n')

# ---- body: walk paragraph blocks ----
pending_table = []
for p in body.split('\n\n'):
    p = p.strip()
    if not p or p == '---':
        continue
    if p.startswith('## '):                              # \section
        out.append('\\section{' + esc_text(re.sub(r'^\d+\.\s*', '', p[3:])) + '}')
    elif p.startswith('### '):                           # \subsection
        out.append('\\subsection{' + esc_text(re.sub(r'^\d+\.\d+\s*', '', p[4:])) + '}')
    elif p.startswith('$$'):                             # display equation, numbered by \tag{n} in the source
        eq = p.strip('$ \n'); tag = re.search(r'\\tag\{(\d+)\}', eq)
        eq = re.sub(r'\s*\\tag\{\d+\}', '', eq).strip()
        out.append('\\begin{equation}\n' + (f'\\label{{eq:{tag.group(1)}}}\n' if tag else '') + eq + '\n\\end{equation}')
    elif re.match(r'\*\*Table \d+', p):                  # table caption; the table follows
        pending_table[:] = [re.match(r'\*\*Table (\d+)', p).group(1), caption_text(p)]
    elif p.startswith('|'):                              # markdown table
        if pending_table:
            out.append(md_table(p, pending_table[1], 'tab:' + pending_table[0], star=True).replace('table*', 'table'))
            pending_table.clear()
        else:
            out.append(md_table(p))
    elif p.startswith('- '):                             # bullet list
        out.append(md_list(p))
    elif re.match(r'\*\*Fig\. \d+', p):                  # figure caption block
        out.append(make_figure(re.match(r'\*\*Fig\. (\d+)', p).group(1), esc_text(caption_text(p))))
    else:
        out.append(esc_text(p))

# ---- back matter ----
bm = {}; cur = 'Funding'; bm[cur] = []
for line in ('## Funding\n' + backm).split('\n'):
    m = re.match(r'## (.+)', line)
    if m:
        cur = m.group(1); bm[cur] = []
    elif line.strip() and line.strip() != '---':
        bm[cur].append(line)
out.append('\\begin{backmatter}')
TBD = {'Funding': '[Funding agencies and grant numbers]', 'Acknowledgment': '[Acknowledgment]'}
for sec in ['Funding', 'Acknowledgment', 'Disclosures', 'Data availability', 'Supplemental document']:
    raw = '\n'.join(bm.get(sec, bm.get(sec + 's', []))).strip()
    if raw == '[TBD]':
        raw = TBD.get(sec, raw)
    out.append(f'\\bmsection{{{sec}}}\n{esc_text(raw)}\n')
out.append('\\end{backmatter}')


# ---- references: already in Optica style in the Markdown ----
def bib_item(it):
    it = ' '.join(it.split())
    it = it.replace('&', '\\&').replace('%', '\\%')
    it = re.sub(r'"([^"]+)"', r"``\1''", it)
    it = re.sub(r'\*\*([^*]+)\*\*', r'\\textbf{\1}', it)
    it = re.sub(r'\*([^*]+)\*', r'\\emph{\1}', it)
    return it


items = re.findall(r'^\d+\.\s+(.+?)(?=\n\n\d+\.|\n*$)', refs_md.strip(), re.S | re.M)
out.append('\\begin{thebibliography}{%d}' % len(items))
for i, it in enumerate(items, 1):
    out.append(f'\\bibitem{{r{i}}} {bib_item(it)}\n')
out.append('\\end{thebibliography}\n\n\\end{document}')

(Path(OUT) / 'main.tex').write_text('\n\n'.join(out), encoding='utf-8')
print(f'main.tex written ({len(items)} references)')


# ---------------- supplement ----------------
ssrc = Path('supplementary.md').read_text(encoding='utf-8')

sout = [r'''% Supplement 1 -- generated from supplementary.md by md2optica.py
\documentclass[10pt]{article}
\usepackage[margin=2.5cm]{geometry}
\usepackage{graphicx,amsmath,url,caption}
\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\title{Supplement 1:\\''' + esc_text(title) + r'''}
\author{}
\date{}
\begin{document}
\maketitle
\sloppy
''']

def scap(b):
    """'**Fig. S1.** text' (or the older '**Fig. S1 | text**') -> '**Fig. S1.** text' without the file marker."""
    b = re.sub(r'\s*\*\(`[^`]+`\)\*\s*$', '', ' '.join(b.split()))
    m = re.match(r'\*\*((?:Fig\.|Table) S?\d+) \|\s*(.+?)\*\*\s*(.*)$', b)
    return f'**{m.group(1)}.** {m.group(2)} {m.group(3)}'.strip() if m else b


blocks = ssrc.split('\n\n')
pend = None                                               # caption waiting for its table
for b in blocks:
    b = b.strip()
    if not b or b == '---' or b.startswith('# '):
        continue
    if b.startswith('## '):
        sout.append('\\section*{' + esc_text(b[3:]) + '}')
    elif re.match(r'\*\*Table S?\d+', b):
        pend = scap(b)
    elif re.match(r'\*\*Fig\. S\d+', b):
        sn = re.match(r'\*\*Fig\. (S\d+)', b).group(1)
        sout.append(f'\\begin{{figure}}[htbp]\n\\centering\n'
                    f'\\includegraphics[width=0.85\\linewidth]{{{FIGFILES[sn]}}}\n'
                    f'\\caption*{{' + esc_text(scap(b)) + '}\n\\end{figure}')
    elif b.startswith('|'):
        tex = md_table(b)
        if pend:
            tex = tex.replace('\\footnotesize', '\\footnotesize\n\\caption*{' + esc_text(pend) + '}', 1)
            pend = None
        sout.append(tex)
    elif b.startswith('$$'):
        sout.append('\\begin{equation*}\n' + re.sub(r'\s*\\tag\{[^}]*\}', '', b.strip('$ \n')) + '\n\\end{equation*}')
    elif b.startswith('- '):
        sout.append(md_list(b))
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

sout.append('\\end{document}')
stext = '\n\n'.join(sout)
# 보충 문서에는 본문 그림 label 이 없으므로 \ref 를 일반 텍스트로 되돌린다
stext = re.sub(r'Fig\.~\\ref\{fig:(\d)\}', r'Fig. \1', stext)
stext = re.sub(r'Table~\\ref\{tab:(\d)\}', r'Table \1', stext)
(Path(OUT) / 'supplement.tex').write_text(stext, encoding='utf-8')
print('supplement.tex written')

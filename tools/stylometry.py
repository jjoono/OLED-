"""Stylometric comparison of a manuscript PDF with reference journal articles (PDF).

    python3 tools/stylometry.py REF_DIR MANUSCRIPT.pdf [more.pdf ...]

REF_DIR holds the reference PDFs (kept outside the repository).  Every PDF goes through the same pipeline:
pdftotext -> body between "1. Introduction" and the back matter -> prose blocks only (captions, headings, tables,
figure lettering and display equations are separated) -> counts per 1000 prose words.  Prints the reference mean
and standard deviation, each manuscript value and its z-score."""
import glob, os, re, subprocess, sys
import numpy as np


def text_of(pdf):
    t = subprocess.run(['pdftotext', pdf, '-'], capture_output=True, text=True).stdout
    t = t.replace('\f', '\n\n')
    t = re.sub(r'^\s*Research Article\s*$', '', t, flags=re.M)
    t = re.sub(r'^\s*Vol\. \d+, No\. \d+ / .*$', '', t, flags=re.M)
    t = re.sub(r'^#\d+.*$|^Journal © .*$|^https://doi\.org/\S+\s*$|^Received .*published.*$', '', t, flags=re.M)
    return t


def split(t):
    ab = re.search(r'Abstract:?\s*(.*?)(?:©|\n\s*1\.\s+Introduction)', t, re.S)
    i = re.search(r'^\s*1\.\s+Introduction', t, re.M)
    body = t[i.end():] if i else t
    refs = ''
    j = [m.start() for m in re.finditer(r'^\s*References\s*$', body, re.M)]
    if j: refs = body[j[-1]:]; body = body[:j[-1]]
    k = re.search(r'^\s*Funding\b\.?', body, re.M)
    if k: body = body[:k.start()]
    return (ab.group(1) if ab else ''), body, refs


CALL = r'(?:Figs?\.|Figures?|Eqs?\.|Equations?|Tables?|Sections?|Secs?\.|Refs?\.|Notes?|Supplement|Dataset|Code|Visualization)'


def analyse(pdf):
    t = text_of(pdf); ab, body, refs = split(t)
    blocks = [b for b in re.split(r'\n\s*\n', body) if b.strip()]
    prose, caps, heads = [], [], []
    for b in blocks:
        one = ' '.join(b.split())
        lines = [l for l in b.split('\n') if l.strip()]
        if re.match(r'(Fig\.|Table) \d+\.', one): caps.append(one); continue
        if re.match(r'\d+(\.\d+)?\.?\s+[A-Z]', one) and len(one.split()) < 16: heads.append(one); continue
        alpha = sum(c.isalpha() or c == ' ' for c in one) / max(1, len(one))
        if len(one.split()) >= 12 and alpha > 0.72 and np.mean([len(l) for l in lines]) > 45: prose.append(b)
    # join lines; a hyphen at a line end is a soft hyphen of the typesetter
    txt = '\n'.join(re.sub(r'(\w)-\n(\w)', r'\1\2', b).replace('\n', ' ') for b in prose)
    words = re.findall(r"[A-Za-z][A-Za-z'’]*", txt); W = len(words); k = 1000.0 / W
    cites = re.findall(r'\[\d+(?:[,–\-]\s*\d+)*\]', txt)
    t2 = re.sub(r'\[\d+(?:[,–\-]\s*\d+)*\]', ' ', txt)
    calls = re.findall(CALL + r'\s*\(?S?\d+\)?(?:\([a-z]\))?(?:\s*(?:–|-|and|,)\s*S?\d*\(?[a-z0-9]*\)?)*', t2)
    t3 = re.sub(CALL + r'\s*\(?S?\d+\)?(?:\([a-z]\))?(?:\s*(?:–|-|and|,)\s*S?\d*\(?[a-z0-9]*\)?)*', ' CALL ', t2)
    sents = [s for s in re.split(r'(?<=[.!?])\s+(?=[A-Z])', t3.replace('\n', ' ')) if len(s.split()) >= 3]
    sl = np.array([len(re.findall(r"[A-Za-z][A-Za-z'’]*", s)) for s in sents])
    nums = re.findall(r'(?<![\w.])[−-]?\d+(?:[.,]\d+)*', t3)
    low = [w.lower() for w in words]
    m = {
        'prose words': W,
        'abstract words': len(ab.split()),
        'figures': len(set(re.findall(r'^\s*Fig\. (\d+)\.', body, re.M))),
        'tables': len(set(re.findall(r'^\s*Table (\d+)\.', body, re.M))),
        'references': max([int(x) for x in re.findall(r'^\s*(\d+)\.\s', refs, re.M)] or [0]),
        'caption words (mean)': np.mean([len(c.split()) for c in caps]) if caps else 0,
        'headings': len(heads),
        'sentence length (mean words)': sl.mean(),
        'sentences > 35 words (%)': 100 * (sl > 35).mean(),
        'sentences < 12 words (%)': 100 * (sl < 12).mean(),
        'word length (mean letters)': np.mean([len(w) for w in words]),
        'numbers /1000w': len(nums) * k,
        'percent signs /1000w': t3.count('%') * k,
        'parentheses /1000w': t3.count('(') * k,
        'citation brackets /1000w': len(cites) * k,
        'figure/table/eq callouts /1000w': len(calls) * k,
        'commas /1000w': t3.count(',') * k,
        'semicolons /1000w': t3.count(';') * k,
        'colons /1000w': t3.count(':') * k,
        'hyphenated words /1000w': len(re.findall(r'[A-Za-z]-[A-Za-z]', t3)) * k,
        'en/em dashes /1000w': (t3.count('–') + t3.count('—')) * k,
        'slashes /1000w': t3.count('/') * k,
        'plus-minus, approx, =,<,> /1000w': sum(t3.count(c) for c in '±≈~=<>≤≥') * k,
        'we/our /1000w': sum(w in ('we', 'our') for w in low) * k,
        'acronyms /1000w': len(re.findall(r'\b[A-Z]{2,}s?\b', t3)) * k,
        'however/therefore/thus/moreover /1000w': sum(w in ('however', 'therefore', 'thus', 'moreover', 'furthermore', 'hence', 'consequently', 'additionally') for w in low) * k,
        'was/were /1000w': sum(w in ('was', 'were') for w in low) * k,
        'is/are /1000w': sum(w in ('is', 'are') for w in low) * k,
        'which/that /1000w': sum(w in ('which', 'that') for w in low) * k,
    }
    return m


if __name__ == '__main__':
    ref = sorted(glob.glob(os.path.join(sys.argv[1], '*.pdf'))); ms = sys.argv[2:]
    R = [analyse(f) for f in ref]; keys = list(R[0])
    M = [analyse(f) for f in ms]
    print(f'{len(ref)} reference articles; manuscripts: ' + ', '.join(os.path.basename(f) for f in ms))
    print(f"{'metric':42s} {'ref mean':>9s} {'sd':>7s} {'min':>7s} {'max':>7s} | " + ' | '.join(f"{'value':>7s} {'z':>5s}" for _ in ms))
    for kx in keys:
        v = np.array([r[kx] for r in R], float); mu, sd = v.mean(), v.std(ddof=1)
        row = ' | '.join(f"{m[kx]:7.1f} {((m[kx] - mu) / sd if sd > 0 else 0):+5.1f}" for m in M)
        flag = '  <--' if any(abs((m[kx] - mu) / sd) > 1 for m in M[-1:]) and sd > 0 else ''
        print(f'{kx:42s} {mu:9.1f} {sd:7.1f} {v.min():7.1f} {v.max():7.1f} | {row}{flag}')

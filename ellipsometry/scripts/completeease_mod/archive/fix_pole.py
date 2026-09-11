"""Repair the two generator scripts: the earlier heredoc patch injected literal
tabs/newlines into the UV-pole regex and broke the source.  Replace that whole
block with a plain string substitution (no regex, no escaping hazards)."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import io

CLEAN = (
    "        # UV/IR poles OFF: this model's Einf already absorbs the far-UV part\n"
    "        txt = txt.replace('267.55497728270007\\tT\\t', '0.0\\tF\\t', 1)\n"
    "        txt = txt.replace('11.767447725672241\\tT\\t', '11.767\\tF\\t', 1)\n"
)

for path in [_os.path.join(ELLIPS_OUT, r'make_mod_safe.py'),
             _os.path.join(ELLIPS_OUT, r'make_mod_v3.py')]:
    src = io.open(path, encoding='utf-8').read()
    lines = src.split('\n')
    start = next(i for i, l in enumerate(lines) if 'UV/IR poles OFF' in l)
    end = next(i for i in range(start, len(lines)) if "out=os.path.join(OUTDIR" in lines[i])
    fixed = lines[:start] + CLEAN.rstrip('\n').split('\n') + lines[end:]
    io.open(path, 'w', encoding='utf-8').write('\n'.join(fixed))
    print('repaired %s  (removed %d broken lines)' % (path.split('\\')[-1], end - start))

# syntax check
import py_compile
for path in [_os.path.join(ELLIPS_OUT, r'make_mod_safe.py'),
             _os.path.join(ELLIPS_OUT, r'make_mod_v3.py')]:
    py_compile.compile(path, doraise=True)
    print('syntax OK:', path.split('\\')[-1])

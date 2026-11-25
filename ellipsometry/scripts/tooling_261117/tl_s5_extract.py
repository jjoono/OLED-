"""Pull the model out of a CompleteEASE snapshot (.SEsnap) as a .mod file.

A snapshot is a zip of _DataFile.SE (obfuscated, unreadable without
CompleteEASE) and _FitLog, which carries the full model in .mod syntax, nested
two tabs deeper under 'start_Model1'.  Stripping those two tabs gives a file
whose section skeleton matches a CompleteEASE-saved .mod line for line (checked
against ADNC_ITO_GenOSC_v2_claude.mod; only the substrate block differs, by
design)."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, zipfile

SNAP = _os.path.join(ELLIPS_DATA, 'BPhen_S12_tooling.SEsnap')
OUT = _os.path.join(ELLIPS_OUT, 'S5_template.mod')

def extract(snap=SNAP, out=OUT):
    with zipfile.ZipFile(snap) as z:
        log = z.read('_FitLog').decode('latin-1')
    i = log.index('ALPHA_MULTILAYER\r\n')
    j = log.index('\t\tend_Model Structure\r\n') + len('\t\tend_Model Structure\r\n')
    lines = [l[2:] if l.startswith('\t\t') else l for l in log[i:j].split('\r\n')]
    text = '\r\n' + '\r\n'.join(lines)
    assert all(ord(c) < 128 for c in text)
    open(out, 'w', encoding='utf-8-sig', newline='').write(text)
    return out

if __name__ == '__main__':
    print('->', extract())

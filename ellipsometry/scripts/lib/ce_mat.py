"""Pull CompleteEASE's OWN SI_JAW / NTVE_JAW optical constants straight out of
the .mod file (Mat Table arrays are base64 of gzip of big-endian float32),
so my fit runs in exactly CompleteEASE's frame."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import re, base64, gzip, zlib, numpy as np

MOD = _os.path.join(ELLIPS_DATA, r'VendorA_ITO_GenOSC_6functions.mod')

def decode(b64):
    raw = base64.b64decode(b64)
    try:
        d = gzip.decompress(raw)
    except Exception:
        d = zlib.decompress(raw)
    for dt in ('>f4', '<f4', '>f8', '<f8'):
        a = np.frombuffer(d, dtype=dt)
        if len(a) and np.all(np.isfinite(a)) and np.abs(a).max() < 1e6:
            return a, dt, len(d)
    return np.frombuffer(d, '>f4'), '?', len(d)

def mats(path):
    t = open(path, 'rb').read().decode('utf-8-sig')
    out = {}
    for m in re.finditer(r"'Layer'\t'([A-Z_0-9]+)'.*?start_Mat Table\r\n\t+(\d+)\t\r\n(.*?)end_Mat Table", t, re.S):
        name, npt, body = m.group(1), int(m.group(2)), m.group(3)
        arr = {}
        for tag in ('Wvl', 'e1', 'e2'):
            b = re.search(r"start_%s Array\r\n\t+'([^']+)'" % tag, body)
            if b:
                a, dt, nb = decode(b.group(1))
                arr[tag] = a
                print('  %-10s %-4s n=%4d dtype=%-4s bytes=%d  range %.4g .. %.4g'
                      % (name, tag, len(a), dt, nb, a.min(), a.max()))
        out[name] = (npt, arr)
    return out

if __name__ == '__main__':
    M = mats(MOD)
    np.savez(_os.path.join(ELLIPS_OUT, r'ce_mat.npz'),
             **{'%s_%s' % (k, t): v for k, (n, a) in M.items() for t, v in a.items()})
    print('\nkeys:', list(M))

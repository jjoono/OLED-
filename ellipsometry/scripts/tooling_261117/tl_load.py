"""Loader for tooling.xlsx: one sheet per sample, AnE Psi/Delta interleaved."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, os, re, numpy as np
REPO = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'joint_nk_260910')
sys.path.insert(0, REPO); sys.path.insert(0, os.path.join(REPO, '..', 'lib'))
XL = _os.path.join(ELLIPS_DATA, 'tooling.xlsx')
# another export of the same layout: pass xl=<path> (cached as <stem>.npz)
_cache = lambda xl: _os.path.join(ELLIPS_OUT, _os.path.splitext(_os.path.basename(xl))[0] + '.npz')
CACHE = _cache(XL)

def _read(xl=None):
    import openpyxl
    xl = xl or XL
    wb = openpyxl.load_workbook(xl, data_only=True, read_only=True)
    out = {}
    for sh in wb.sheetnames:
        rows = list(wb[sh].iter_rows(values_only=True))
        h = None
        for i, r in enumerate(rows[:10]):
            for j, v in enumerate(r):
                if isinstance(v, str) and 'Psi (' in v:
                    h, c0, angs = i, j, [float(x) for x in re.findall(r'[\d.]+', v.split('(')[1])]
        rec = []
        for r in rows[h + 1:]:
            try:
                wl = float(r[0])
                psi = [float(r[c0 + 2 * a]) for a in range(len(angs))]
                dl = [float(r[c0 + 1 + 2 * a]) for a in range(len(angs))]
            except (TypeError, ValueError, IndexError):
                continue
            rec.append([wl] + psi + dl)
        out[sh] = np.array(rec)
        out['ANG_' + sh] = np.array(angs)
    np.savez_compressed(_cache(xl), **out)

def sheets(xl=None):
    c = _cache(xl or XL)
    if not os.path.exists(c): _read(xl)
    z = np.load(c)
    return [k for k in z.files if not k.startswith('ANG_')]

def se(sh, lo=245.0, hi=1690.0, step=1, xl=None):
    c = _cache(xl or XL)
    if not os.path.exists(c): _read(xl)
    z = np.load(c); A = z[sh]; ang = z['ANG_' + sh]; na = len(ang)
    m = (A[:, 0] >= lo) & (A[:, 0] <= hi); A = A[m][::step]
    return A[:, 0], A[:, 1:1 + na], A[:, 1 + na:1 + 2 * na], ang

def ncs(P, D):
    p = np.deg2rad(P); d = np.deg2rad(D)
    return np.cos(2 * p), np.sin(2 * p) * np.cos(d), np.sin(2 * p) * np.sin(d)

if __name__ == '__main__':
    for sh in sheets():
        wl, P, D, ang = se(sh)
        print('%-5s n=%d  %.1f-%.1f nm  angles %s' % (sh, len(wl), wl.min(), wl.max(), ang))

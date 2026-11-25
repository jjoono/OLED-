"""Loader for GraHIL_perov.xlsx: 'layer1' = GraHIL on Si, 'layer2' = perovskite
NC film on GraHIL on Si.  65 / 70 / 75 deg, 192.4-1688.3 nm.

Each sheet carries Psi/Delta, N/C/S, intensity and % depolarization per angle.
layer1's Psi/Delta block holds five scans (65, 65, 70, 70, 75 deg: two repeats
at 65 and 70); they are matched to angles through the N,C,S block and the
repeats averaged, their difference kept as a noise estimate.  N, C, S returned
here are the normalised ones rebuilt from Psi/Delta; the depolarization is
returned separately (layer2's exported N,C,S are normalised already).
Negative depolarization (layer2 beyond ~1000 nm, down to -50 %) is unphysical:
those points are instrument artefacts, not sample."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, re, numpy as np
_HERE = _os.path.dirname(_os.path.abspath(__file__))
for p in ('../tooling_261117', '../joint_nk_260910', '../lib'):
    sys.path.insert(0, _os.path.join(_HERE, p))
XL = _os.environ.get('GRAHIL_XLSX', _os.path.join(ELLIPS_DATA, 'GraHIL_perov.xlsx'))
CACHE = _os.path.join(ELLIPS_OUT, 'grahil_perov.npz')
ANG = np.array([65.0, 70.0, 75.0])


def _read():
    import openpyxl
    wb = openpyxl.load_workbook(XL, data_only=True, read_only=True)
    out = {}
    for sh in wb.sheetnames:
        rows = [r for r in wb[sh].iter_rows(values_only=True)]
        hdr = rows[2]; body = [r for r in rows[3:] if r[0] is not None]
        col = {}
        for j, v in enumerate(hdr):
            if isinstance(v, str): col.setdefault(v.split('(')[0].strip(), []).append(j)
        def block(j0):
            """consecutive filled columns starting at j0"""
            n = 0
            while j0 + n < len(body[0]) and body[0][j0 + n] is not None: n += 1
            return np.array([[float(r[j0 + k]) for k in range(n)] for r in body])
        wl = block(col['Wavelength'][0])[:, 0]
        pd = block(col['Psi'][0]); ncs = block(col['N'][0])
        psi_all, del_all = pd[:, 0::2], pd[:, 1::2]
        N = ncs[:, 0::3]
        # the Psi/Delta block may hold repeat scans (layer1: 65, 65, 70, 70, 75): match each to an angle through N
        P = np.sqrt(ncs[:, 0::3]**2 + ncs[:, 1::3]**2 + ncs[:, 2::3]**2)
        which = [int(np.argmin([np.mean(np.abs(N[:, a] / P[:, a] - np.cos(2 * np.deg2rad(psi_all[:, j])))) for a in range(3)])) for j in range(psi_all.shape[1])]
        psi = np.zeros((len(wl), 3)); dl = np.zeros((len(wl), 3)); rp = np.full((len(wl), 3), np.nan); rd = np.full((len(wl), 3), np.nan)
        for a in range(3):
            js = [j for j, w in enumerate(which) if w == a]
            ref = del_all[:, js[0]]
            dd = [ref + ((del_all[:, j] - ref + 180) % 360 - 180) for j in js]
            psi[:, a] = np.mean([psi_all[:, j] for j in js], 0); dl[:, a] = np.mean(dd, 0)
            if len(js) > 1:
                rp[:, a] = np.abs(psi_all[:, js[0]] - psi_all[:, js[1]]); rd[:, a] = np.abs(dd[0] - dd[1])
        out[sh + '_wl'] = wl; out[sh + '_psi'] = psi; out[sh + '_del'] = dl
        out[sh + '_rep_psi'] = rp; out[sh + '_rep_del'] = rd; out[sh + '_nrep'] = np.array([which.count(a) for a in range(3)])
        out[sh + '_I'] = np.stack([block(j)[:, 0] for j in col['Intensity']], 1)
        out[sh + '_dep'] = np.stack([block(j)[:, 0] for j in col['% Depolarization']], 1)
    np.savez_compressed(CACHE, **out)


def load(sh, lo=190.0, hi=1700.0, step=1):
    """-> dict wl, psi, del (deg, repeat scans averaged), dep (% depolarization), I, rep_psi, rep_del
    (|difference| of repeat scans, NaN where there is one scan), N, C, S (normalised, from Psi/Delta), ang."""
    if not _os.path.exists(CACHE): _read()
    z = np.load(CACHE); wl = z[sh + '_wl']; m = (wl >= lo) & (wl <= hi)
    d = {k: z[sh + '_' + k][m][::step] for k in ('wl', 'psi', 'del', 'I', 'dep', 'rep_psi', 'rep_del')}
    p = np.deg2rad(d['psi']); dl = np.deg2rad(d['del'])
    d['N'] = np.cos(2 * p); d['C'] = np.sin(2 * p) * np.cos(dl); d['S'] = np.sin(2 * p) * np.sin(dl)
    d['ang'] = ANG; d['nrep'] = z[sh + '_nrep']
    return d


if __name__ == '__main__':
    for sh in ('layer1', 'layer2'):
        d = load(sh); wl = d['wl']
        print('== %s  %d pts  %.1f-%.1f nm   scans per angle %s' % (sh, len(wl), wl.min(), wl.max(), d['nrep']))
        for a, b in ((192, 250), (250, 400), (400, 700), (700, 1000), (1000, 1700)):
            k = (wl >= a) & (wl < b)
            print('   %4d-%4d nm: depol %% mean %s  min %s  max %s | repeat |dPsi| %s |dDelta| %s' % (a, b, np.round(d['dep'][k].mean(0), 2), np.round(d['dep'][k].min(0), 1),
                  np.round(d['dep'][k].max(0), 1), np.round(np.nanmean(d['rep_psi'][k], 0), 3), np.round(np.nanmean(d['rep_del'][k], 0), 3)))

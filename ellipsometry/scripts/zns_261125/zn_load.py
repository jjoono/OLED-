"""Loader for ZnS_fitting_Sanghee.xlsx: two sheets, each ZnS on a Si wafer (stated to be the same
thickness), 65 / 70 / 75 deg, 192.4-1688.3 nm.  Columns per sheet: Psi/Delta interleaved by angle,
Re/Im(rho), N/C/S, intensity, % depolarization.  N, C, S returned are rebuilt from Psi/Delta."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
import numpy as np

XL = _os.environ.get('ZNS_XLSX', _os.path.join(ELLIPS_DATA, 'ZnS_fitting_Sanghee.xlsx'))
ANG = np.array([65.0, 70.0, 75.0]); SHEETS = ('Sheet1', 'Sheet2')
_C = {}


def raw(sh):
    if sh not in _C:
        import openpyxl
        ws = openpyxl.load_workbook(XL, data_only=True, read_only=True)[sh]
        a = np.array([[np.nan if c is None else c for c in r[:36]] for r in ws.iter_rows(min_row=4, values_only=True) if isinstance(r[0], (int, float))], float)
        _C[sh] = a
    return _C[sh]


def load(sh, lo=0.0, hi=1e9, step=1):
    a = raw(sh); m = (a[:, 0] >= lo) & (a[:, 0] <= hi); a = a[m][::step]
    psi, dl = a[:, [1, 3, 5]], a[:, [2, 4, 6]]; p, d = np.deg2rad(psi), np.deg2rad(dl)
    return dict(wl=a[:, 0], ang=ANG, psi=psi, **{'del': dl}, N=np.cos(2 * p), C=np.sin(2 * p) * np.cos(d), S=np.sin(2 * p) * np.sin(d),
                I=a[:, 28:31], dep=a[:, 33:36])


if __name__ == '__main__':
    for sh in SHEETS:
        d = load(sh); w = d['wl']
        print(sh, len(w), 'points %.1f-%.1f nm' % (w[0], w[-1]))
        for lo, hi in ((192, 300), (300, 400), (400, 700), (700, 1000), (1000, 1200), (1200, 1700)):
            m = (w >= lo) & (w < hi); print('   %4d-%4d nm  depolarization mean %6.2f %%  min %6.2f  max %6.2f | psi 70 deg %5.1f-%5.1f' % (lo, hi, d['dep'][m].mean(), d['dep'][m].min(), d['dep'][m].max(), d['psi'][m, 1].min(), d['psi'][m, 1].max()))
    a, b = load('Sheet1'), load('Sheet2')
    for lo, hi in ((192, 340), (340, 1000), (1000, 1700)):
        m = (a['wl'] >= lo) & (a['wl'] < hi)
        print('Sheet1 - Sheet2, %d-%d nm: rms dN %.4f dC %.4f dS %.4f' % (lo, hi, *[np.sqrt(np.mean((a[k][m] - b[k][m])**2)) for k in 'NCS']))

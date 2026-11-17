"""Fit the 5 ITO/IZO samples in CompleteEASE's OWN frame, full range 192-1688 nm.

  ambient / roughness (Bruggeman 50% film + 50% void) / Gen-Osc film
          / NTVE_JAW 3.0 nm / SI_JAW substrate

Substrate + native-oxide optical constants are the ones decoded out of the .mod
file itself, so a .mod written from this fit reproduces the data on Generate
alone -- no Fit needed.
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import numpy as np, openpyxl, json
from scipy.optimize import least_squares
import ellipsometry_fit as ef
import ce_osc as osc

XLSX = _os.path.join(ELLIPS_DATA, r'se추출.xlsx')
ANG0 = np.array([65.0, 70.0, 75.0])
D_OX = 3.0                                  # nm, NTVE_JAW (fixed, per user: <5 nm)

_m = np.load(_os.path.join(ELLIPS_OUT, r'ce_mat.npz'))
def _mat(tag, wl):
    w = _m[tag + '_Wvl'] / 10.0
    e = np.interp(wl, w, _m[tag + '_e1']) + 1j * np.interp(wl, w, _m[tag + '_e2'])
    N = np.sqrt(e.astype(complex))
    return np.where(N.imag < 0, -N, N)

def load(sheet):
    wb = openpyxl.load_workbook(XLSX, data_only=True, read_only=True)
    ws = wb[sheet]
    rows = []
    for r in ws.iter_rows(min_row=4, values_only=True):
        try: rows.append([float(v) for v in r[:7]])
        except (TypeError, ValueError): pass
    A = np.array(rows)
    return A[:, 0], A[:, 1:7:2], A[:, 2:7:2]

def ncs(psi_deg, del_deg):
    p = np.deg2rad(psi_deg); d = np.deg2rad(del_deg)
    return np.cos(2*p), np.sin(2*p)*np.cos(d), np.sin(2*p)*np.sin(d)

# ---- parameter vector -------------------------------------------------------
# 0 d  1 rough  2 dth  3 Einf  4 rho  5 tau
# 6 TL_A 7 TL_C 8 TL_E0 9 TL_Eg
# 10 Guv_A 11 Guv_Br 12 Guv_En
# 13 Gsb_A 14 Gsb_Br 15 Gsb_En
LO = np.array([30., 0.0, -1.0, 0.5, 5e-5, 1.0,  10., 0.05, 3.2, 2.8,  0.0, 0.5, 4.5,  0.0, 0.2, 1.0])
HI = np.array([90., 15.,  1.0, 4.0, 5e-2, 25.,  600., 5.0, 5.5, 3.9, 60.0, 6.0, 9.0,  5.0, 4.0, 3.2])

def film_N(p, wl):
    E = 1239.841984 / wl
    e1, e2 = osc.drude_rt(E, p[4], p[5])
    a1, a2 = osc.tauc_lorentz(E, p[6], p[7], p[8], p[9]); e1 += a1; e2 += a2
    a1, a2 = osc.gaussian(E, p[10], p[11], p[12]);        e1 += a1; e2 += a2
    a1, a2 = osc.gaussian(E, p[13], p[14], p[15]);        e1 += a1; e2 += a2
    N = np.sqrt((p[3] + e1 + 1j*e2).astype(complex))
    return np.where(N.imag < 0, -N, N)

def forward(p, wl, N_ox, N_si):
    Nf = film_N(p, wl)
    Nr = ef.bruggeman_ema50(Nf, np.ones_like(Nf))
    amb = np.ones_like(Nf)
    out = []
    for a in ANG0 + p[2]:
        rp, rs = ef._tmm(wl, [amb, Nr, Nf, N_ox, N_si],
                         [p[1], p[0], D_OX], a)
        rho = rp / rs
        psi = np.arctan(np.abs(rho)); dl = np.angle(rho)
        out.append((np.cos(2*psi), np.sin(2*psi)*np.cos(dl), np.sin(2*psi)*np.sin(dl)))
    return out

def make_resid(wl, Pm, Dm, N_ox, N_si):
    M = [ncs(Pm[:, i], Dm[:, i]) for i in range(3)]
    def r(p):
        C = forward(p, wl, N_ox, N_si)
        out = []
        for (nm, cm, sm), (nc, cc, sc) in zip(M, C):
            out += [nc-nm, cc-cm, sc-sm]
        return np.concatenate(out)
    return r

def mse(res, npar):
    n = len(res) // 3
    return 1000.0 * np.sqrt(np.sum(res**2) / (3*n - npar))

def fit_sheet(sheet, p0, nstart=6, seed=0):
    wl, Pm, Dm = load(sheet)
    N_ox, N_si = _mat('NTVE_JAW', wl), _mat('SI_JAW', wl)
    r = make_resid(wl, Pm, Dm, N_ox, N_si)
    rng = np.random.default_rng(seed)
    best = None
    for k in range(nstart):
        s = np.array(p0, float)
        if k:
            s = s * (1.0 + 0.25*rng.standard_normal(len(s)))
            s[2] = p0[2] + 0.3*rng.standard_normal()
        s = np.clip(s, LO + 1e-9, HI - 1e-9)
        try:
            b = least_squares(r, s, bounds=(LO, HI), method='trf',
                              x_scale='jac', max_nfev=6000)
        except Exception:
            continue
        if best is None or b.cost < best.cost:
            best = b
    return wl, Pm, Dm, best, mse(best.fun, len(LO)), r

if __name__ == '__main__':
    P0 = [51., 2., 0., 1.0, 5e-4, 5.5, 150., 1.2, 4.2, 3.30, 6.0, 2.5, 6.0, 0.05, 1.0, 1.9]
    wl, Pm, Dm, b, M, r = fit_sheet('#1', P0, nstart=8)
    nm = ['d','rough','dth','Einf','rho','tau','TL_A','TL_C','TL_E0','TL_Eg',
          'Guv_A','Guv_Br','Guv_En','Gsb_A','Gsb_Br','Gsb_En']
    print('#1  MSE = %.3f' % M)
    for n_, v, lo, hi in zip(nm, b.x, LO, HI):
        flag = '  <-- AT BOUND' if (abs(v-lo) < 1e-6*max(1,abs(lo)) or abs(v-hi) < 1e-6*max(1,abs(hi))) else ''
        print('  %-7s %12.6g   [%g, %g]%s' % (n_, v, lo, hi, flag))
    np.save('ce_fit_p1.npy', b.x)

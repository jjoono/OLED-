"""The lab's practice for thin films: the lab model with its n,k fixed, fitted
for thickness with the UV left out.  (Bphen_S12's own 13.8 nm sample matches
above ~450 nm and misses below.)

Every piece, every lower cutoff of the fit window (245 = the full range, 400,
450, 500 nm; upper end 1688 nm), thickness free (and on glass the backside
fraction).  Lab models from Tooling_Evaporator.zip via tl_lab_refit:
NPB_tooling_S8_S13 (uniaxial), TAPC_S11_tooling, BPhen_S12_tooling; glass and
backside as in the TAPC model for every glass piece.  Si oxide: 1.45 nm (the
NPB / TAPC models), and for Bphen also 0 (the Bphen model's own).
Src10 (a real ~10 nm NPB film, SRC10_XLSX) is the thin-film control.

For each fit: d, the MSE inside the window, and the MSE over the full range
with that d (how badly the excluded UV disagrees).  For the 10 nm pieces the
same window is also fitted with no film (Si: a Si/void EMA layer; glass: bare
glass) - inside the window, can the data tell a film from no film?"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_lab_refit as R, tl_load as L, tl_ce_genosc as C, tl_bphen10 as B, tl_glass_global as GG
import jnk_model as M
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
SRC10 = _os.environ.get('SRC10_XLSX', _os.path.join(ELLIPS_DATA, 'Src10_10nm.xlsx'))
CUTS = (245.0, 400.0, 450.0, 500.0)
HI = 1690.0
PIECES = [('Bphen', s) for s in ('2_21', '2_22', '2_23', '2_24', '2_31', '2_32', '2_33', '2_34')] + \
         [('NPB', s) for s in ('1_21', '1_22', '1_23', '1_31', '1_32', '1_33', 'Src10')] + [('TAPC', '1_24'), ('TAPC', '1_34')]
SNAP = {'Bphen': 'BPhen_S12_tooling', 'TAPC': 'TAPC_S11_tooling', 'NPB': 'NPB_tooling_S8_S13'}
R.NOMINAL['Src10'] = 10


def make_set(mat, sh, lo, dox_si):
    glass_m = R.parse('TAPC_S11_tooling')['glass']
    g = R.is_glass(sh) if sh != 'Src10' else False
    wl, P, D, ang = L.se('Sheet1', lo, HI, 2, xl=SRC10) if sh == 'Src10' else L.se(sh, lo, HI, 2)
    E = R.HC / wl
    return dict(sh=sh, glass=g, dox=0.0 if g else dox_si, wl=wl, E=E, ang=ang, meas=L.ncs(P, D),
                eg=R.eps(glass_m, E) if g else None, back=g)


def nofilm(s):
    """Si: Si/void EMA layer on the oxide; glass: bare glass with its backside."""
    if s['glass']:
        r = least_squares(lambda x: R.resid(R.ncs_glass(s['wl'], s['ang'], np.ones(len(s['wl'])), np.ones(len(s['wl'])), 0.0, s['eg'], x[0], True), s),
                          [0.8], bounds=([0.0], [1.0]))
        return float(M.mse_ncs(r.fun, 1)), {'f': float(r.x[0])}
    one = np.ones(len(s['wl']), complex); ox = C.N_ce('NTVE_JAW', s['wl'])**2; si = C.N_ce('SI_JAW', s['wl'])**2
    def res(x):
        lay = [one] + ([ox] if s['dox'] > 0 else []) + [B.bruggeman_eps(si, one, 1 - x[1]), si]
        ds = ([s['dox']] if s['dox'] > 0 else []) + [x[0]]
        return R.resid(B.ncs_regions(s['wl'], s['ang'], [(1.0, lay, lay, ds)], GG.CONV['normalised']), s)
    best = None
    for t0, f0 in ((5.0, 0.3), (12.0, 0.3), (3.0, 0.6), (15.0, 0.7)):
        r = least_squares(res, [t0, f0], bounds=([0.0, 0.0], [40.0, 0.95]))
        if best is None or r.cost < best.cost: best = r
    return float(M.mse_ncs(best.fun, 2)), {'t': float(best.x[0]), 'fv': float(best.x[1])}


def job(mat, sh, dox_si):
    film = R.parse(SNAP[mat])['film']; out = []
    full = make_set(mat, sh, CUTS[0], dox_si)
    for lo in CUTS:
        s = make_set(mat, sh, lo, dox_si)
        r = R.fit_thickness(film, s, bounded=False)
        d = float(r.x[0]); f = float(r.x[1]) if s['back'] else 0.0
        mse_full = float(M.mse_ncs(R.resid(R.model_ncs(film, full, d, f), full), len(r.x)))
        row = dict(mat=mat, sh=sh, dox_si=dox_si, lo=lo, d=d, f=f if s['back'] else None,
                   mse_window=float(M.mse_ncs(r.fun, len(r.x))), mse_full_with_this_d=mse_full)
        if R.NOMINAL[sh] == 10:
            row['nofilm_mse_window'], row['nofilm'] = nofilm(s)
        out.append(row)
    return out


def main():
    from multiprocessing import Pool
    jobs = [(m, s, 1.45) for m, s in PIECES] + [('Bphen', s, 0.0) for m, s in PIECES if m == 'Bphen']
    with Pool(4) as pool:
        res = [r for rr in pool.starmap(job, jobs) for r in rr]
    json.dump(res, open(W + 'lab_visible.json', 'w'), indent=1)
    print('%-6s %-6s %-5s %4s | %s' % ('mat', 'piece', 'NTVE', 'nom', ' | '.join('from %3.0f nm: d (MSE win / full)%s' % (c, '') for c in CUTS)))
    for (m, s, dox) in jobs:
        rows = [r for r in res if r['sh'] == s and r['mat'] == m and r['dox_si'] == dox]
        line = '%-6s %-6s %4.2f %4d | ' % (m, s, dox, R.NOMINAL[s]) + ' | '.join('%6.2f (%5.1f / %5.1f)' % (r['d'], r['mse_window'], r['mse_full_with_this_d']) for r in rows)
        if 'nofilm_mse_window' in rows[0]:
            line += '   || no film MSE in window: ' + ' / '.join('%.1f' % r['nofilm_mse_window'] for r in rows)
        print(line)


if __name__ == '__main__':
    main()

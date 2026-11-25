"""ZnS on Si, two pieces said to be the same thickness (expected 100-250 nm).

    air / roughness (ZnS + void 50:50) d_r / ZnS d / native oxide (NTVE_JAW, D_OX) / Si (SI_JAW)

Stage A  transparent range (430-1000 nm), Sellmeier ZnS: every thickness from 40 to 400 nm tried in turn
         (the n.d alias check), each piece on its own.
Stage C  205-1000 nm, Kramers-Kronig-consistent dispersion (Gaussians on a fixed grid from 3.0 eV up, UV
         pole, eps_inf), started from stage A.  Variants:
           one piece at a time / both pieces with one dispersion (thickness and roughness per piece)
           graded     (eps - 1) running linearly through the depth: top (1 - g), bottom (1 + g) x mean
           nonuniform lateral thickness spread, fitted against N,C,S x sqrt(1 - depolarization)
           oxide      native-oxide thickness free
Beyond 1000 nm the wafer is transparent (back reflection): not fitted.  -> zns_fits.json  (zn_report.py writes the n,k table)"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
_H = _os.path.dirname(_os.path.abspath(__file__)); sys.path.insert(0, _H); sys.path.insert(0, _os.path.join(_H, '../grahil_perov_261124'))
import zn_load as L, gp_model as P
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
LO, HI, STEP = 205.0, 1000.0, 2
D_OX = 1.5
_c = np.r_[np.arange(3.0, 4.41, 0.07), np.arange(4.6, 7.41, 0.2)]
CEN, WID = _c, 1.6 * np.gradient(_c); E_UV = 10.0; SMOOTH = 0.003; M_GRADE = 8


def setup(sheets, lo=LO, hi=HI, step=STEP):
    D = [L.load(s, lo, hi, step) for s in sheets]; wl = D[0]['wl']
    return dict(D=D, wl=wl, ang=L.ANG, air=P.n_air(wl), ox=P.n_ox(wl), si=P.n_si(wl), kb=P.KKBasis(wl, CEN, WID, E_UV), n=len(sheets))


def stack(S, N, d, d_r, g=0.0, dox=D_OX, sig=0.0):
    """N: ZnS index (mean).  -> N, C, S of one piece"""
    chi = N**2 - 1; m = M_GRADE if g != 0 else 1
    lay = [np.sqrt(1 + (1 + g * (2 * (i + 0.5) / m - 1)) * chi) for i in range(m)]               # top first
    top = lay[0]
    Ls = [S['air'], P.bruggeman(top, S['air'], 0.5)] + lay + [S['ox'], S['si']]; ds = [d_r] + [d / m] * m + [dox]
    return P.ncs(S['wl'], S['ang'], Ls, ds, sig, vary=list(range(1, m + 1)))


# ------------------------------------------------------------------ stage A
def stage_a(sheet):
    S = setup([sheet], 430.0, HI, 3); wl = S['wl']; meas = P.data_ncs(S['D'][0]); out = []
    sel = lambda x: np.sqrt(x[0] + x[1] * wl**2 / (wl**2 - x[2]**2) + 0j)
    for d in np.arange(40.0, 401.0, 4.0):
        best = None
        for x0 in ([3.0, 2.2, 250.0, 3.0], [1.5, 3.2, 280.0, 8.0], [2.0, 1.5, 200.0, 1.0]):
            r = least_squares(lambda x: P.resid(stack(S, sel(x), d, x[3]), meas), x0, bounds=([1, 0, 120, 0], [6, 6, 330, 40]), x_scale='jac', max_nfev=80)
            if best is None or r.cost < best.cost: best = r
        out.append(dict(d=float(d), mse=float(P.mse_ncs(best.fun, 5)), x=[float(v) for v in best.x], n550=float(np.sqrt(best.x[0] + best.x[1] * 550.0**2 / (550.0**2 - best.x[2]**2)))))
    k = int(np.argmin([o['mse'] for o in out])); b = out[k]
    r = least_squares(lambda x: P.resid(stack(S, sel(x[1:4]), x[0], x[4]), meas), [b['d']] + b['x'], bounds=([30, 1, 0, 120, 0], [420, 6, 6, 330, 40]), x_scale='jac', max_nfev=300)
    return dict(sheet=sheet, scan=out, d=float(r.x[0]), sellmeier=[float(v) for v in r.x[1:4]], d_r=float(r.x[4]), mse=float(P.mse_ncs(r.fun, 5)),
                n550=float(np.sqrt(r.x[1] + r.x[2] * 550.0**2 / (550.0**2 - r.x[3]**2))))


# ------------------------------------------------------------------ stage C
def init_basis(kb, sel_x):
    """amplitudes reproducing the Sellmeier index of stage A in the transparent range (start only)"""
    from scipy.optimize import lsq_linear
    wl = kb.wl; m = wl > 430; n = np.sqrt(sel_x[0] + sel_x[1] * wl[m]**2 / (wl[m]**2 - sel_x[2]**2)); A = np.c_[np.ones(m.sum()), kb.pole[m], kb.G[m].real]
    D2 = np.c_[np.zeros((kb.n - 2, 2)), np.diff(np.eye(kb.n), 2, axis=0)] * 0.3
    r = lsq_linear(np.r_[A, D2], np.r_[n**2, np.zeros(kb.n - 2)], bounds=(np.r_[1.0, 0.0, np.zeros(kb.n)], np.r_[4.0, 400.0, np.full(kb.n, 30.0)]))
    return r.x


def fit(sheets, a_res, opts=(), x_from=None):
    S = setup(sheets); kb = S['kb']; nb = kb.n; ns = S['n']; dep = 'nonuniform' in opts
    meas = [P.data_ncs(d, with_dep=dep) for d in S['D']]
    names = ['d%d' % (i + 1) for i in range(ns)] + ['d_r%d' % (i + 1) for i in range(ns)] + [o for o in ('g', 'sig', 'dox') if {'g': 'graded', 'sig': 'nonuniform', 'dox': 'oxide'}[o] in opts]
    lo = [30.0] * ns + [0.0] * ns + [{'g': -0.5, 'sig': 0.0, 'dox': 0.0}[k] for k in names[2 * ns:]]
    hi = [420.0] * ns + [40.0] * ns + [{'g': 0.5, 'sig': 0.2, 'dox': 8.0}[k] for k in names[2 * ns:]]
    x0 = [a['d'] for a in a_res] + [max(a['d_r'], 1.0) for a in a_res] + [{'g': 0.02, 'sig': 0.01, 'dox': D_OX}[k] for k in names[2 * ns:]]
    np_ = len(names); xb = init_basis(kb, a_res[0]['sellmeier']) if x_from is None else x_from
    lo = np.r_[lo, 1.0, 0.0, np.zeros(nb)]; hi = np.r_[hi, 4.0, 400.0, np.full(nb, 30.0)]; x0 = np.clip(np.r_[x0, xb], lo + 1e-9, hi - 1e-9)
    D2 = np.diff(np.eye(nb), 2, axis=0)
    def parts(x):
        v = dict(zip(names, x[:np_])); N = kb.N(x[np_ + 2:], x[np_ + 1], x[np_]); return v, N
    def res(x, reg=True):
        v, N = parts(x)
        r = [P.resid(stack(S, N, v['d%d' % (i + 1)], v['d_r%d' % (i + 1)], v.get('g', 0.0), v.get('dox', D_OX), v.get('sig', 0.0)), meas[i]) for i in range(ns)]
        return np.r_[np.concatenate(r), SMOOTH * (D2 @ x[np_ + 2:])] if reg else r
    r = least_squares(res, x0, bounds=(lo, hi), x_scale='jac', max_nfev=150)
    v, N = parts(r.x); rr = res(r.x, False); wl = S['wl']
    J = r.jac; s2 = np.sum(r.fun**2) / max(r.fun.size - r.x.size, 1); ci = 1.645 * np.sqrt(np.clip(np.diag(np.linalg.pinv(J.T @ J) * s2), 0, None))
    band = lambda q, p_, e: float(1e3 * np.sqrt(np.mean(q.reshape(3, len(wl), 3)[:, (wl >= p_) & (wl < e), :]**2)))
    return dict(sheets=list(sheets), opts=list(opts), names=names, x=[float(t) for t in r.x], struct={k: float(t) for k, t in v.items()}, ci={k: float(c) for k, c in zip(names, ci[:np_])},
                mse=[float(P.mse_ncs(q, len(r.x) / ns)) for q in rr], mse_all=float(P.mse_ncs(np.concatenate(rr), len(r.x))),
                bands=[{'205-330': band(q, 205, 330), '330-430': band(q, 330, 430), '430-1000': band(q, 430, 1001)} for q in rr],
                nk={str(w): [float(np.interp(w, wl, N.real)), float(np.interp(w, wl, N.imag))] for w in (250, 300, 340, 360, 400, 450, 550, 633, 800, 1000)})


def nk_of(r, wl):
    kb = P.KKBasis(wl, CEN, WID, E_UV); k = len(r['names']); x = np.array(r['x']); return kb.N(x[k + 2:], x[k + 1], x[k])


def job(a):
    return fit(*a)


def main():
    from multiprocessing import Pool
    with Pool(4) as pool:
        A = pool.map(stage_a, L.SHEETS)
        for a in A:
            sc = sorted(a['scan'], key=lambda o: o['mse']); loc = [o for i, o in enumerate(a['scan']) if 0 < i < len(a['scan']) - 1 and o['mse'] < a['scan'][i - 1]['mse'] and o['mse'] < a['scan'][i + 1]['mse']]
            print('stage A %s: best d %.1f nm, roughness %.1f nm, n(550) %.3f, MSE %.2f | local minima of the thickness scan: %s' % (a['sheet'], a['d'], a['d_r'], a['n550'], a['mse'], ', '.join('%.0f nm (MSE %.1f, n %.2f)' % (o['d'], o['mse'], o['n550']) for o in sorted(loc, key=lambda o: o['mse'])[:5])))
        jobs = [([s], [a], ()) for s, a in zip(L.SHEETS, A)] + [(list(L.SHEETS), A, o) for o in ((), ('graded',), ('nonuniform',), ('oxide',), ('graded', 'nonuniform'), ('graded', 'nonuniform', 'oxide'))]
        R = pool.map(job, jobs)
    json.dump(dict(stage_a=A, fits=R), open(W + 'zns_fits.json', 'w'), indent=1)
    print('%-34s %-13s | %-17s %-13s | %6s %6s %6s | %s' % ('fit', 'MSE (pieces)', 'thickness (nm)', 'roughness', 'g', 'sig %', 'oxide', 'n,k at 340 / 400 / 550 / 633; bands of piece 1'))
    for r in R:
        v = r['struct']; ns = len(r['sheets'])
        print('%-34s %-13s | %-17s %-13s | %6.3f %6.2f %6.2f | %s; %s' % ('+'.join(r['sheets']) + (' ' + '+'.join(r['opts']) if r['opts'] else ''), ' / '.join('%.2f' % m for m in r['mse']),
              ' / '.join('%.2f' % v['d%d' % (i + 1)] for i in range(ns)), ' / '.join('%.2f' % v['d_r%d' % (i + 1)] for i in range(ns)), v.get('g', 0), 100 * v.get('sig', 0), v.get('dox', D_OX),
              '  '.join('%.3f %.3f' % tuple(r['nk'][w]) for w in ('340', '400', '550', '633')), ' '.join('%.1f' % b for b in r['bands'][0].values())))


if __name__ == '__main__':
    main()

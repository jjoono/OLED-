"""Bphen 10 nm (2_23, 2_24) the CompleteEASE way: fit the visible only, then
expand into the UV with a B-spline film.  NTVE_JAW fixed at 1.5 nm.

Controls run through the identical procedure: Src10 (a real ~10 nm NPB film,
SRC10_XLSX) and 2_21 (Bphen ~45 nm, same run as 2_23 / 2_24).

Step A, visible only (420-1000 nm, above the Bphen onset and Si E1):
  * the material's Gen-Osc fixed, d free
  * a transparent Cauchy film, d free, and its MSE profile against d -
    on a ~10 nm film the visible fixes roughly (n^2 - 1) d, not d
  * no film: a Si/void EMA layer (section 3's reading), for comparison
Step B, wavelength expansion to 245 nm at the step-A Gen-Osc thickness, the
film a B-spline in photon energy (nodes every 0.2 eV, 1.1-5.3 eV):
  * free:  n and k each a cubic B-spline, k >= 0 (CompleteEASE's plain B-spline)
  * KK:    eps2 from non-negative Gaussians on the same 0.2 eV grid, eps1 their
           exact KK partner plus a UV pole (CompleteEASE's KK-consistent option)
  The range opens 20 nm at a time, each step refitted from the last; a weak
  second-difference penalty carries nodes that have no data yet.
What decides the question is not the final MSE (a B-spline fits almost
anything) but whether the extracted k shows the material's bands or the
substrate's (Si E1 3.4 eV / 365 nm, E2 4.25 eV / 292 nm)."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, zlib, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_ce_genosc as C, tl_glass_global as GG, tl_bphen10 as B, tl_genosc as GO, tl_ce_thickness as T
import jnk_model as M, ce_osc as osc
from scipy.optimize import least_squares
from scipy.interpolate import BSpline

W = _os.path.join(ELLIPS_OUT, '')
SRC10 = _os.environ.get('SRC10_XLSX', _os.path.join(ELLIPS_DATA, 'Src10_10nm.xlsx'))
HC = 1239.841984
DOX = 1.5                                   # nm, NTVE_JAW fixed
VIS = (420.0, 1000.0)
PIECES = {'2_23': 'Bphen', '2_24': 'Bphen', 'Src10': 'NPB', '2_21': 'Bphen'}
KNOTS = np.arange(1.1, 5.31, 0.2)
_t = np.r_[[KNOTS[0]] * 3, KNOTS, [KNOTS[-1]] * 3]
NB = len(_t) - 4
GC = np.arange(1.3, 5.51, 0.2)              # Gaussian centres of the KK basis, eV
GBR = 0.35                                  # their FWHM, eV
PEN = 0.02                                  # second-difference penalty weight
EXPAND = [400.0, 380.0, 360.0, 340.0, 320.0, 300.0, 280.0, 260.0, 245.0]


def basis(wl):
    E = np.clip(HC / wl, KNOTS[0], KNOTS[-1])
    return BSpline.design_matrix(E, _t, 3).toarray()


def data(piece, lo, hi):
    sh, xl = ('Sheet1', SRC10) if piece == 'Src10' else (piece, None)
    wl, P, D, ang = L.se(sh, lo, hi, 2, xl=xl)
    return wl, L.ncs(P, D), ang


def ncs_film(wl, ang, ef, d, dox=DOX, extra=None):
    air = np.ones(len(wl), complex); ox = C.N_ce('NTVE_JAW', wl)**2; si = C.N_ce('SI_JAW', wl)**2
    lay, ds = ([air, ef, ox, si], [d, dox]) if extra is None else extra(air, ox, si)
    return B.ncs_regions(wl, ang, [(1.0, lay, lay, ds)], GG.CONV['normalised'])


def resid(model, meas):
    return np.concatenate([(model[k] - meas[k]).ravel() for k in range(3)])


def bands(r, wl, ang):
    rr = r.reshape(3, len(wl), len(ang))
    f = lambda a, b: float(1e3 * np.sqrt(np.mean(rr[:, (wl >= a) & (wl < b), :]**2))) if np.any((wl >= a) & (wl < b)) else None
    return {'245-350': f(245, 350), '350-420': f(350, 420), '420-1000': f(420, 1001)}


def multistart(res, x0, lo, hi, seed, n=12):
    rng = np.random.default_rng(seed); best = None
    lo, hi = np.array(lo, float), np.array(hi, float)
    for k in range(n):
        s = np.array(x0, float) if k == 0 else lo + rng.random(len(lo)) * (hi - lo)
        r = least_squares(res, np.clip(s, lo + 1e-6, hi - 1e-6), bounds=(lo, hi), max_nfev=500)
        if best is None or r.cost < best.cost: best = r
    return best


# ------------------------------------------------------------------ step A
def step_a(piece):
    mat = PIECES[piece]; wl, meas, ang = data(piece, *VIS)
    p = np.array(T.disp(mat)['p'])
    ef = GO.eps_genosc(p, wl)
    out = {}
    r = multistart(lambda x: resid(ncs_film(wl, ang, ef, x[0]), meas), [10.0], [0.0], [80.0], 1)
    out['genosc'] = dict(d=float(r.x[0]), mse=float(M.mse_ncs(r.fun, 1)))
    cau = lambda x: (x[1] + x[2] / (wl / 1000)**2 + 0j)**2
    r = multistart(lambda x: resid(ncs_film(wl, ang, cau(x), x[0]), meas), [10.0, 1.7, 0.01], [0.0, 1.3, 0.0], [80.0, 2.6, 0.1], 2)
    J = r.jac; cov = np.linalg.pinv(J.T @ J) * (r.fun @ r.fun) / (len(r.fun) - 3)
    out['cauchy'] = dict(d=float(r.x[0]), A=float(r.x[1]), B=float(r.x[2]), mse=float(M.mse_ncs(r.fun, 3)),
                         d_ci90=float(1.645 * np.sqrt(cov[0, 0])), A_ci90=float(1.645 * np.sqrt(cov[1, 1])))
    prof = []
    for dd in (3.0, 5.0, 7.0, 9.0, 12.0, 16.0, 20.0):
        rr = multistart(lambda x: resid(ncs_film(wl, ang, cau([dd] + list(x)), dd), meas), [1.7, 0.01], [1.0, 0.0], [3.5, 0.2], 3, n=4)
        prof.append(dict(d=dd, A=float(rr.x[0]), mse=float(M.mse_ncs(rr.fun, 2))))
    out['cauchy_profile'] = prof
    def nofilm(x):
        ext = lambda air, ox, si: ([air, ox, B.bruggeman_eps(si, air, 1 - x[1]), si], [DOX, x[0]])
        return resid(ncs_film(wl, ang, None, 0, extra=ext), meas)
    r = multistart(nofilm, [8.0, 0.3], [0.0, 0.0], [30.0, 0.9], 4)
    out['no film'] = dict(t=float(r.x[0]), fv=float(r.x[1]), mse=float(M.mse_ncs(r.fun, 2)))
    return out


# ------------------------------------------------------------------ step B
def eps_free(c, Bm):
    n = Bm @ c[:NB]; k = Bm @ c[NB:]
    return (n + 1j * k)**2

def eps_kk(c, wl):
    E = HC / wl; e = 1.0 + c[0] / (c[1]**2 - E**2) + 0j
    for A, En in zip(c[2:], GC):
        a, b = osc.gaussian(E, A, GBR, En); e = e + a + 1j * b
    return e


def expand(piece, d, kind):
    mat = PIECES[piece]; p = np.array(T.disp(mat)['p'])
    wl_all = np.linspace(245.0, 1000.0, 400)
    # start: the material's Gen-Osc projected onto the basis
    if kind == 'free':
        N0 = np.sqrt(GO.eps_genosc(p, wl_all).astype(complex)); Bm = basis(wl_all)
        c = np.r_[np.linalg.lstsq(Bm, N0.real, rcond=None)[0], np.clip(np.linalg.lstsq(Bm, N0.imag, rcond=None)[0], 0, None)]
        lo = np.r_[np.full(NB, 1.0), np.zeros(NB)]; hi = np.r_[np.full(NB, 3.5), np.full(NB, 3.0)]
        D2 = np.diff(np.eye(NB), 2, axis=0); pen = lambda c: PEN * np.r_[D2 @ c[:NB], D2 @ c[NB:]]
    else:
        c = np.r_[p[0], p[1], np.full(len(GC), 0.01)]
        lo = np.r_[0.0, 6.0, np.zeros(len(GC))]; hi = np.r_[1000.0, 15.0, np.full(len(GC), 10.0)]
        D2 = np.diff(np.eye(len(GC)), 2, axis=0); pen = lambda c: PEN * (D2 @ c[2:])
    c = np.clip(c, lo + 1e-9, hi - 1e-9)
    steps = []
    for lo_wl in [VIS[0]] + EXPAND:
        wl, meas, ang = data(piece, lo_wl, VIS[1])
        Bm = basis(wl) if kind == 'free' else None
        ef = (lambda c: eps_free(c, Bm)) if kind == 'free' else (lambda c: eps_kk(c, wl))
        res = lambda c: np.r_[resid(ncs_film(wl, ang, ef(c), d), meas), pen(c)]
        r = least_squares(res, c, bounds=(lo, hi), max_nfev=300, x_scale='jac'); c = r.x
        rd = resid(ncs_film(wl, ang, ef(c), d), meas)
        steps.append(dict(lo=lo_wl, mse=float(M.mse_ncs(rd, len(c)))))
    wl, meas, ang = data(piece, 245.0, VIS[1])
    epsf = eps_free(c, basis(wl_all)) if kind == 'free' else eps_kk(c, wl_all)
    Nf = np.sqrt(epsf.astype(complex)); Nf = np.where(Nf.imag < 0, -Nf, Nf)
    rd = resid(ncs_film(wl, ang, eps_free(c, basis(wl)) if kind == 'free' else eps_kk(c, wl), d), meas)
    Nref = np.sqrt(GO.eps_genosc(p, wl_all).astype(complex))
    kk = Nf.imag; pk = [float(wl_all[i]) for i in range(1, len(kk) - 1) if kk[i] > kk[i - 1] and kk[i] >= kk[i + 1] and kk[i] > 0.02]
    pick = lambda arr, w: float(np.interp(w, wl_all, arr))
    return dict(piece=piece, kind=kind, d=d, steps=steps, mse=float(M.mse_ncs(rd, len(c))), bands=bands(rd, wl, ang),
                k_peaks_nm=pk, nk={str(int(w)): [pick(Nf.real, w), pick(Nf.imag, w)] for w in (260, 280, 300, 330, 365, 400, 450, 550, 700)},
                nk_ref={str(int(w)): [pick(Nref.real, w), pick(Nref.imag, w)] for w in (260, 280, 300, 330, 365, 400, 450, 550, 700)},
                wl=wl_all.tolist(), n=Nf.real.tolist(), k=Nf.imag.tolist())


def genosc_full(piece, d):
    """The step-A Gen-Osc film carried into the UV unchanged."""
    p = np.array(T.disp(PIECES[piece])['p']); wl, meas, ang = data(piece, 245.0, VIS[1])
    rd = resid(ncs_film(wl, ang, GO.eps_genosc(p, wl), d), meas)
    return dict(mse=float(M.mse_ncs(rd, 1)), bands=bands(rd, wl, ang))


def run_piece(piece):
    a = step_a(piece); d = a['genosc']['d']
    return dict(piece=piece, mat=PIECES[piece], step_a=a, genosc_full=genosc_full(piece, d),
                free=expand(piece, d, 'free'), kk=expand(piece, d, 'kk'))


def main():
    from multiprocessing import Pool
    with Pool(4) as pool:
        res = pool.map(run_piece, list(PIECES))
    json.dump(res, open(W + 'bphen10_bspline.json', 'w'), indent=1)
    for r in res:
        a = r['step_a']
        print('\n== %s (%s), NTVE %.1f nm' % (r['piece'], r['mat'], DOX))
        print('  visible %d-%d nm:  %s Gen-Osc d = %.2f nm MSE %.2f | Cauchy d = %.2f +- %.2f, A = %.3f +- %.3f, MSE %.2f | no film MSE %.2f (t %.1f, void %.2f)'
              % (VIS[0], VIS[1], r['mat'], a['genosc']['d'], a['genosc']['mse'], a['cauchy']['d'], a['cauchy']['d_ci90'],
                 a['cauchy']['A'], a['cauchy']['A_ci90'], a['cauchy']['mse'], a['no film']['mse'], a['no film']['t'], a['no film']['fv']))
        print('  Cauchy MSE vs fixed d: ' + '  '.join('%g nm %.2f (A %.2f)' % (q['d'], q['mse'], q['A']) for q in a['cauchy_profile']))
        g = r['genosc_full']
        print('  Gen-Osc carried to 245 nm unchanged: MSE %.2f | rms x1e3 %s' % (g['mse'], g['bands']))
        for kind in ('free', 'kk'):
            e = r[kind]
            print('  B-spline %-4s expanded to 245 nm: MSE %.2f | rms x1e3 %s | k peaks at %s nm'
                  % (kind, e['mse'], e['bands'], ', '.join('%.0f' % x for x in e['k_peaks_nm']) or 'none'))
            print('      n,k (spline / %s Gen-Osc): ' % r['mat'] + '  '.join('%s: %.2f,%.2f / %.2f,%.2f' % (w, *e['nk'][w], *e['nk_ref'][w]) for w in e['nk']))


if __name__ == '__main__':
    main()

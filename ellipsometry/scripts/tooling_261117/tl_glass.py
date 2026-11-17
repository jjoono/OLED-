"""Glass pieces: identity and thickness with the Si-derived dispersions held fixed.

Front stack air / film / glass coherent.  The glass index is unknown (no bare
reference in this set), so a 2-term Cauchy is fitted with it.  The back surface
is optional: f is the fraction of the back-surface beam that reaches the
detector, and its contributions add incoherently (Mueller sum in N,C,S).
f = 0 means the backside was roughened or index-matched."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_osc as O
import jnk_ref as R, jnk_model as M
from tmm_fix import tmm
from scipy.optimize import least_squares
W = _os.path.join(ELLIPS_OUT, '')
OSC = json.load(open(W + 'osc.json'))

def n_glass(wl, A, B):
    return A + B / (wl / 1000.0)**2 + 0j

def ncs_glass(wl, Nf, d, ng, ang, f=0.0, nb=20):
    """N,C,S of film-on-thick-glass with a fraction f of the incoherent backside beams."""
    air = R.n_air(wl); ng = np.asarray(ng, complex) * np.ones(len(wl))
    out = []
    for th in ang:
        rp, rs = tmm(wl, [air, Nf, ng], [d], th)
        Ip = np.abs(rp)**2; Is = np.abs(rs)**2; X = np.conj(rp) * rs
        if f > 0:
            thg = np.rad2deg(np.arcsin(np.sin(np.deg2rad(th)) / ng.real))
            tdp, tds = M.t_amp(wl, [air, Nf, ng], [d], th)
            tup, tus = M.t_amp(wl, [ng, Nf, air], [d], thg)
            rip, ris = tmm(wl, [ng, Nf, air], [d], thg)          # front stack seen from inside
            rbp, rbs = tmm(wl, [ng, air], [], thg)                # bare back surface
            ap = tdp * rbp * tup; as_ = tds * rbs * tus
            for _ in range(nb):
                Ip = Ip + f * np.abs(ap)**2; Is = Is + f * np.abs(as_)**2; X = X + f * np.conj(ap) * as_
                ap = ap * rip * rbp; as_ = as_ * ris * rbs
        I = Ip + Is
        out.append(((Is - Ip) / I, 2 * X.real / I, 2 * X.imag / I))
    N = np.stack([o[0] for o in out], 1); C = np.stack([o[1] for o in out], 1); S = np.stack([o[2] for o in out], 1)
    return N, C, S

def as_reported(N, C, S):
    """What a Psi/Delta export carries: Psi = acos(N)/2, Delta = atan2(S, C)."""
    dl = np.arctan2(S, C); s2 = np.sqrt(np.clip(1 - N**2, 0, 1))
    return N, s2 * np.cos(dl), s2 * np.sin(dl)

def fit(sh, mat, with_back=False, lo=245.0, hi=1000.0, d0=None):
    wl, P, D, ang = L.se(sh, lo, hi, step=2)
    meas = L.ncs(P, D)
    Nf = O.N_of(np.array(OSC[mat]['q']), wl)
    def res(x):
        d, A, B = x[:3]; f = x[3] if with_back else 0.0
        N, C, S = as_reported(*ncs_glass(wl, Nf, d, n_glass(wl, A, B), ang, f))
        return np.concatenate([(N - meas[0]).ravel(), (C - meas[1]).ravel(), (S - meas[2]).ravel()])
    best = None
    starts = [d0] if d0 else list(np.arange(2.0, 120.0, 4.0))
    for ds in starts:
        x0 = [ds, 1.52, 0.004] + ([0.5] if with_back else [])
        lo_b = [0.5, 1.40, -0.01] + ([0.0] if with_back else [])
        hi_b = [300.0, 1.70, 0.03] + ([1.0] if with_back else [])
        r = least_squares(res, x0, bounds=(lo_b, hi_b), max_nfev=400 if not d0 else 2000)
        if best is None or r.cost < best.cost: best = r
    if not d0:     # polish the winner
        best = least_squares(res, best.x, bounds=(lo_b, hi_b), max_nfev=3000)
    return best.x, M.mse_ncs(best.fun, len(best.x))

if __name__ == '__main__':
    from multiprocessing import Pool
    shs = [s for s in L.sheets() if s.split('_')[1].startswith('3')]
    jobs = [(sh, m) for sh in shs for m in ('NPB', 'TAPC', 'Bphen')]
    with Pool(4) as p:
        res = p.starmap(fit, jobs)
    tab = {}
    for (sh, m), (x, mse) in zip(jobs, res):
        tab.setdefault(sh, {})[m] = dict(x=[float(v) for v in x], mse=float(mse))
    print('front surface only (f = 0), glass Cauchy fitted, film dispersion fixed from Si')
    print('%-5s | %-22s | %-22s | %-22s | %s' % ('sheet', 'NPB  d / MSE', 'TAPC  d / MSE', 'Bphen  d / MSE', 'glass n(550) at best'))
    for sh in shs:
        t = tab[sh]; b = min(t, key=lambda m: t[m]['mse'])
        x = t[b]['x']
        print('%-5s | %s | %s' % (sh, ' | '.join('%7.2f / %8.2f%s' % (t[m]['x'][0], t[m]['mse'], ' *' if m == b else '  ')
                                      for m in ('NPB', 'TAPC', 'Bphen')),
              '%.4f' % n_glass(np.array([550.]), x[1], x[2]).real[0]))
    json.dump(tab, open(W + 'glass_f0.json', 'w'), indent=1)

"""Was the glass piece measured from the side WITHOUT the film?

Bphen allowed on the front surface (d_front), on the back surface (d_back) or
both; the back film is seen only through the incoherent back-surface beams, at
the fraction f fitted for the other glass pieces (and once with f free).
If the deposit sits on the back, d_back lands near the deposited thickness and
the front goes to zero."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, zlib, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_genosc as GO, tl_glass_global as GG
import jnk_model as M
from tmm_fix import tmm
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')

def ncs_two_sided(wl, Nf, d_front, d_back, ng, ang, f, nb=12):
    air = np.ones(len(wl), complex); ng = ng * np.ones(len(wl))
    fl = lambda d: ([Nf], [d]) if d > 1e-3 else ([], [])
    out = []
    for th in ang:
        Lf, df = fl(d_front); Lb, db = fl(d_back)
        rp, rs = tmm(wl, [air] + Lf + [ng], df, th)
        Ip = np.abs(rp)**2; Is = np.abs(rs)**2; X = np.conj(rp) * rs
        thg = np.rad2deg(np.arcsin(np.sin(np.deg2rad(th)) / ng.real))
        tdp, tds = M.t_amp(wl, [air] + Lf + [ng], df, th)
        tup, tus = M.t_amp(wl, [ng] + Lf[::-1] + [air], df[::-1], thg)
        rip, ris = tmm(wl, [ng] + Lf[::-1] + [air], df[::-1], thg)     # front, from inside
        rbp, rbs = tmm(wl, [ng] + Lb + [air], db, thg)                 # back surface, with its film
        ap = tdp * rbp * tup; as_ = tds * rbs * tus
        for _ in range(nb):
            Ip = Ip + f * np.abs(ap)**2; Is = Is + f * np.abs(as_)**2; X = X + f * np.conj(ap) * as_
            ap = ap * rip * rbp; as_ = as_ * ris * rbs
        I = Ip + Is
        out.append(GG.norm((Is - Ip) / I, 2 * X.real / I, 2 * X.imag / I))
    return tuple(np.stack([o[k] for o in out], 1) for k in range(3))

CASES = {'front only':      (['df'],            [3.0],        [0.0],       [60.0]),
         'back only':       (['db'],            [10.0],       [0.0],       [200.0]),
         'front + back':    (['df', 'db'],      [1.0, 10.0],  [0.0, 0.0],  [60.0, 200.0]),
         'back only, f free': (['db', 'f'],     [10.0, 0.46], [0.0, 0.0],  [200.0, 1.0])}

def fit(sheet, case, n_start=16):
    CE = json.load(open(W + 'ce_genosc.json')); GJ = json.load(open(W + 'ce_glass.json'))
    p = np.array(CE['Bphen']['pick']['p'])
    wl, P, D, ang = L.se(sheet, 300.0, 1000.0, 3); meas = L.ncs(P, D)
    A, B = GJ['cauchy']; ng = A + B / (wl / 1000)**2 + 0j; Nf = GO.N_of(p, wl)
    names, x0, lo, hi = CASES[case]; lo, hi = np.array(lo, float), np.array(hi, float)
    def res(x):
        v = dict(zip(names, x))
        N, C_, S_ = ncs_two_sided(wl, Nf, v.get('df', 0.0), v.get('db', 0.0), ng, ang, v.get('f', GJ['f']))
        return np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    rng = np.random.default_rng(zlib.crc32(repr((sheet, case)).encode())); best = None
    for k in range(n_start):
        s = np.array(x0, float) if k == 0 else lo + rng.random(len(lo)) * (hi - lo)
        r = least_squares(res, np.clip(s, lo + 1e-6, hi - 1e-6), bounds=(lo, hi), max_nfev=600)
        if best is None or r.cost < best.cost: best = r
    return dict(sheet=sheet, case=case, x=dict(zip(names, [float(v) for v in best.x])), mse=float(M.mse_ncs(best.fun, len(names))))

def main():
    from multiprocessing import Pool
    jobs = [(sh, c) for sh in ('2_33', '2_34', '2_31') for c in CASES]
    with Pool(4) as pool:
        res = pool.starmap(fit, jobs)
    json.dump(res, open(W + 'bphen10_back.json', 'w'), indent=1)
    for r in res:
        print('%s  %-18s MSE %6.2f | %s' % (r['sheet'], r['case'], r['mse'], '  '.join('%s=%.3g' % kv for kv in r['x'].items())))

if __name__ == '__main__':
    main()

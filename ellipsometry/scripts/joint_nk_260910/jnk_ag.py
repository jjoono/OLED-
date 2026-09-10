"""Stage 2 - Ag on the two seeds, from the SE data alone.

Stack: air / roughness (Bruggeman 50 % Ag + void) / Ag / seed / oxide / Si,
Ag = Einf + Drude(RT) + 3 Gaussians over 260-1080 nm.  Two Gaussians cannot
follow Ag's interband structure at 260-350 nm, and that band carries thickness
information, so it is described rather than discarded.

Two steps, as in the 260819 pass:
  scan  Ag thickness pinned at the deposited value, one common seed thickness
        scanned for all samples of a seed - the Ag-covered pieces need a
        slightly thicker underlayer than the bare-seed fit gives.
  fit   everything free at that seed thickness.
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import os, sys, json, time, numpy as np
from scipy.optimize import least_squares
import jnk_data as D, jnk_ref as R, jnk_model as M, jnk_seed as S

WL_LO, WL_HI = 260.0, 1080.0
D_OX = 2.0
SEEDJ = _os.path.join(ELLIPS_OUT, 'jnk_seed.json')
OUT   = _os.path.join(ELLIPS_OUT, 'jnk_ag.json')
SCAN  = _os.path.join(ELLIPS_OUT, 'jnk_ag_seedscan.json')

#              d_Ag rough Einf  rho    tau    G1A  G1Br G1En   G2A  G2Br G2En   G3A  G3Br G3En
LO = np.array([0.5, 0.0,  0.5,  3e-7,  0.4,   0.0, 0.15, 3.5,  0.0, 0.15, 1.3,  0.0, 0.15, 4.2])
HI = np.array([30., 8.0,  6.0,  1e-1,  50.,   30., 4.00, 6.5,  30., 3.00, 3.5,  30., 4.00, 7.5])
NM = ['d_Ag', 'rough', 'Einf', 'rho', 'tau',
      'G1A', 'G1Br', 'G1En', 'G2A', 'G2Br', 'G2En', 'G3A', 'G3Br', 'G3En']
P0 = np.array([6.0, 1.0, 3.5, 3e-6, 12.0, 3.0, 1.0, 4.2, 0.5, 1.0, 2.5, 3.0, 1.5, 5.5])

HB_EVS, EPS0 = 6.582119569e-16, 8.8541878128e-12


def agN(p, wl):
    return M.osc_N(p[2], [('drude', p[3], p[4]),
                          ('gauss', p[5], p[6], p[7]),
                          ('gauss', p[8], p[9], p[10]),
                          ('gauss', p[11], p[12], p[13])], wl)


def plasma_eV(p):
    """Screened-free-carrier plasma energy from the CompleteEASE Drude pair."""
    return np.sqrt(HB_EVS**2 / (EPS0 * (p[3] / 100.0) * (p[4] * 1e-15)))


def seed_params():
    r = json.load(open(SEEDJ))
    out = {}
    for k, v in r.items():
        rec = [x for x in v['scan'] if abs(x['d_ox'] - D_OX) < 1e-9][0]
        out[k] = (rec['p'][0], np.array(rec['p'][1:]))
    return out


def residual_fun(sheet, d_seed, ps, pin_d=None, step=2):
    wl, P, Dl, _ = D.se(sheet, WL_LO, WL_HI, step)
    meas = D.ncs(P, Dl)
    Ns = S.seed_N(ps, wl)
    ox, si = R.n_sio2(wl), R.n_si(wl)
    amb = R.n_air(wl)
    def r(p):
        d_ag = pin_d if pin_d is not None else p[0]
        Na = agN(p, wl)
        Nr = M.bruggeman(Na, amb, 0.5)
        layers = [amb, Nr, Na, Ns, ox, si]
        ds = [p[1], d_ag, d_seed, D_OX]
        N, C, S_ = M.se_ncs(wl, layers, ds, D.ANG)
        return np.concatenate([(N - meas[0]).ravel(), (C - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    return r, wl


def multistart(r, p0, nstart, rng, pin_first=False):
    lo, hi = LO.copy(), HI.copy()
    if pin_first:
        lo[0], hi[0] = p0[0] - 1e-6, p0[0] + 1e-6
    best = None
    for k in range(nstart):
        if k == 0:
            s = p0
        elif k % 5 == 4:
            s = lo + rng.random(len(lo)) * (hi - lo)
        else:
            s = p0 * (1 + 0.30 * rng.standard_normal(len(p0)))
        s = np.clip(s, lo + 1e-9, hi - 1e-9)
        try:
            b = least_squares(r, s, bounds=(lo, hi), x_scale='jac', max_nfev=1200)
        except Exception:
            continue
        if best is None or b.cost < best.cost:
            best = b
    return best


def _scan_job(a):
    seed, d_seed, sh, nom, ps, nstart = a
    r, _ = residual_fun(sh, d_seed, ps, pin_d=float(nom))
    p0 = P0.copy(); p0[0] = nom
    b = multistart(r, p0, nstart, np.random.default_rng(hash((sh, int(d_seed * 100))) % 2**31),
                   pin_first=True)
    return seed, d_seed, sh, float(M.mse_ncs(b.fun, len(P0) - 1))


def scan_seed(nstart=4, offsets=(-1.0, 0.0, 1.0, 2.0, 3.0), nproc=4):
    from multiprocessing import Pool
    SP = seed_params()
    out = {}
    jobs = []
    for seed in ('HATCN', 'MoOx'):
        d0, ps = SP[seed]
        for off in offsets:
            for smp in [s for s in D.SAMPLES if s[0] == seed and s[2] > 0]:
                jobs.append((seed, d0 + off, smp[3], smp[2], ps, nstart))
    with Pool(nproc) as pool:
        res = pool.map(_scan_job, jobs)
    print('seed-thickness scan, Ag pinned at the deposited thickness')
    for seed in ('HATCN', 'MoOx'):
        d0, _ = SP[seed]
        noms = [s[2] for s in D.SAMPLES if s[0] == seed and s[2] > 0]
        print('%-6s %8s %s' % (seed, 'd_seed', '  '.join('%5d' % n for n in noms)))
        rows = []
        for off in offsets:
            d_seed = d0 + off
            mses = [m for (sd, dd, sh, m) in res if sd == seed and abs(dd - d_seed) < 1e-9]
            rows.append(dict(d_seed=float(d_seed), mse=mses, total=float(np.mean(mses))))
            print('%-6s %8.2f %s   mean %6.2f'
                  % ('', d_seed, '  '.join('%5.2f' % x for x in mses), np.mean(mses)), flush=True)
        best = min(rows, key=lambda x: x['total'])
        out[seed] = dict(d_bare=d0, scan=rows, d_best=best['d_seed'])
        print('  -> %s best d_seed = %.2f nm (bare-seed fit gave %.2f)\n' % (seed, best['d_seed'], d0))
    json.dump(out, open(SCAN, 'w'), indent=1)
    print('saved ->', SCAN)
    return out


def _fit_job(a):
    seed, nom, sh, tag, d_seed, ps, nstart = a
    r, _ = residual_fun(sh, d_seed, ps)
    p0 = P0.copy(); p0[0] = nom
    b = multistart(r, p0, nstart, np.random.default_rng(hash(sh) % 2**31))
    wp = np.array([450., 550., 633., 800., 1000.])
    N = agN(b.x, wp)
    return sh, dict(sheet=sh, seed=seed, ag_nom=nom, tag=tag, d_seed=float(d_seed),
                    mse=float(M.mse_ncs(b.fun, len(P0))), p=[float(v) for v in b.x],
                    n=[float(v) for v in N.real], k=[float(v) for v in N.imag],
                    wp=float(plasma_eV(b.x)), wl_probe=[float(v) for v in wp])


def fit_all(nstart=14, nproc=4):
    from multiprocessing import Pool
    SP = seed_params()
    sc = json.load(open(SCAN)) if os.path.exists(SCAN) else None
    jobs = []
    for seed, ds_nom, nom, sh, tag, lab in D.SAMPLES:
        if nom == 0:
            continue
        d0, ps = SP[seed]
        jobs.append((seed, nom, sh, tag, sc[seed]['d_best'] if sc else d0, ps, nstart))
    with Pool(nproc) as pool:
        res = pool.map(_fit_job, jobs)
    out = dict(res)
    print('%-6s %-6s %4s %6s %7s %7s %7s %8s %8s %7s'
          % ('sheet', 'seed', 'nom', 'MSE', 'd_Ag', 'rough', 'n633', 'k633', 'hw_p eV', 'rho'))
    for sh, v in res:
        x = v['p']
        print('%-6s %-6s %4d %6.2f %7.2f %7.2f %7.3f %8.3f %8.2f %7.2e'
              % (sh, v['seed'], v['ag_nom'], v['mse'], x[0], x[1], v['n'][2], v['k'][2], v['wp'], x[3]))
    json.dump(out, open(OUT, 'w'), indent=1)
    print('saved ->', OUT)
    return out


if __name__ == '__main__':
    import os
    t0 = time.time()
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    if what in ('scan', 'all'):
        scan_seed()
    if what in ('fit', 'all'):
        fit_all()
    print('elapsed %.1f s' % (time.time() - t0))

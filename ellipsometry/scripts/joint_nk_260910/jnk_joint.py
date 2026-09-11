"""Stage 4 - one dielectric function for both measurements.

The Si piece and the glass piece of a given run were coated together, so the Ag
on them should be the same material.  Fit one Einf + Drude + 3 Gaussians to the
Si-piece Psi/Delta AND the glass-piece absolute T and R at once, giving each
piece its own thickness and roughness (they are different pieces, and the
substrate underneath differs).

Three weightings are run for every sample:
  se      T/R ignored, and the glass piece then forced to the Si piece's
          thickness and roughness - i.e. what the SE solution predicts for the
          other substrate if the two pieces really got the same film
  joint   both blocks carry the same total weight
  tr      SE ignored - the T/R data alone, with the SE solution as the start

The three answers together are the result: where they agree the dispersion is
determined by the data, and where they do not, the gap is the systematic error
of a specular-only measurement.  The residual of each block is always reported
against both, so a weighting that buys one block at the other's expense is
visible instead of hidden.
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import os, sys, json, time, numpy as np
from scipy.optimize import least_squares
import jnk_data as D, jnk_ref as R, jnk_model as M, jnk_seed as S, jnk_ag as A
import jnk_tr as TR, jnk_glass as G

SIG_NCS = 0.002        # scatter of the N,C,S residual of the Ag-free fits
SIG_TR  = 0.003        # 0.3 %p, the bare-glass model-vs-measurement gap

AGJ = _os.path.join(ELLIPS_OUT, 'jnk_ag.json')
OUT = _os.path.join(ELLIPS_OUT, 'jnk_joint.json')

# d_Si rough_Si d_gl rough_gl Einf rho tau G1A G1Br G1En G2A G2Br G2En G3A G3Br G3En
LO = np.array([0.5, 0.0, 0.5, 0.0, 0.5, 3e-7, 0.4, 0.0, 0.15, 3.5, 0.0, 0.15, 1.3, 0.0, 0.15, 4.2])
HI = np.array([30., 8.0, 30., 8.0, 6.0, 1e-1, 50., 30., 4.00, 6.5, 30., 3.00, 3.5, 30., 4.00, 7.5])
NM = ['d_Si', 'rg_Si', 'd_gl', 'rg_gl', 'Einf', 'rho', 'tau',
      'G1A', 'G1Br', 'G1En', 'G2A', 'G2Br', 'G2En', 'G3A', 'G3Br', 'G3En']


def agN(p, wl):
    """p[4:] is the same oscillator vector jnk_ag uses from index 2 on."""
    return M.osc_N(p[4], [('drude', p[5], p[6]),
                          ('gauss', p[7], p[8], p[9]),
                          ('gauss', p[10], p[11], p[12]),
                          ('gauss', p[13], p[14], p[15])], wl)


def blocks(sh, rec, ps, step_se=2):
    """Pre-compute everything that does not depend on the parameters."""
    wl_s, P, Dl, _ = D.se(sh, A.WL_LO, A.WL_HI, step_se)
    meas = D.ncs(P, Dl)
    Ns_s = S.seed_N(ps, wl_s)
    ox, si, amb_s = R.n_sio2(wl_s), R.n_si(wl_s), R.n_air(wl_s)

    wl_t, Tm, Rm = TR.measured(rec['tag'])
    Ns_t = S.seed_N(ps, wl_t)
    ng, tau_g = G.n_glass(wl_t), G.tau_glass(wl_t)
    return dict(wl_s=wl_s, meas=meas, Ns_s=Ns_s, ox=ox, si=si, amb_s=amb_s,
                wl_t=wl_t, Tm=Tm, Rm=Rm, Ns_t=Ns_t, ng=ng, tau_g=tau_g,
                d_seed=rec['d_seed'])


def res_se(p, B):
    Na = agN(p, B['wl_s'])
    Nr = M.bruggeman(Na, B['amb_s'], 0.5)
    layers = [B['amb_s'], Nr, Na, B['Ns_s'], B['ox'], B['si']]
    ds = [p[1], p[0], B['d_seed'], A.D_OX]
    N, C, Sx = M.se_ncs(B['wl_s'], layers, ds, D.ANG)
    return np.concatenate([(N - B['meas'][0]).ravel(), (C - B['meas'][1]).ravel(),
                           (Sx - B['meas'][2]).ravel()])


def res_tr(p, B):
    Na = agN(p, B['wl_t'])
    layers, ds = TR.stack_of(Na, p[2], p[3], B['Ns_t'], B['d_seed'], B['wl_t'])
    T, Rr = M.stack_TR(B['wl_t'], layers, ds, B['ng'], tau=B['tau_g'])
    out = []
    if B['Tm'] is not None:
        out.append(T - B['Tm'])
    if B['Rm'] is not None:
        out.append(Rr - B['Rm'])
    return np.concatenate(out) if out else np.zeros(0)


def make_res(B, w_se, w_tr):
    """Blocks scaled by 1/(sigma sqrt(N)) so 'both count the same' is literal."""
    n_se = res_se(np.array(LO), B).size
    n_tr = res_tr(np.array(LO), B).size
    a = w_se / (SIG_NCS * np.sqrt(max(n_se, 1)))
    b = w_tr / (SIG_TR * np.sqrt(max(n_tr, 1)))
    def r(p):
        parts = []
        if w_se > 0:
            parts.append(a * res_se(p, B))
        if w_tr > 0:
            parts.append(b * res_tr(p, B))
        return np.concatenate(parts)
    return r


def score(p, B):
    rs = res_se(p, B)
    rt = res_tr(p, B)
    n_t = B['wl_t'].size
    dT = np.mean(np.abs(rt[:n_t])) if B['Tm'] is not None else np.nan
    dR = np.mean(np.abs(rt[-n_t:])) if B['Rm'] is not None else np.nan
    return M.mse_ncs(rs, 16), 100 * dT, 100 * dR


def _job(a):
    sh, rec, ps, nstart = a
    B = blocks(sh, rec, ps)
    p_se = np.array(rec['p'])
    p0 = np.concatenate([[p_se[0], p_se[1], p_se[0], p_se[1]], p_se[2:]])
    p0 = np.clip(p0, LO + 1e-9, HI - 1e-9)
    rng = np.random.default_rng(hash(sh) % 2**31)
    out = {}
    for tag, (w_se, w_tr) in (('se', (1.0, 0.0)), ('joint', (1.0, 1.0)), ('tr', (0.0, 1.0))):
        r = make_res(B, w_se, w_tr)
        best = None
        for k in range(nstart if tag != 'se' else max(nstart // 2, 2)):
            s = p0 if k == 0 else np.clip(p0 * (1 + 0.25 * rng.standard_normal(len(p0))),
                                          LO + 1e-9, HI - 1e-9)
            try:
                b = least_squares(r, s, bounds=(LO, HI), x_scale='jac', max_nfev=1500)
            except Exception:
                continue
            if best is None or b.cost < best.cost:
                best = b
        x = best.x.copy()
        if tag == 'se':
            # T/R carried no weight, so the glass-side geometry never moved
            # under a gradient; scoring it where the optimiser happened to
            # leave it would be noise.  Pin it to the Si piece instead.
            x[2], x[3] = x[0], x[1]
        mse, dT, dR = score(x, B)
        wp = np.array([450., 550., 633., 700., 800.])
        N = agN(x, wp)
        out[tag] = dict(p=[float(v) for v in x], mse=float(mse), dT=float(dT), dR=float(dR),
                        n=[float(v) for v in N.real], k=[float(v) for v in N.imag],
                        wl_probe=[float(v) for v in wp])
        if tag == 'se':
            p0 = np.clip(x, LO + 1e-9, HI - 1e-9)          # warm-start the others
    out['sheet'] = sh; out['seed'] = rec['seed']; out['ag_nom'] = rec['ag_nom']
    out['tag'] = rec['tag']; out['d_seed'] = rec['d_seed']
    return sh, out


def main(nstart=8, nproc=4):
    from multiprocessing import Pool
    AG = json.load(open(AGJ))
    SP = A.seed_params()
    jobs = [(sh, rec, SP[rec['seed']][1], nstart) for sh, rec in AG.items()
            if rec['tag'] and TR.measured(rec['tag'])[0] is not None]
    with Pool(nproc) as pool:
        res = pool.map(_job, jobs)
    out = dict(res)
    json.dump(out, open(OUT, 'w'), indent=1)
    print('one Ag dielectric function per sample, fitted three ways')
    print('%-6s %-6s %4s | %-22s | %-22s | %-22s'
          % ('sheet', 'seed', 'nom', 'SE only  MSE/dT/dR', 'joint    MSE/dT/dR', 'T/R only MSE/dT/dR'))
    for sh, v in res:
        row = []
        for t in ('se', 'joint', 'tr'):
            row.append('%6.2f %6.2f %6.2f' % (v[t]['mse'], v[t]['dT'], v[t]['dR']))
        print('%-6s %-6s %4d | %s | %s | %s' % (sh, v['seed'], v['ag_nom'], *row))
    print('\nn and k at 550 nm')
    print('%-6s %-6s %4s | %-7s %-7s %-7s | %-7s %-7s %-7s | %-7s %-7s'
          % ('sheet', 'seed', 'nom', 'n SE', 'n join', 'n TR', 'k SE', 'k join', 'k TR',
             'd_Si', 'd_gl'))
    for sh, v in res:
        i = 1        # 550 nm in wl_probe
        print('%-6s %-6s %4d | %-7.3f %-7.3f %-7.3f | %-7.3f %-7.3f %-7.3f | %-7.2f %-7.2f'
              % (sh, v['seed'], v['ag_nom'],
                 v['se']['n'][i], v['joint']['n'][i], v['tr']['n'][i],
                 v['se']['k'][i], v['joint']['k'][i], v['tr']['k'][i],
                 v['joint']['p'][0], v['joint']['p'][2]))
    print('saved ->', OUT)


def rescore(path=None):
    """Recompute every stored block score from the stored parameters.

    Cheap - no fitting - so a change to the scoring convention does not cost a
    refit.  Rewrites the file in place.
    """
    path = path or OUT
    JO = json.load(open(path))
    AG = json.load(open(AGJ))
    SP = A.seed_params()
    for sh, v in JO.items():
        B = blocks(sh, AG[sh], SP[v['seed']][1])
        for tag in ('se', 'joint', 'tr'):
            x = np.array(v[tag]['p'])
            if tag == 'se':
                x[2], x[3] = x[0], x[1]
            mse, dT, dR = score(x, B)
            v[tag].update(p=[float(t) for t in x], mse=float(mse),
                          dT=float(dT), dR=float(dR))
    json.dump(JO, open(path, 'w'), indent=1)
    return JO


if __name__ == '__main__':
    t0 = time.time()
    if len(sys.argv) > 1 and sys.argv[1] == 'rescore':
        JO = rescore()
        print('%-6s %-6s %4s | %-22s | %-22s | %-22s'
              % ('sheet', 'seed', 'nom', 'SE only  MSE/dT/dR', 'joint    MSE/dT/dR',
                 'T/R only MSE/dT/dR'))
        for sh, v in JO.items():
            print('%-6s %-6s %4d | %s' % (sh, v['seed'], v['ag_nom'],
                  ' | '.join('%6.2f %6.2f %6.2f' % (v[t]['mse'], v[t]['dT'], v[t]['dR'])
                             for t in ('se', 'joint', 'tr'))))
    else:
        main()
    print('elapsed %.1f s' % (time.time() - t0))

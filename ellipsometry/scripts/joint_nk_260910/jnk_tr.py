"""Stage 3 - what the absolute T/R on glass says on its own.

Glass piece stack: air / roughness (Bruggeman 50 % Ag + void) / Ag / seed /
1 mm soda-lime (n and single-pass tau from jnk_glass) / air.

Two independent readings of the same films, neither of which uses the Si data:

  d      thickness inverted from T alone and from T+R, with the SE dispersion
         held fixed - the classic cross-check.
  nk     model-free (n, k) per wavelength at the SE thickness.  Two observables,
         two unknowns, so this is an inversion and not a fit: it is exactly the
         quantity that disagreed with SE for the thicker films, because a
         2-observable inversion has no way to tell absorption from light
         scattered out of the specular beam.

Useful window is 420-780 nm: soda-lime starts absorbing below ~400 nm and the
spectrophotometer changes source/detector near 800 nm, where T drops ~0.7 %p
and 1-T-R jumps from 1.0 to 1.5 %.
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import json, numpy as np
from scipy.optimize import least_squares, minimize_scalar
import jnk_data as D, jnk_ref as R, jnk_model as M, jnk_seed as S, jnk_ag as A, jnk_glass as G

WL_LO, WL_HI = 420.0, 780.0
STEP = 2
AGJ = _os.path.join(ELLIPS_OUT, 'jnk_ag.json')
OUT = _os.path.join(ELLIPS_OUT, 'jnk_tr.json')


def measured(tag, step=STEP):
    """T and R on a common grid; either channel may be missing."""
    wT, T = D.tr(tag, 'T', WL_LO, WL_HI)
    wR, Rm = D.tr(tag, 'R', WL_LO, WL_HI)
    wl = wT if wT is not None else wR
    if wl is None:
        return None, None, None
    wl = wl[::step]
    Ti = np.interp(wl, wT, T) if wT is not None else None
    Ri = np.interp(wl, wR, Rm) if wR is not None else None
    return wl, Ti, Ri


def stack_of(Na, d_ag, rough, Ns, d_seed, wl):
    amb = R.n_air(wl)
    Nr = M.bruggeman(Na, amb, 0.5)
    layers, ds = [], []
    if rough > 1e-6:
        layers.append(Nr); ds.append(rough)
    layers += [Na, Ns]; ds += [d_ag, d_seed]
    return layers, ds


def model_TR(Na, d_ag, rough, Ns, d_seed, wl):
    layers, ds = stack_of(Na, d_ag, rough, Ns, d_seed, wl)
    return M.stack_TR(wl, layers, ds, G.n_glass(wl), tau=G.tau_glass(wl))


def invert_d(rec, ps, wl, Tm, Rm):
    """Thickness from T alone and from T+R, SE dispersion and roughness fixed."""
    p = np.array(rec['p'])
    Na = A.agN(p, wl); Ns = S.seed_N(ps, wl)
    def cost(d, use_r):
        T, Rr = model_TR(Na, d, p[1], Ns, rec['d_seed'], wl)
        e = 0.0
        if Tm is not None:
            e = e + np.mean((T - Tm)**2)
        if use_r and Rm is not None:
            e = e + np.mean((Rr - Rm)**2)
        return e
    dT = minimize_scalar(lambda d: cost(d, False), bounds=(1.0, 25.0), method='bounded').x if Tm is not None else np.nan
    dTR = minimize_scalar(lambda d: cost(d, True), bounds=(1.0, 25.0), method='bounded').x
    T, Rr = model_TR(Na, dTR, p[1], Ns, rec['d_seed'], wl)
    resT = np.mean(np.abs(T - Tm)) if Tm is not None else np.nan
    resR = np.mean(np.abs(Rr - Rm)) if Rm is not None else np.nan
    return dT, dTR, resT, resR


def invert_nk(rec, ps, wl, Tm, Rm):
    """Model-free (n, k) per wavelength at the SE geometry."""
    p = np.array(rec['p'])
    Ns = S.seed_N(ps, wl)
    N_se = A.agN(p, wl)
    n_out = np.full(len(wl), np.nan); k_out = np.full(len(wl), np.nan)
    for i, w in enumerate(wl):
        wi = np.array([w])
        def f(v):
            Na = np.array([v[0] + 1j * abs(v[1])])
            T, Rr = model_TR(Na, p[0], p[1], Ns[i:i + 1], rec['d_seed'], wi)
            return [T[0] - Tm[i], Rr[0] - Rm[i]]
        x0 = [max(N_se[i].real, 0.01), max(N_se[i].imag, 0.01)]
        try:
            s = least_squares(f, x0, bounds=([0.001, 0.001], [4.0, 12.0]), xtol=1e-12, ftol=1e-12)
            if np.max(np.abs(s.fun)) < 2e-3:
                n_out[i], k_out[i] = s.x[0], abs(s.x[1])
        except Exception:
            pass
    return n_out, k_out


def main():
    AG = json.load(open(AGJ))
    SP = A.seed_params()
    out = {}
    print('thickness from the glass-side T/R, SE dispersion fixed   (window %g-%g nm)'
          % (WL_LO, WL_HI))
    print('%-6s %-6s %4s %6s %8s %8s %8s | %8s %8s'
          % ('sheet', 'seed', 'nom', 'd_SE', 'd(T)', 'd(T+R)', 'd-d_SE', '|dT| %p', '|dR| %p'))
    for sh, rec in AG.items():
        tag = rec['tag']
        if not tag:
            continue
        wl, Tm, Rm = measured(tag)
        if wl is None or Rm is None:
            continue
        ps = SP[rec['seed']][1]
        dT, dTR, resT, resR = invert_d(rec, ps, wl, Tm, Rm)
        out.setdefault(sh, {}).update(tag=tag, d_se=rec['p'][0], d_T=float(dT), d_TR=float(dTR),
                                      res_T=float(resT), res_R=float(resR))
        print('%-6s %-6s %4d %6.2f %8.2f %8.2f %+8.2f | %8.3f %8.3f'
              % (sh, rec['seed'], rec['ag_nom'], rec['p'][0], dT, dTR, dTR - rec['p'][0],
                 100 * resT, 100 * resR), flush=True)

    print('\nmodel-free (n, k) from T and R at the SE thickness, vs the SE dispersion')
    print('%-6s %-6s %4s | %-15s | %-15s | %-15s'
          % ('sheet', 'seed', 'nom', 'n 550  SE / TR', 'k 550  SE / TR', 'n 700  SE / TR'))
    for sh, rec in AG.items():
        tag = rec['tag']
        if not tag:
            continue
        wl, Tm, Rm = measured(tag)
        if wl is None or Rm is None:
            continue
        ps = SP[rec['seed']][1]
        n_tr, k_tr = invert_nk(rec, ps, wl, Tm, Rm)
        N_se = A.agN(np.array(rec['p']), wl)
        out.setdefault(sh, {}).update(wl=[float(v) for v in wl],
                                      n_tr=[float(v) for v in n_tr],
                                      k_tr=[float(v) for v in k_tr],
                                      n_se=[float(v) for v in N_se.real],
                                      k_se=[float(v) for v in N_se.imag])
        g = lambda a, w: a[np.argmin(abs(wl - w))]
        print('%-6s %-6s %4d | %6.3f / %6.3f | %6.3f / %6.3f | %6.3f / %6.3f'
              % (sh, rec['seed'], rec['ag_nom'],
                 g(N_se.real, 550), g(n_tr, 550), g(N_se.imag, 550), g(k_tr, 550),
                 g(N_se.real, 700), g(n_tr, 700)), flush=True)
    json.dump(out, open(OUT, 'w'), indent=1)
    print('saved ->', OUT)


if __name__ == '__main__':
    main()

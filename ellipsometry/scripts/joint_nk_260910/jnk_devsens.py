"""How far do these constants carry into a device calculation?

The joint fit reproduces both measurements it was given.  That is not the same
as reproducing a device, so this script propagates the remaining ambiguity into
the three quantities a TMM / CPS calculation actually asks for, at a fixed
thickness so only the dispersion varies:

  eps1 = n^2 - k^2   sets where the SPP sits
  eps2 = 2 n k       sets how lossy it is - and n is the poorly determined half
  A_Ag(theta)        electrode absorption inside a device-like stack
  Im[(eps_m - eps_d)/(eps_m + eps_d)]   near-field dipole quenching by the metal
  n_eff, L_spp       the short-range SPP of the embedded film

Read the SE-vs-joint gap as the honest bracket: both describe the ellipsometry
to MSE 2-3, only the joint one also reproduces the absolute T/R.  The T/R-only
column is shown to make the point that intensity alone does not determine n,k -
it is not a candidate answer.
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import json, numpy as np
from scipy.optimize import root
import jnk_joint as J, jnk_model as M

WL = 550.0
K0 = 2 * np.pi / WL
N_ORG, N_SEED, N_CAP, D_CAP = 1.80, 1.85, 2.10, 65.0
SOL = ('se', 'joint', 'tr')


def _load():
    JO = json.load(open(_os.path.join(ELLIPS_OUT, 'jnk_joint.json')))
    rows = sorted((v['ag_nom'], sh) for sh, v in JO.items() if v['seed'] == 'HATCN')
    return JO, rows


def eps_of(JO, sh, sol):
    return J.agN(np.array(JO[sh][sol]['p']), np.array([WL]))[0]**2


def a_ag(JO, sh, sol, ang, d_ag):
    """Ag absorptance in organic / seed / Ag / capping / matched organic.

    Everything but the Ag is lossless, so A = 1 - R - T is exactly the
    electrode loss.  Unpolarised, angle measured inside the organic.
    """
    wl = np.array([WL])
    Na = J.agN(np.array(JO[sh][sol]['p']), wl)
    layers = [np.array([N_ORG + 0j]), np.array([N_SEED + 0j]), Na,
              np.array([N_CAP + 0j]), np.array([N_ORG + 0j])]
    ds = [JO[sh]['d_seed'], d_ag, D_CAP]
    rp, rs = M.tmm(wl, layers, ds, ang)
    tp, ts = M.t_amp(wl, layers, ds, ang)
    # incidence and exit media are the same, so the flux factor is 1
    ap = 1 - abs(rp[0])**2 - abs(tp[0])**2
    a_s = 1 - abs(rs[0])**2 - abs(ts[0])**2
    return 100 * (ap + a_s) / 2


def spp(em, d, e1=N_ORG**2, e2=N_CAP**2):
    """Short-range TM mode of insulator / metal / insulator -> n_eff, L_spp(nm)."""
    def kap(b2, e):
        k = np.sqrt(b2 - e * K0**2 + 0j)
        return k if k.real >= 0 else -k

    def f(x):
        b2 = ((x[0] + 1j * x[1]) * K0)**2
        km, k1, k2 = kap(b2, em), kap(b2, e1), kap(b2, e2)
        r = (np.exp(-2 * km * d)
             - ((em * k1 + e1 * km) / (em * k1 - e1 * km))
             * ((em * k2 + e2 * km) / (em * k2 - e2 * km)))
        return [r.real, r.imag]

    best = None
    for g in (2.5, 3.5, 5.0, 8.0, 12.0, 18.0):
        s = root(f, [g, 0.05], tol=1e-12)
        if s.success and s.x[0] > max(N_ORG, N_CAP) and s.x[1] > 0:
            if best is None or s.x[0] > best[0]:
                best = (s.x[0], s.x[1])
    if best is None:
        return np.nan, np.nan
    return best[0], WL / (4 * np.pi * best[1])


def main():
    JO, rows = _load()
    ed = N_ORG**2
    tri = lambda f: '%8.3f %8.3f %8.3f' % tuple(f(s) for s in SOL)
    gap = lambda f: abs(f('joint') - f('se')) / max(abs(f('joint')), 1e-12) * 100

    print('Ag on HATCN at %g nm.  Thickness held at the joint fit\'s Si-piece value '
          'for all three\ndispersions, so only the optical constants differ.\n' % WL)

    print('%4s %6s | %-26s %7s | %-26s %7s'
          % ('Ag', 'd(nm)', 'eps1   SE / joint / TR', 'SE-join', 'eps2   SE / joint / TR', 'SE-join'))
    for nom, sh in rows:
        d = JO[sh]['joint']['p'][0]
        e1 = lambda s: eps_of(JO, sh, s).real
        e2 = lambda s: eps_of(JO, sh, s).imag
        print('%4d %6.2f | %-26s %6.0f%% | %-26s %6.0f%%'
              % (nom, d, tri(e1), gap(e1), tri(e2), gap(e2)))

    print('\nelectrode absorption, organic(%.1f) / seed / Ag / capping %.0f nm / matched %.1f'
          % (N_ORG, D_CAP, N_ORG))
    for ang in (0.0, 60.0):
        print('  %2.0f deg inside the organic' % ang)
        print('  %4s | %-26s %7s' % ('Ag', 'A(%)   SE / joint / TR', 'SE-join'))
        for nom, sh in rows:
            d = JO[sh]['joint']['p'][0]
            f = lambda s: a_ag(JO, sh, s, ang, d)
            print('  %4d | %-26s %6.0f%%' % (nom, tri(f), gap(f)))

    print('\nnear-field dipole quenching, Im[(eps_m - eps_d)/(eps_m + eps_d)], eps_d = %.2f' % ed)
    print('%4s | %-26s %7s' % ('Ag', 'FOM    SE / joint / TR', 'SE-join'))
    for nom, sh in rows:
        f = lambda s: ((eps_of(JO, sh, s) - ed) / (eps_of(JO, sh, s) + ed)).imag
        print('%4d | %-26s %6.0f%%' % (nom, tri(f), gap(f)))

    print('\nshort-range SPP of the embedded film')
    print('%4s | %-26s | %-26s' % ('Ag', 'n_eff  SE / joint / TR', 'L_spp (nm)  SE / joint / TR'))
    for nom, sh in rows:
        d = JO[sh]['joint']['p'][0]
        r3 = [spp(eps_of(JO, sh, s), d) for s in SOL]
        print('%4d | %8.2f %8.2f %8.2f   | %8.0f %8.0f %8.0f'
              % (nom, *[a for a, _ in r3], *[b for _, b in r3]))

    print('\nNeither data set reaches the evanescent region: SE at 45-65 deg in air '
          'probes u = sin(theta)\nup to 0.91 and the T/R is at normal incidence, '
          'while the CPS integral is dominated by u > 1.\nEverything above u = 1 '
          'is the oscillator model extrapolating.')


if __name__ == '__main__':
    main()

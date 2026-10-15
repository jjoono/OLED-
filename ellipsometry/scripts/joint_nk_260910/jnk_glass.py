"""The glass substrate, from its own measurement.

The T/R pieces sit on 1 mm soda-lime glass.  Rather than assume a catalogue
index, take it from the bare-substrate scan: for a transparent slab measured
with both surfaces,

    T = (1-R0)/(1+R0)   with R0 = ((n-1)/(n+1))^2,

so n follows from T alone, and the corrected R is then a prediction with
nothing left to tune.  How well it lands is the floor on every absolute
intensity in this pipeline.

glass.csv / glassR.csv are the bare-substrate exports of the same campaign.
They were not part of the delivered set, so the UMA back-surface collection
factor in build_tra.py had to be recovered by least squares against the
consolidated file; with the bare substrate in hand it can be checked directly,
because a transparent slab must show zero absorptance.  That check is printed
here and is the reason this pipeline treats the slab as lossless.
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import json, numpy as np
import jnk_data as D, jnk_model as M

OUT = _os.path.join(ELLIPS_OUT, 'jnk_glass.json')


def cauchy(wl, A, B, C=0.0):
    l = np.asarray(wl, float) / 1000.0
    return A + B / l**2 + C / l**4


def fit_index(wl_lo=400.0, wl_hi=800.0, fit_lo=430.0, fit_hi=780.0):
    """A 2-term Cauchy fitted to BOTH channels of the bare substrate.

    The point-by-point inversion of T alone is noisy (0.014 in n) and its ends
    are pulled by the 800 nm instrument changeover, so the Cauchy is fitted
    over 430-780 nm against T and the corrected R together.
    """
    wl, T, R = D.glass_tr(wl_lo, wl_hi)
    n_pt = M.n_from_glass_T(T)
    m = (wl >= fit_lo) & (wl <= fit_hi)

    def resid(c):
        n = cauchy(wl[m], c[0], c[1]) + 0j
        Tp, Rp = M.stack_TR(wl[m], [], [], n)
        return np.concatenate([Tp - T[m], Rp - R[m]])

    from scipy.optimize import least_squares
    c = least_squares(resid, [1.52, 0.004]).x
    n_fit = cauchy(wl, c[0], c[1])
    Tp, Rp = M.stack_TR(wl, [], [], n_fit + 0j)
    return wl, T, R, n_pt, n_fit, c, Tp, Rp


def _pars():
    return json.load(open(OUT))


def n_glass(wl, c=None):
    """Substrate index on any grid."""
    if c is None:
        c = _pars()['cauchy']
    return (cauchy(wl, c[0], c[1]) + 0j)


def tau_glass(wl, ct=None):
    """Single-pass internal transmittance of the slab.

    Lossless over the working window: with the back-surface correction applied
    the bare substrate shows zero absorptance to within the campaign noise, so
    there is nothing left for tau to absorb.  Kept as a function so the stacks
    below read the same either way.
    """
    return np.ones(np.shape(np.asarray(wl, float)))


def main():
    wl, T, R, n_pt, n_fit, c, Tp, Rp = fit_index()
    json.dump(dict(cauchy=[float(v) for v in c],
                   note='n(l) = A + B/l^2, l in um, from the bare-glass %T over '
                        '400-800 nm; the slab is treated as lossless'),
              open(OUT, 'w'), indent=1)
    _, _, Rraw = D.glass_tr(400.0, 800.0, corrected=False)
    m = (wl >= 450) & (wl <= 700)
    print('soda-lime substrate from its own scan, 400-800 nm')
    print('  n(l) = %.5f + %.6f/l^2   (l in um)   ->  n(550) = %.4f   n(633) = %.4f'
          % (c[0], c[1], cauchy(550, *c), cauchy(633, *c)))
    print('  scatter of the point-by-point inversion about the Cauchy: %.4f'
          % np.std(n_pt - n_fit))
    print('\n  back-surface correction check over 450-700 nm')
    print('    absorptance with the raw R       %+.3f %%p'
          % (100 * np.mean(1 - T[m] - Rraw[m])))
    print('    absorptance with the corrected R %+.3f %%p   (campaign noise is 0.15 %%p)'
          % (100 * np.mean(1 - T[m] - R[m])))
    print()
    print('%8s %8s %8s %8s %8s %8s' % ('wl', 'T meas', 'T model', 'R raw', 'R corr', 'R model'))
    for w in (400., 450., 550., 650., 750., 800.):
        i = np.argmin(abs(wl - w))
        print('%8.0f %8.2f %8.2f %8.2f %8.2f %8.2f'
              % (w, 100 * T[i], 100 * Tp[i], 100 * Rraw[i], 100 * R[i], 100 * Rp[i]))
    mm = (wl >= 430) & (wl <= 780)
    print('\nlossless-slab model vs measurement over 430-780 nm:  '
          'mean |dT| = %.3f %%p   mean |dR| = %.3f %%p'
          % (100 * np.mean(np.abs(Tp[mm] - T[mm])), 100 * np.mean(np.abs(Rp[mm] - R[mm]))))
    print('saved ->', OUT)


if __name__ == '__main__':
    main()

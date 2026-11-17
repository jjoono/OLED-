"""Model-free (n, k) per wavelength at the Cauchy thickness, for the thick Si pieces.
No oscillator count assumed: the absorption bands come straight out of the data."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L
import jnk_ref as R, jnk_model as M
from scipy.optimize import least_squares
W = _os.path.join(ELLIPS_OUT, '')
C = json.load(open(W + 'cauchy_si.json'))
D_OX = 2.0

def invert(sh):
    d = C[sh]['d']; A, B, Cc = C[sh]['cauchy']
    wl, P, D, ang = L.se(sh, 245, 1000)
    meas = L.ncs(P, D)
    order = np.argsort(-wl)                 # start in the transparent red, walk into the UV
    n_out = np.full(len(wl), np.nan); k_out = np.full(len(wl), np.nan); res = np.full(len(wl), np.nan)
    x = [A + B / (wl[order[0]] / 1000)**2, 1e-3]
    for i in order:
        w = wl[i:i + 1]
        mi = [m[i:i + 1] for m in meas]
        def f(v):
            N = np.array([v[0] + 1j * v[1]])
            n_, c_, s_ = M.se_ncs(w, [R.n_air(w), N, R.n_sio2(w), R.n_si(w)], [d, D_OX], ang)
            return np.concatenate([(n_ - mi[0]).ravel(), (c_ - mi[1]).ravel(), (s_ - mi[2]).ravel()])
        best = None
        for x0 in (x, [x[0], x[1] + 0.1], [x[0] + 0.15, x[1] + 0.05], [x[0] - 0.15, x[1] + 0.05]):
            r = least_squares(f, x0, bounds=([1.0, 0.0], [3.5, 2.5]))
            # continuity: prefer the root nearest the previous wavelength unless clearly worse
            sc = r.cost + 0.02 * ((r.x[0] - x[0])**2 + (r.x[1] - x[1])**2)
            if best is None or sc < best[0]:
                best = (sc, r)
        r = best[1]
        n_out[i], k_out[i] = r.x; res[i] = np.sqrt(np.mean(r.fun**2)); x = list(r.x)
    return wl, n_out, k_out, res

out = {}
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for sh in ['1_21', '1_22', '1_24', '2_21', '2_22']:
    wl, n, k, res = invert(sh)
    out[sh] = dict(wl=wl.tolist(), n=n.tolist(), k=k.tolist(), res=res.tolist())
    ax[0].plot(wl, n, label='%s (%.1f nm)' % (sh, C[sh]['d'])); ax[1].plot(wl, k, label=sh)
    m = wl < 450
    pk = wl[m][np.argmax(k[m])]
    on = wl[m][np.where(k[m] > 0.02)[0].max()] if (k[m] > 0.02).any() else np.nan
    sec = [(wl[j], k[j]) for j in range(1, len(wl) - 1) if m[j] and k[j] > k[j - 1] and k[j] >= k[j + 1] and k[j] > 0.05]
    print('%-5s d=%5.1f  k_max %.3f at %5.1f nm   onset(k>0.02) %5.1f nm   n550 %.3f   local k maxima: %s   max res %.4f'
          % (sh, C[sh]['d'], np.nanmax(k[m]), pk, on, np.interp(550, wl[::-1], n[::-1]),
             ', '.join('%.0f nm (%.2f)' % (a, b) for a, b in sorted(sec)), np.nanmax(res)))
for a, lab in zip(ax, ('n', 'k')):
    a.set_xlabel('wavelength (nm)'); a.set_ylabel(lab); a.grid(alpha=.3); a.legend(fontsize=8); a.set_xlim(245, 700)
fig.tight_layout(); fig.savefig(W + 'pointnk_si.png', dpi=120)
json.dump(out, open(W + 'pointnk_si.json', 'w'))

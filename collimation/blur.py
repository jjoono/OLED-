"""Lateral spread (blur / pixel crosstalk) of the scatterer + DBR collimator.

Random walk: every traversal of the gap t between device reflector and DBR moves a photon
sideways by t*tan(theta) in a random azimuth, theta drawn from the randomized (Lambertian)
distribution in the medium.  At each DBR encounter the photon escapes with probability T(theta);
otherwise it returns, survives the round trip with probability q, and is re-randomized.
Returns the radii enclosing 50 % and 90 % of the escaped power, in units of t.
"""
import numpy as np
import model as Mo


def spread(T, q, n=400000, seed=0, max_pass=400):
    rng = np.random.default_rng(seed)
    cdf = np.cumsum(Mo.WN); cdf /= cdf[-1]
    x = np.zeros(n); y = np.zeros(n); alive = np.ones(n, bool); out_r = []
    first = True
    for k in range(max_pass):
        idx = np.flatnonzero(alive)
        if idx.size == 0:
            break
        th = np.searchsorted(cdf, rng.random(idx.size))
        th = np.clip(th, 0, Mo.NB - 1)
        phi = rng.random(idx.size) * 2 * np.pi
        step = np.tan(Mo.THR[th])
        if not first:                       # downward + upward traversal of the gap
            th2 = np.clip(np.searchsorted(cdf, rng.random(idx.size)), 0, Mo.NB - 1)
            phi2 = rng.random(idx.size) * 2 * np.pi
            x[idx] += np.tan(Mo.THR[th2]) * np.cos(phi2); y[idx] += np.tan(Mo.THR[th2]) * np.sin(phi2)
        x[idx] += step * np.cos(phi); y[idx] += step * np.sin(phi)
        esc = rng.random(idx.size) < T[th]
        out_r.append(np.hypot(x[idx[esc]], y[idx[esc]]))
        alive[idx[esc]] = False
        rest = idx[~esc]
        alive[rest[rng.random(rest.size) > q]] = False      # absorbed on the round trip
        first = False
    r = np.concatenate(out_r)
    return float(np.percentile(r, 50)), float(np.percentile(r, 90)), r.size / n


if __name__ == '__main__':
    TF = Mo.fresnel_T()
    print('bare (no DBR) Fresnel exit, q=0.97:', spread(TF, 0.97))
    for beta in (10, 20, 30):
        Tid = ((Mo.TH_AIR < beta) & Mo.ESC).astype(float)
        for q in (0.85, 0.95, 0.97):
            print(f'ideal cone {beta} deg, q={q}:', np.round(spread(Tid, q), 3))

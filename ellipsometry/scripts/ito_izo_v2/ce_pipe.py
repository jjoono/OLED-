"""Full pipeline per sample, in CompleteEASE's frame, over the valid range.
  1  geometry scan  (chain residual, 400-1300 nm)
  2  model-free chain n,k at that geometry, 400-1689 nm
  3  Gen-Osc fitted to the chain n,k          -> well-conditioned start
  4  refine all 13 parameters against raw NCS -> final
"""
import numpy as np, json, time
from scipy.optimize import least_squares
import ce_fit as cf, ce_chain as cc, ce_fit2 as f2, ce_osc as osc, ellipsometry_fit as ef

FLO, FHI = 400.0, 1689.0
SHEETS = ['#1', '#2', '#3', '#4', '#5']

def geo_scan(sheet):
    b, _ = __import__('ce_geo').run(sheet, [44., 46., 48., 50., 52.], [0.0, 1.5, 3.0],
                                    [0.2, 0.4, 0.6, 0.8, 1.0], stride=12)
    return b[1], b[2], b[3], b[0]

def stage_osc(wl, n_ch, k_ch, w):
    """fit Einf+Drude+TL+Gaussian to model-free n,k (well conditioned)"""
    def r(p):
        N = f2.film_N(np.r_[0, 0, 0, p], wl)
        return np.concatenate([(N.real - n_ch) * w, (N.imag - k_ch) * w * 3.0])
    lo, hi = f2.LO[3:], f2.HI[3:]
    best = None
    rng = np.random.default_rng(3)
    p0 = np.array([2.0, 6e-4, 6.0, 150., 1.0, 4.1, 3.45, 0.05, 1.0, 1.9])
    for k in range(18):
        s = p0 if k == 0 else np.clip(p0 * (1 + 0.4 * rng.standard_normal(len(p0))), lo + 1e-9, hi - 1e-9)
        if k % 5 == 4: s = lo + rng.random(len(lo)) * (hi - lo)
        try: b = least_squares(r, s, bounds=(lo, hi), x_scale='jac', max_nfev=1500)
        except Exception: continue
        if best is None or b.cost < best.cost: best = b
    return best.x

def run(sheet):
    t0 = time.time()
    d0, rg0, dth0, rms0 = geo_scan(sheet)
    wl, Pm, Dm = cf.load(sheet)
    m = (wl >= FLO) & (wl <= FHI); wlv, Pv, Dv = wl[m], Pm[m], Dm[m]
    ox, si = cf._mat('NTVE_JAW', wlv), cf._mat('SI_JAW', wlv)
    ch = cc.chain(wlv, Pv, Dv, ox, si, d0, rg0, dth0, stride=2)
    good = (ch[:, 0] <= 1250) & (ch[:, 3] < 5 * np.median(ch[:, 3]))
    w = np.where(good, 1.0, 0.15)
    p_osc = stage_osc(ch[:, 0], ch[:, 1], ch[:, 2], w)

    p = np.r_[d0, rg0, dth0, p_osc]
    lo = f2.LO.copy(); hi = f2.HI.copy()
    lo[0], hi[0] = max(30., d0 - 7), d0 + 7
    lo[1], hi[1] = 0.0, 8.0
    lo[2], hi[2] = dth0 - 0.35, dth0 + 0.35
    r = f2.resid_factory(wlv, Pv, Dv, ox, si)
    b = least_squares(r, np.clip(p, lo + 1e-9, hi - 1e-9), bounds=(lo, hi),
                      x_scale='jac', max_nfev=2500)
    M = cf.mse(b.fun, 13)
    at = [f2.NAMES[i] for i, v in enumerate(b.x)
          if abs(v - lo[i]) < 1e-6 * max(1, abs(lo[i])) or abs(v - hi[i]) < 1e-6 * max(1, abs(hi[i]))]
    return dict(sheet=sheet, p=b.x.tolist(), mse=float(M), geo0=[d0, rg0, dth0, rms0],
                at_bound=at, secs=round(time.time() - t0, 1))

if __name__ == '__main__':
    R = {}
    for s in SHEETS:
        r = run(s); R[s] = r
        print('%s  MSE=%7.3f  d=%.2f rough=%.2f dth=%+.3f  scan(d=%.0f,rg=%.1f,dth=%+.1f rms=%.5f)  bound:%s  %.0fs'
              % (s, r['mse'], r['p'][0], r['p'][1], r['p'][2],
                 r['geo0'][0], r['geo0'][1], r['geo0'][2], r['geo0'][3], r['at_bound'] or '-', r['secs']), flush=True)
        print('    ' + '  '.join('%s=%.5g' % (n, v) for n, v in zip(f2.NAMES, r['p'])), flush=True)
        json.dump(R, open('ce_pipe_result.json', 'w'), indent=1)
    print('done')

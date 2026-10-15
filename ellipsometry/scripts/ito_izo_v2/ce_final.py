"""Final Gen-Osc fit over the VALID window only.

  <340 nm  : no homogeneous (n,k) can describe the data (proven by exhaustive
             grid at every roughness 0-15 nm) - excluded
  >1080 nm : Si turns transparent at 1107 nm, wafer BACKSIDE reflections enter
             and the front-surface model breaks - excluded

Seeded from the validated chain->KK solution (genosc_params.json) at the
geometry confirmed against measured transmittance.
"""
import numpy as np, json, math
from scipy.optimize import least_squares
import ce_fit as cf, ce_fit2 as f2

LO_W, HI_W = 340.0, 1080.0
HBAR_EVS, EPS0, HB_FS = 6.582119569e-16, 8.8541878128e-12, 0.6582119569
GEO = {'#1': (50.0, 5.0, 0.38), '#2': (50.5, 3.0, 0.37),
       '#3': (42.0, 2.0, 0.18), '#4': (42.0, 2.0, -0.14), '#5': (42.0, 2.0, 0.75)}
KEY = {'#1': '1', '#2': '2', '#3': '3', '#4': '4', '#5': '5'}

def seed(sh):
    q = json.load(open('genosc_params.json'))[KEY[sh]]['p']
    tau = HB_FS / 0.10
    rho = HBAR_EVS**2 / (q[1]*EPS0*(tau*1e-15)) * 100.0
    d, rg, dth = GEO[sh]
    return np.array([d, rg, dth, q[0], rho, tau, q[2], q[4], q[3], q[5],
                     q[6], q[8]*2*math.sqrt(math.log(2.0)), q[7]])

def run(sh):
    wl, Pm, Dm = cf.load(sh)
    m = (wl >= LO_W) & (wl <= HI_W); wl, Pm, Dm = wl[m], Pm[m], Dm[m]
    ox, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
    r = f2.resid_factory(wl, Pm, Dm, ox, si)
    p0 = seed(sh)
    lo, hi = f2.LO.copy(), f2.HI.copy()
    lo[0], hi[0] = p0[0]-4, p0[0]+4          # geometry already validated vs transmittance
    lo[1], hi[1] = 0.0, 10.0
    lo[2], hi[2] = p0[2]-0.35, p0[2]+0.35
    lo[3], hi[3] = 1.0, 4.5                  # Einf must stay physical
    lo[5], hi[5] = 2.0, 12.0                 # Drude tau
    lo[7], hi[7] = 0.15, 4.0                 # TL Br
    best = None
    rng = np.random.default_rng(7)
    for k in range(10):
        s = p0 if k == 0 else np.clip(p0*(1+0.20*rng.standard_normal(len(p0))), lo+1e-9, hi-1e-9)
        try: b = least_squares(r, np.clip(s, lo+1e-9, hi-1e-9), bounds=(lo, hi),
                               x_scale='jac', max_nfev=1200)
        except Exception: continue
        if best is None or b.cost < best.cost: best = b
    at = [f2.NAMES[i] for i, v in enumerate(best.x)
          if abs(v-lo[i]) < 1e-6*max(1, abs(lo[i])) or abs(v-hi[i]) < 1e-6*max(1, abs(hi[i]))]
    return dict(sheet=sh, p=best.x.tolist(), mse=float(cf.mse(best.fun, 13)),
                mse_seed=float(cf.mse(r(p0), 13)), at_bound=at)

if __name__ == '__main__':
    R = {}
    for sh in ['#1', '#2', '#3', '#4', '#5']:
        x = run(sh); R[sh] = x
        print('%s  MSE %7.3f  (seed %7.3f)   d=%.2f rough=%.2f dth=%+.3f   bound:%s'
              % (sh, x['mse'], x['mse_seed'], x['p'][0], x['p'][1], x['p'][2], x['at_bound'] or '-'), flush=True)
        print('    ' + '  '.join('%s=%.5g' % (n, v) for n, v in zip(f2.NAMES, x['p'])), flush=True)
        json.dump(R, open('ce_final_result.json', 'w'), indent=1)
    print('done')

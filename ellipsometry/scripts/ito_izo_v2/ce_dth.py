"""Is the +0.6 deg angle offset real, or is it standing in for the fixed 3 nm
native oxide / zero roughness?  Angle offset, oxide thickness, roughness and
film thickness all trade against each other, so scan them jointly and look at
the model-free chain residual (no oscillator model involved)."""
import numpy as np, ce_fit as cf, ce_v11 as v11, ellipsometry_fit as ef, tmm_fix
from scipy.optimize import least_squares

sh = '#1'
wl, Pm, Dm = cf.load(sh)
m = (wl >= 340) & (wl <= 1080); wl, Pm, Dm = wl[m], Pm[m], Dm[m]
ox0, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
M = [cf.ncs(Pm[:, i], Dm[:, i]) for i in range(3)]

def chain_res(d, rough, dth, d_ox, stride=14):
    idx = np.arange(len(wl)-1, -1, -stride); g = np.array([2.0, 0.05]); tot = []
    for i in idx:
        w1 = wl[i:i+1]; o1, s1 = ox0[i:i+1], si[i:i+1]
        meas = np.array([[M[a][c][i] for c in range(3)] for a in range(3)]).ravel()
        def res(v):
            Nf = np.array([complex(v[0], max(v[1], 0.0))])
            Nr = ef.bruggeman_ema50(Nf, np.ones_like(Nf)); amb = np.ones_like(Nf); o = []
            for an in cf.ANG0 + dth:
                rp, rs = tmm_fix.tmm(w1, [amb, Nr, Nf, o1, s1], [rough, d, d_ox], an)
                r = np.conj(rp/rs); ps, dl = np.arctan(np.abs(r)), np.angle(r)
                o += [np.cos(2*ps)[0], (np.sin(2*ps)*np.cos(dl))[0], (np.sin(2*ps)*np.sin(dl))[0]]
            return np.array(o) - meas
        b = least_squares(res, g, bounds=([0.05, 0.0], [6.0, 8.0]), xtol=1e-12, ftol=1e-12)
        g = b.x; tot.append(2*b.cost/9)
    return np.sqrt(np.mean(tot))

print('#1: best chain residual, optimising film thickness for each (angle offset, oxide)')
print('    oxide->   1.0nm    2.0nm    3.0nm    4.0nm    5.0nm')
for dth in [0.0, 0.2, 0.4, 0.6, 0.8]:
    row = []
    for dox in [1.0, 2.0, 3.0, 4.0, 5.0]:
        best = min(chain_res(d, 0.0, dth, dox) for d in [53., 54., 55., 56., 57., 58.])
        row.append(best)
    print('  dth%+.1f  ' % dth + ''.join('%8.5f ' % v for v in row))
print()
print('and with roughness instead (oxide fixed 3 nm):')
print('    rough->   0.0nm    1.0nm    2.0nm    3.0nm    5.0nm')
for dth in [0.0, 0.2, 0.4, 0.6, 0.8]:
    row = []
    for rg in [0.0, 1.0, 2.0, 3.0, 5.0]:
        best = min(chain_res(d, rg, dth, 3.0) for d in [52., 54., 55., 56., 57.])
        row.append(best)
    print('  dth%+.1f  ' % dth + ''.join('%8.5f ' % v for v in row))

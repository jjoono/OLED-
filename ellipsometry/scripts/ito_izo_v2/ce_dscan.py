"""Wide thickness / angle-offset scan with the CORRECTED transfer matrix.
The old scan topped out at 54 nm and the fit sat on that edge, and the lab
operator's CompleteEASE fit had claimed 64 nm - so cover the whole range."""
import numpy as np, sys
import ce_fit as cf, ce_v11 as v11

sh = sys.argv[1] if len(sys.argv) > 1 else '#1'
wl, Pm, Dm = cf.load(sh)
m = (wl >= 340) & (wl <= 1080); wl, Pm, Dm = wl[m], Pm[m], Dm[m]
ox, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
res = {}
for dth in [0.0, 0.2, 0.4, 0.6, 0.8]:
    row = []
    for d in range(42, 73, 2):
        c = v11.chain(wl, Pm, Dm, ox, si, float(d), 0.0, dth, stride=16)
        row.append(np.sqrt(np.mean(c[:, 3]**2)))
    res[dth] = row
ds = list(range(42, 73, 2))
print('%s  chain rms (rough=0), corrected TMM' % sh)
print('  d(nm) ' + ''.join('%9.1f' % d for d in ds))
for dth, row in res.items():
    print('  dth%+.1f' % dth + ''.join('%9.5f' % v for v in row))
best = min(((v, dth, ds[i]) for dth, row in res.items() for i, v in enumerate(row)))
print('  best: rms=%.5f at d=%d nm, dth=%+.1f' % (best[0], best[2], best[1]))

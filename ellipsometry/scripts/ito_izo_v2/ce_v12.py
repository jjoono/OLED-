"""Final pipeline, corrected TMM: wide+fine geometry scan, then oscillator fit."""
import numpy as np, json, time
from scipy.optimize import least_squares
import ce_fit as cf, ce_v11 as v11

LOW, HIGH = 340.0, 1080.0

def scan(wl, Pm, Dm, ox, si, ds, rgs, dths, stride):
    best = None
    for d in ds:
        for rg in rgs:
            for dth in dths:
                c = v11.chain(wl, Pm, Dm, ox, si, d, rg, dth, stride)
                s = np.sqrt(np.mean(c[:, 3]**2))
                if best is None or s < best[0]: best = (s, d, rg, dth)
    return best

R = {}
for sh in ['#1', '#2', '#3', '#4', '#5']:
    t0 = time.time()
    wl, Pm, Dm = cf.load(sh)
    m = (wl >= LOW) & (wl <= HIGH); wl, Pm, Dm = wl[m], Pm[m], Dm[m]
    ox, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
    c1 = scan(wl, Pm, Dm, ox, si, [float(d) for d in range(40, 74, 2)], [0.0], [0.0, 0.4, 0.8], 18)
    d0 = c1[1]
    c2 = scan(wl, Pm, Dm, ox, si,
              [d0 + x for x in (-2., -1., -0.5, 0., 0.5, 1., 2.)],
              [0.0, 1.5, 3.0, 5.0], [0.0, 0.2, 0.4, 0.6, 0.8], 12)
    rms, d, rg, dth = c2
    b, _ = v11.fit_osc(sh, wl, Pm, Dm, ox, si, d, dth, nstart=24)
    M = cf.mse(b.fun, len(v11.LO) + 2)
    ch = v11.chain(wl, Pm, Dm, ox, si, d, b.x[0], dth, stride=3)
    at = [i for i, v in enumerate(b.x)
          if abs(v - v11.LO[i]) < 1e-6*max(1, abs(v11.LO[i])) or abs(v - v11.HI[i]) < 1e-6*max(1, abs(v11.HI[i]))]
    NMS = ['rough','Einf','rho','tau','TL_A','TL_Br','TL_Eo','TL_Eg','G_A','G_Br','G_En']
    R[sh] = dict(sheet=sh, d=d, dth=dth, rough_scan=rg, chain_rms=float(rms), mse=float(M),
                 p=[float(v) for v in b.x], at_bound=[NMS[i] for i in at],
                 chain=[[float(v) for v in row] for row in ch])
    print('%s  d=%.1f rough=%.2f dth=%+.2f  chain_rms=%.5f  MSE=%6.2f  bound:%s  %.0fs'
          % (sh, d, b.x[0], dth, rms, M, R[sh]['at_bound'] or '-', time.time()-t0), flush=True)
    print('    Einf=%.3f rho=%.4g tau=%.2f | TL A=%.1f Br=%.3f Eo=%.3f Eg=%.3f | G A=%.4f Br=%.3f En=%.3f'
          % tuple(b.x[1:]), flush=True)
    json.dump(R, open('ce_v12_result.json', 'w'), indent=1)
print('done')

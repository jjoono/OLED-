"""Model-free per-wavelength (n,k) extraction over the FULL 192-1688 nm range,
in CompleteEASE's frame, with a geometry scan over (d, roughness, angle offset)."""
import numpy as np, ce_fit as cf, ellipsometry_fit as ef
from scipy.optimize import least_squares

D_OX = cf.D_OX
ANG0 = cf.ANG0

def chain(wl, Pm, Dm, N_ox, N_si, d, rough, dth, stride=1, n0=2.0, k0=0.05):
    idx = np.arange(len(wl)-1, -1, -stride)          # long -> short wavelength
    M = [cf.ncs(Pm[:, i], Dm[:, i]) for i in range(3)]
    ang = ANG0 + dth
    out = np.empty((len(idx), 4))
    g = np.array([n0, k0])
    for t, i in enumerate(idx):
        w1 = wl[i:i+1]
        ox, si = N_ox[i:i+1], N_si[i:i+1]
        meas = np.array([[M[a][c][i] for c in range(3)] for a in range(3)]).ravel()
        def res(v):
            Nf = np.array([complex(v[0], max(v[1], 0.0))])
            Nr = ef.bruggeman_ema50(Nf, np.ones_like(Nf))
            amb = np.ones_like(Nf)
            o = []
            for a in ang:
                rp, rs = ef._tmm(w1, [amb, Nr, Nf, ox, si], [rough, d, D_OX], a)
                r = rp/rs
                p, dl = np.arctan(np.abs(r)), np.angle(r)
                o += [np.cos(2*p)[0], (np.sin(2*p)*np.cos(dl))[0], (np.sin(2*p)*np.sin(dl))[0]]
            return np.array(o) - meas
        b = least_squares(res, g, bounds=([0.05, 0.0], [6.0, 8.0]), xtol=1e-12, ftol=1e-12)
        g = b.x
        out[t] = (wl[i], b.x[0], b.x[1], np.sqrt(2*b.cost/9))
    return out[::-1]

def score(wl, Pm, Dm, N_ox, N_si, d, rough, dth, stride=8):
    c = chain(wl, Pm, Dm, N_ox, N_si, d, rough, dth, stride)
    return np.sqrt(np.mean(c[:, 3]**2)), c

if __name__ == '__main__':
    import sys, json
    sheet = sys.argv[1] if len(sys.argv) > 1 else '#1'
    wl, Pm, Dm = cf.load(sheet)
    N_ox, N_si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
    best = None
    print('coarse scan (stride 12)')
    for d in [46., 49., 51., 53., 56.]:
        for rg in [0.0, 1.5, 3.0]:
            for dth in [0.0, 0.4, 0.8]:
                s, _ = score(wl, Pm, Dm, N_ox, N_si, d, rg, dth, stride=12)
                if best is None or s < best[0]:
                    best = (s, d, rg, dth); print('   d=%.1f rg=%.1f dth=%+.2f  rms=%.5f  *' % (d, rg, dth, s))
    print('coarse best: rms=%.5f  d=%.1f rough=%.1f dth=%+.2f' % best)
    json.dump(dict(sheet=sheet, rms=best[0], d=best[1], rough=best[2], dth=best[3]),
              open('ce_geo_%s.json' % sheet.strip('#'), 'w'))

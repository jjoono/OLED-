"""Chain extraction with a coarse (n,k) GRID initialiser at every wavelength,
so no wavelength can be trapped in the k=0 local minimum, then a local refine.
The per-wavelength residual doubles as a map of where a homogeneous single-layer
model is capable of describing the data at all."""
import numpy as np, ce_fit as cf, ellipsometry_fit as ef
from scipy.optimize import least_squares

NG = np.linspace(0.10, 4.20, 66)
KG = np.linspace(0.00, 3.00, 51)

def extract(wl, Pm, Dm, N_ox, N_si, d, rough, dth, stride=1):
    M = [cf.ncs(Pm[:, i], Dm[:, i]) for i in range(3)]
    ang = cf.ANG0 + dth
    NN, KK = np.meshgrid(NG, KG, indexing='ij')
    Gf = (NN + 1j*KK).ravel()
    Gr = ef.bruggeman_ema50(Gf, np.ones_like(Gf))
    amb = np.ones_like(Gf)
    idx = np.arange(0, len(wl), stride)
    out = np.empty((len(idx), 4))
    for t, i in enumerate(idx):
        meas = np.array([[M[a][c][i] for c in range(3)] for a in range(3)]).ravel()
        # --- grid ---
        w = np.full(Gf.shape, wl[i]); oxv = np.full(Gf.shape, N_ox[i]); siv = np.full(Gf.shape, N_si[i])
        tot = np.zeros(Gf.shape)
        for a, an in enumerate(ang):
            rp, rs = ef._tmm(w, [amb, Gr, Gf, oxv, siv], [rough, d, cf.D_OX], an)
            r = rp/rs; ps, dl = np.arctan(np.abs(r)), np.angle(r)
            for c, v in enumerate([np.cos(2*ps), np.sin(2*ps)*np.cos(dl), np.sin(2*ps)*np.sin(dl)]):
                tot += (v - meas[3*a+c])**2
        g = np.array([NN.ravel()[np.argmin(tot)], KK.ravel()[np.argmin(tot)]])
        # --- refine ---
        w1 = wl[i:i+1]; ox1, si1 = N_ox[i:i+1], N_si[i:i+1]
        def res(v):
            Nf = np.array([complex(v[0], max(v[1], 0.0))])
            Nr = ef.bruggeman_ema50(Nf, np.ones_like(Nf)); a1 = np.ones_like(Nf)
            o = []
            for an in ang:
                rp, rs = ef._tmm(w1, [a1, Nr, Nf, ox1, si1], [rough, d, cf.D_OX], an)
                r = rp/rs; ps, dl = np.arctan(np.abs(r)), np.angle(r)
                o += [np.cos(2*ps)[0], (np.sin(2*ps)*np.cos(dl))[0], (np.sin(2*ps)*np.sin(dl))[0]]
            return np.array(o) - meas
        b = least_squares(res, g, bounds=([0.05, 0.0], [6.0, 8.0]), xtol=1e-13, ftol=1e-13)
        out[t] = (wl[i], b.x[0], b.x[1], np.sqrt(2*b.cost/9))
    return out

if __name__ == '__main__':
    wl, Pm, Dm = cf.load('#1')
    ox, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
    c = extract(wl, Pm, Dm, ox, si, 51.0, 3.0, 0.0, stride=4)
    np.save('ce_chain1.npy', c)
    print(' wl        n        k      resid')
    for r in c[::6]:
        print('%7.1f %8.3f %8.4f  %8.5f' % tuple(r))
    print()
    for lo, hi in [(192,250),(250,300),(300,350),(350,400),(400,500),(500,900),
                   (900,1200),(1200,1400),(1400,1550),(1550,1689)]:
        m = (c[:,0] >= lo) & (c[:,0] < hi)
        if m.sum(): print('  rms %4d-%4d = %.5f  (%d pts)' % (lo, hi, np.sqrt(np.mean(c[m,3]**2)), m.sum()))

"""Exhaustive (n,k) grid at selected wavelengths: is the UV/NIR misfit a local
minimum, or genuinely unreachable by ANY homogeneous (n,k)?"""
import numpy as np, ce_fit as cf, ellipsometry_fit as ef

wl, Pm, Dm = cf.load('#1')
ox, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
M = [cf.ncs(Pm[:, i], Dm[:, i]) for i in range(3)]

def best_nk(i, d, rough, dth, ngrid=260, kgrid=200):
    w1 = wl[i:i+1]
    meas = np.array([[M[a][c][i] for c in range(3)] for a in range(3)]).ravel()
    ns = np.linspace(0.10, 4.5, ngrid); ks = np.linspace(0.0, 4.0, kgrid)
    NN, KK = np.meshgrid(ns, ks, indexing='ij')
    Nf = (NN + 1j*KK).ravel()
    Nr = ef.bruggeman_ema50(Nf, np.ones_like(Nf))
    amb = np.ones_like(Nf); w = np.full(Nf.shape, wl[i])
    oxv = np.full(Nf.shape, ox[i]); siv = np.full(Nf.shape, si[i])
    tot = 0.0
    for a, ang in enumerate(cf.ANG0 + dth):
        rp, rs = ef._tmm(w, [amb, Nr, Nf, oxv, siv], [rough, d, cf.D_OX], ang)
        r = rp/rs; p, dl = np.arctan(np.abs(r)), np.angle(r)
        for c, v in enumerate([np.cos(2*p), np.sin(2*p)*np.cos(dl), np.sin(2*p)*np.sin(dl)]):
            tot = tot + (v - meas[3*a+c])**2
    j = np.argmin(tot)
    return NN.ravel()[j], KK.ravel()[j], np.sqrt(tot[j]/9)

print('exhaustive grid, d=51 rough=3 dth=0:')
print('   wl        n       k     resid')
for t in [200, 250, 300, 350, 400, 600, 1000, 1400, 1650]:
    i = np.argmin(abs(wl-t))
    n, k, r = best_nk(i, 51.0, 3.0, 0.0)
    print('  %6.1f  %7.3f %7.3f  %8.5f' % (wl[i], n, k, r))

print('\nsame, but d=46 rough=1.5 dth=0.8:')
for t in [200, 250, 300, 350, 400, 600, 1000, 1400, 1650]:
    i = np.argmin(abs(wl-t))
    n, k, r = best_nk(i, 46.0, 1.5, 0.8)
    print('  %6.1f  %7.3f %7.3f  %8.5f' % (wl[i], n, k, r))

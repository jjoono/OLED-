import numpy as np, tmodel as T, hi175 as H
from multiprocessing import Pool
NK = {8: (H.load(8), 7.87), 12: (H.load(12), 10.71)}
LF = np.arange(560, 721, 4)
EFF = H.eff(tuple(LF))
def fine(a):
    ag, z, dE, dH = a; nk, d = NK[ag]
    w = np.array([T.S[T.idx(l)] for l in LF]); w /= w.sum()
    eo = np.array([T.eta_out(l, dE, dH, 10, d, *EFF[l], agtop=nk, tfe='al2o3', al=0.2, zns=z)[0] for l in LF])
    return a, T.ETA_RAD * float(w @ eo), H.eqe(dE, dH, d, nk, dict(al=0.2, zns=z)), eo
if __name__ == '__main__':
    J = [(8, 0, 131.0, 78.1), (8, 100, 177.5, 222.3), (8, 120, 191.6, 215.2), (8, 130, 200.0, 281.8), (12, 100, 191.9, 221.6), (12, 0, 131.9, 91.6)]
    with Pool(4) as p:
        for a, ef, e6, eo in p.map(fine, J):
            print(a, 'fine(41 lam) %.3f  6-lam %.3f' % (ef, e6), ' eta_out(lam) min/max %.2f/%.2f' % (eo.min(), eo.max()))

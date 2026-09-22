import numpy as np, scipy.io as sio, sys, csv
from multiprocessing import Pool
sys.path.insert(0, '/home/user/OLED-/sim/fig3c'); sys.path.insert(0, '/home/user/OLED-/sim/mla')
import cps2, materials as M, series
AL = 0.958336812040867 + 6.68678039868499j
TH = np.arange(90) + 0.5
W = np.cos(np.radians(TH)) * np.sin(np.radians(TH))
B = sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat')['BSDF_MLA']
NS = np.round(np.arange(1.30, 2.001, 0.05), 2)
def bsdf(n):
    k = int(np.argmin(abs(NS - n))); b = B[:, :, k]
    return b[:90, :].sum(0), b[90:180, :][::-1, :]
def stack(refl, n_sub, lam):
    return cps2.Stack(lam, (1.8, 1.8), 20.0, 10.0, above=[(1.8, 1.8, 200.0), (refl, refl, 100.0)],
                      below=[(1.8, 1.8, 200.0), (M.ITO, M.ITO, 50.0)], n_sub=n_sub, h=2/3)
LAM = np.arange(490.0, 610.1, 5.0)
WG = np.exp(-4 * np.log(2) * ((LAM - 550.0) / 70.0) ** 2); WG /= WG.sum()      # Gaussian, FWHM 70 nm
WU = np.ones_like(LAM) / LAM.size                                              # uniform 490-610
def one(args):
    name, refl, n = args
    BT, BR = bsdf(n)
    es, po, pp, A = [], [], [], []
    for lam in LAM:
        S = stack(refl, n, lam); r = cps2.solve_pol(S, npts=12000)
        R = cps2.stack_reflectance(S, TH); P = cps2.sub_angular(S, TH)
        e, _ = series.eta_ext(BT, BR, R, P, n_term=100)
        esub = r['air'] + r['sub']
        es.append(esub); po.append(esub * e); pp.append(BT @ (P / P.sum())); A.append(1 - np.sum(R * W) / np.sum(W))
    es, po, pp, A = map(np.array, (es, po, pp, A))
    out = [name, '%.2f' % n]
    i0 = int(np.argmin(abs(LAM - 550)))
    for w in (None, WG, WU):
        if w is None: out += ['%.6f' % es[i0], '%.6f' % (po[i0] / es[i0]), '%.6f' % po[i0], '%.6f' % pp[i0], '%.6f' % A[i0]]
        else: out += ['%.6f' % (w @ es), '%.6f' % ((w @ po) / (w @ es)), '%.6f' % (w @ po), '%.6f' % ((w @ (pp * es)) / (w @ es)), '%.6f' % (w @ A)]
    print(out, flush=True); return out
if __name__ == '__main__':
    jobs = [(nm, rf, float(n)) for nm, rf in (('Al', AL), ('Ag', M.AG)) for n in NS]
    with Pool(4) as p: rows = p.map(one, jobs)
    hdr = ['reflector', 'n_sub'] + [f'{k}_{t}' for t in ('mono550', 'gauss70', 'uniform') for k in ('eta_sub', 'eta_ext', 'EQE', 'p_Psub', 'Aprime')]
    with open('fig2def_specavg.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(hdr); w.writerows(rows)
    print('done')

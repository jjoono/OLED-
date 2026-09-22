import numpy as np, scipy.io as sio, sys, csv
sys.path.insert(0, '/home/user/OLED-/sim/fig3c'); sys.path.insert(0, '/home/user/OLED-/sim/mla')
import cps2, materials as M, series
AL = 0.958336812040867 + 6.68678039868499j
TH = np.arange(90) + 0.5; W = np.cos(np.radians(TH)) * np.sin(np.radians(TH))
B = sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat')['BSDF_MLA']
NS = np.round(np.arange(1.30, 2.001, 0.05), 2)
def bsdf(n):                       # linear interpolation between the LightTools slices
    j = min(int((n - 1.30) / 0.05 + 1e-9), 13); f = (n - NS[j]) / 0.05
    b = (1 - f) * B[:, :, j] + f * B[:, :, j + 1]
    return b[:90, :].sum(0), b[90:180, :][::-1, :]
def stack(refl, n_sub):
    return cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0, above=[(1.8, 1.8, 200.0), (refl, refl, 100.0)],
                      below=[(1.8, 1.8, 200.0), (M.ITO, M.ITO, 50.0)], n_sub=n_sub, h=2/3)
rows = []
for nm, refl in (('Al', AL), ('Ag', M.AG)):
    for n in np.round(np.arange(1.30, 2.0001, 0.01), 2):
        S = stack(refl, float(n)); BT, BR = bsdf(float(n)); R = cps2.stack_reflectance(S, TH); P = cps2.sub_angular(S, TH)
        r = cps2.solve_pol(S, npts=12000); es = r['air'] + r['sub']; A = 1 - np.sum(R * W) / np.sum(W)
        pc = np.sum(BT * W) / np.sum(W); e, _ = series.eta_ext(BT, BR, R, P, 100); ecf = pc / (pc + (1 - pc) * A)
        P70 = P[TH > 70].sum() / P.sum()
        rows.append([nm, '%.2f' % n, '%.6f' % es, '%.6f' % r['wg'], '%.6f' % A, '%.6f' % pc, '%.6f' % (BT @ (P / P.sum())), '%.6f' % e, '%.6f' % ecf, '%.6f' % (es * e), '%.6f' % (es * ecf), '%.6f' % P70])
with open('fig2def_fine.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['reflector', 'n_sub', 'eta_sub', 'wg', 'Aprime', 'p_cos_sin', 'p_Psub', 'eta_ext_series', 'eta_ext_closed', 'EQE_series', 'EQE_closed', 'graze70']); w.writerows(rows)
print('done')

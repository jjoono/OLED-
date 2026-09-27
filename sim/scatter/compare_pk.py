import numpy as np, scipy.io as sio, sys, csv
sys.path[:0] = ['/home/user/OLED-/sim/fig3c']
import cps2, materials as M
TH = np.arange(90) + 0.5; g = np.cos(np.radians(TH))*np.sin(np.radians(TH)); g /= g.sum()
def split(B): return B[:90].sum(0), B[90:][::-1]
mla = sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat')['BSDF_MLA']
S = {'MLA (hemisphere)': mla[:, :, 10]}
S.update({f'scattering {k}': v for k, v in np.load('scatter_bsdf_n1.8.npz').items()})
EML = (1.83273, 1.67122); TAPC = (1.69075, 1.66375); B3 = M.ETL['B3PyMPM']
St = cps2.Stack(550.0, EML, 25.0, 12.5, above=[(B3[0], B3[1], 300.0), (M.AG, M.AG, 100.0)],
                below=[(TAPC[0], TAPC[1], 300.0), (M.ITO, M.ITO, 50.0)], n_sub=1.8)
R = np.real(cps2.stack_reflectance(St, TH)); P = np.real(cps2.sub_angular(St, TH)); P /= P.sum()
A_l = 1 - R @ g
rows = [['structure', 'input', 'k', 'p_k', 'p_k/p_Lamb', 'frac>60 arriving']]
print(f"{'structure':26s} p_Lamb | Lambertian input, lossless: p_1..p_6 | device (Ag stack n1.8, 550 nm): p_1..p_6 | eta_ext series / eq3")
for name, B in S.items():
    BT, BR = split(B); pl = BT @ g; out = {}
    for inp, v0, RR in (('Lambertian, R=1', g, np.ones(90)), ('device', P, R)):
        v = v0.copy(); ps = []; e = 0
        for k in range(300):
            U = v.sum(); pk = BT @ v / U; e += BT @ v
            if k < 30: rows.append([name, inp, k + 1, pk, pk/pl, v[TH > 60].sum()/U]); ps.append(pk)
            v = RR * (BR @ v)
        out[inp] = (ps, e)
    e3 = pl/(pl + (1 - pl)*A_l)
    print(f"{name:26s} {pl:.3f} | " + ' '.join(f'{x:.3f}' for x in out['Lambertian, R=1'][0][:6]) + ' | ' + ' '.join(f'{x:.3f}' for x in out['device'][0][:6]) + f" | {out['device'][1]:.4f} / {e3:.4f}")
with open('pk_mla_vs_scattering.csv', 'w', newline='') as f: csv.writer(f).writerows(rows)

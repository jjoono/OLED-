import numpy as np, sys, csv
from multiprocessing import Pool
from scatter_bsdf import bsdf
sys.path[:0] = ['/home/user/OLED-/sim/fig3c']
import cps2, materials as M
TH = np.arange(90) + 0.5; g0 = np.cos(np.radians(TH))*np.sin(np.radians(TH)); g0 /= g0.sum()
EML = (1.83273, 1.67122); TAPC = (1.69075, 1.66375); B3 = M.ETL['B3PyMPM']
St = cps2.Stack(550.0, EML, 25.0, 12.5, above=[(B3[0], B3[1], 300.0), (M.AG, M.AG, 100.0)],
                below=[(TAPC[0], TAPC[1], 300.0), (M.ITO, M.ITO, 50.0)], n_sub=1.8)
R = np.real(cps2.stack_reflectance(St, TH)); P = np.real(cps2.sub_angular(St, TH)); P /= P.sum()
def ext(B):
    BT, BR = B[:90].sum(0), B[90:][::-1]; v = P.copy(); e = 0
    for k in range(400): e += BT @ v; v = R*(BR @ v)
    pl = BT @ g0; A = 1 - R @ g0; return e, pl, pl/(pl + (1 - pl)*A)
def job(a):
    S, g = a; return (S, g) + ext(bsdf(1.8, S, g, N=6000, seed=int(S*100 + g*10)))
if __name__ == '__main__':
    jobs = [(S, g) for S in (0.5, 1, 2, 3, 5, 7, 10, 15, 20, 30) for g in (0.0, 0.5, 0.7, 0.8, 0.9, 0.95)]
    with Pool(4) as p: res = p.map(job, jobs)
    with open('sweep_sg_n1.8.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['S', 'g', 'eta_ext_series', 'p_Lamb', 'eta_ext_eq3']); w.writerows(res)
    import scipy.io as sio
    mla = sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat')['BSDF_MLA'][:, :, 10]
    print('MLA', [round(x, 4) for x in ext(mla)])
    for r in sorted(res, key=lambda r: -r[2])[:8]: print([round(x, 4) for x in r])
    print('worst', [round(x, 4) for x in min(res, key=lambda r: r[2])])

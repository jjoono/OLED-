"""Table 1 over the full emission spectrum (400-700 nm, Ir(ppy)2acac, photon-number weighting) and at 550 nm,
with the dispersive optical constants of the author's library (nk_JH_total.mat) in the Python solver."""
import numpy as np, scipy.io as sio, sys, csv
from multiprocessing import Pool
sys.path.insert(0, '/home/user/OLED-/sim/fig3c'); sys.path.insert(0, '/home/user/OLED-/sim/mla')
import cps2, series
m = sio.loadmat('/home/user/OLED-/sim/table1_author/nk_JH_total.mat', squeeze_me=True, struct_as_record=False)
mat, spec = m['material'], m['spectrum']
LAM = np.arange(400, 701)                                   # library grid 400-800 (401), spectrum grid 400-700 (301)
def lib(name): return np.asarray(getattr(mat, name))[:301]
AG, ITO = lib('l_Ag_McPeak'), lib('l_ITO')
B3o, B3e, EMo, EMe, TPo, TPe = lib('l_B3_o_JO'), lib('l_B3_e_JO'), lib('l_TCTA_B3_o_JO'), lib('l_TCTA_B3_e_JO'), lib('l_TAPC_o_JO'), lib('l_TAPC_e_JO')
I = np.asarray(spec.I_Irppy2acac)[:301]; W = LAM * I / np.sum(LAM * I)   # photon-number weight
B = sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat')['BSDF_MLA'][:, :, 10]; BT = B[:90].sum(0); BR = B[90:180][::-1]
TH = np.arange(90) + 0.5
k550 = ITO[150].imag
def one(args):
    d_etl, k_ito, h, j = args
    lam = float(LAM[j]); ito = complex(ITO[j].real, ITO[j].imag * k_ito / k550)
    S = cps2.Stack(lam, (EMo[j], EMe[j]), 25.0, 12.5, above=[(B3o[j], B3e[j], d_etl), (AG[j], AG[j], 100.0)],
                   below=[(TPo[j], TPe[j], 180.0), (ito, ito, 50.0)], n_sub=1.8, h=h)
    r = cps2.solve_pol(S, npts=8000); R = cps2.stack_reflectance(S, TH); P = cps2.sub_angular(S, TH)
    e, _ = series.eta_ext(BT, BR, R, P, 100); es = r['air'] + r['sub']
    return (d_etl, k_ito, h, j, es, r['spp'], e, es * e)
if __name__ == '__main__':
    cases = [(d, k, h) for d in (200.0, 300.0) for k in (0.0032, 0.002) for h in (2/3, 0.8, 0.9)]
    jobs = [(d, k, h, j) for (d, k, h) in cases for j in range(301)]
    with Pool(4) as p: res = p.map(one, jobs, chunksize=20)
    out = [['d_ETL_nm', 'k_ITO', 'Theta', 'mode', 'eta_sub', 'spp', 'eta_ext', 'EQE']]
    for (d, k, h) in cases:
        rows = sorted([r for r in res if r[0] == d and r[1] == k and r[2] == h], key=lambda r: r[3])
        es = np.array([r[4] for r in rows]); sp = np.array([r[5] for r in rows]); q = np.array([r[7] for r in rows])
        out.append([int(d), k, round(h, 3), 'full_400_700', W @ es, W @ sp, (W @ q) / (W @ es), W @ q])
        r5 = rows[150]; out.append([int(d), k, round(h, 3), 'single_550', r5[4], r5[5], r5[6], r5[7]])
        print(out[-2]); print(out[-1], flush=True)
    with open('table1_full_spectrum.csv', 'w', newline='') as f: csv.writer(f).writerows(out)
    print('done')

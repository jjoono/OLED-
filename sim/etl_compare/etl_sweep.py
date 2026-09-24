"""ETL comparison at 550 nm in the Table 1 stack: anisotropic B3PyMPM vs isotropic low-index 3TPYMB,
plus control ETLs (isotropic n_o, isotropic n_e) to separate the effect of n_o and n_e."""
import numpy as np, scipy.io as sio, sys, csv
from multiprocessing import Pool
sys.path.insert(0, '/home/user/OLED-/sim/fig3c'); sys.path.insert(0, '/home/user/OLED-/sim/mla')
import cps2, series
m = sio.loadmat('/home/user/OLED-/sim/table1_author/nk_JH_total.mat', squeeze_me=True, struct_as_record=False)['material']
g = lambda n: complex(np.asarray(getattr(m, n))[150])
AG, ITO = g('l_Ag_McPeak'), g('l_ITO')
EMo, EMe, TPo, TPe = g('l_TCTA_B3_o_JO'), g('l_TCTA_B3_e_JO'), g('l_TAPC_o_JO'), g('l_TAPC_e_JO')
B3o, B3e, TY = g('l_B3_o_JO'), g('l_B3_e_JO'), g('TPYMB_3')
ETLS = {'B3PyMPM': (B3o, B3e), '3TPYMB': (TY, TY), 'iso_no': (B3o, B3o), 'iso_ne': (B3e, B3e), 'B3_ne_1.5': (B3o, 1.5+0j)}
B = sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat')['BSDF_MLA'][:, :, 10]; BT = B[:90].sum(0); BR = B[90:180][::-1]
TH = np.arange(90) + 0.5
def one(a):
    name, d, h = a; no, ne = ETLS[name]
    S = cps2.Stack(550.0, (EMo, EMe), 25.0, 12.5, above=[(no, ne, d), (AG, AG, 100.0)],
                   below=[(TPo, TPe, 180.0), (ITO, ITO, 50.0)], n_sub=1.8, h=h)
    r = cps2.solve_pol(S, npts=8000); R = cps2.stack_reflectance(S, TH); P = cps2.sub_angular(S, TH)
    e, _ = series.eta_ext(BT, BR, R, P, 100); es = r['air'] + r['sub']
    Ap = 1 - (P / P.sum()) @ R
    return [name, d, round(h, 3), es, r['spp'], r['wg'], r['abs'], Ap, e, es * e]
if __name__ == '__main__':
    jobs = [(n, float(d), h) for n in ETLS for d in range(50, 425, 25) for h in (2/3, 0.8)]
    with Pool(4) as p: res = p.map(one, jobs)
    with open('/home/user/OLED-/sim/etl_compare/etl_sweep.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['ETL', 'd_ETL', 'Theta', 'eta_sub', 'spp', 'wg', 'abs', 'Aprime_Psub', 'eta_ext', 'EQE']); w.writerows(res)

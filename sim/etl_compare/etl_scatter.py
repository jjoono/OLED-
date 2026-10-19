"""Same ETL comparison as etl_sweep.py (550 nm, Ag 100 | ETL d | TCTA:B3PyMPM 25 | TAPC 180 | ITO 50 | n 1.8 substrate),
with the hemisphere MLA replaced by index-matched (n 1.8) volume scattering layers on the substrate
(Monte-Carlo BSDFs sim/scatter/scatter_bsdf_n1.8.npz: optical thickness S 1/3/10, HG g 0.5/0.9, planar top to air).
Also a flat substrate (no extractor) for reference. -> etl_scatter.csv"""
import numpy as np, scipy.io as sio, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla', '/home/user/OLED-/sim/etl_compare']
import cps2, series
from etl_sweep import ETLS, AG, ITO, EMo, EMe, TPo, TPe, TH

def split(B): return B[:90].sum(0), B[90:180][::-1]
EXT = {'MLA': split(sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat')['BSDF_MLA'][:, :, 10])}
EXT.update({k: split(v) for k, v in np.load('/home/user/OLED-/sim/scatter/scatter_bsdf_n1.8.npz').items()})
# flat substrate/air: specular, unpolarised Fresnel
c = np.cos(np.radians(TH)); s2 = (1.8 * np.sin(np.radians(TH))) ** 2; ok = s2 < 1; ct = np.sqrt(np.clip(1 - s2, 0, 1))
rs = ((1.8 * c - ct) / (1.8 * c + ct)) ** 2; rp = ((1.8 * ct - c) / (1.8 * ct + c)) ** 2
Rf = np.where(ok, 0.5 * (rs + rp), 1.0); EXT['flat'] = (1 - Rf, np.diag(Rf))
G = np.cos(np.radians(TH)) * np.sin(np.radians(TH)); G /= G.sum()


def one(a):
    name, d, h = a; no, ne = ETLS[name]
    S = cps2.Stack(550.0, (EMo, EMe), 25.0, 12.5, above=[(no, ne, d), (AG, AG, 100.0)],
                   below=[(TPo, TPe, 180.0), (ITO, ITO, 50.0)], n_sub=1.8, h=h)
    r = cps2.solve_pol(S, npts=8000); R = np.real(cps2.stack_reflectance(S, TH)); P = np.real(cps2.sub_angular(S, TH))
    es = float(np.real(r['air'] + r['sub'])); Pn = P / P.sum()
    out = []
    for k, (BT, BR) in EXT.items():
        e, t = series.eta_ext(BT, BR, R, P, 200)
        out.append([name, d, round(h, 3), k, es, float(Pn[TH > 60].sum()), float(1 - Pn @ R), float(BT @ Pn), float(BT @ G), e, es * e])
    return out


if __name__ == '__main__':
    jobs = [(n, float(d), h) for n in ('B3PyMPM', '3TPYMB', 'iso_no', 'iso_ne') for d in range(50, 425, 25) for h in (2 / 3, 0.8)]
    with Pool(4) as p: res = [x for o in p.map(one, jobs) for x in o]
    with open('etl_scatter.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['ETL', 'd_ETL', 'Theta', 'extractor', 'eta_sub', 'Psub_frac_gt60', 'Aprime_Psub', 'p1_device', 'p1_Lamb', 'eta_ext', 'EQE']); w.writerows(res)

"""ETL comparison (B3PyMPM vs 3TPYMB, as etl_sweep.py / etl_scatter.py) with a hexagonal cone array
(r 1.4 um, h 1.8 um, 0.2 um gap; sim/mla/cone_bsdf.npz) against the close-packed LightTools hemisphere MLA
and a hemisphere array with the same 0.2 um gap. -> etl_cone.csv"""
import numpy as np, scipy.io as sio, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla', '/home/user/OLED-/sim/etl_compare']
import cps2, series
from etl_sweep import ETLS, AG, ITO, EMo, EMe, TPo, TPe, TH

def split(B): return B[:90].sum(0), B[90:180][::-1]
C = np.load('/home/user/OLED-/sim/mla/cone_bsdf.npz')
EXT = {'hemisphere MLA (LT, close-packed)': split(sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat')['BSDF_MLA'][:, :, 10]),
       'hemisphere, 0.2 um gap': split(C['hemi_gap_n1.8']), 'cone r1.4 h1.8, 0.2 um gap': split(C['cone_n1.8'])}


def one(a):
    name, d, h = a; no, ne = ETLS[name]
    S = cps2.Stack(550.0, (EMo, EMe), 25.0, 12.5, above=[(no, ne, d), (AG, AG, 100.0)],
                   below=[(TPo, TPe, 180.0), (ITO, ITO, 50.0)], n_sub=1.8, h=h)
    r = cps2.solve_pol(S, npts=8000); R = np.real(cps2.stack_reflectance(S, TH)); P = np.real(cps2.sub_angular(S, TH))
    es = float(np.real(r['air'] + r['sub'])); Pn = P / P.sum()
    return [[name, d, round(h, 3), k, es, float(BT @ Pn), series.eta_ext(BT, BR, R, P, 200)[0], es * series.eta_ext(BT, BR, R, P, 200)[0]]
            for k, (BT, BR) in EXT.items()]


if __name__ == '__main__':
    jobs = [(n, float(d), h) for n in ('B3PyMPM', '3TPYMB') for d in range(50, 425, 25) for h in (2 / 3, 0.8)]
    with Pool(4) as p: res = [x for o in p.map(one, jobs) for x in o]
    with open('etl_cone.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['ETL', 'd_ETL', 'Theta', 'extractor', 'eta_sub', 'p1_device', 'eta_ext', 'EQE']); w.writerows(res)
    best = {}
    for r in res:
        k = (r[3], r[2], r[0])
        if k not in best or r[-1] > best[k][-1]: best[k] = r
    for th in (0.667, 0.8):
        for ex in EXT:
            a, b = best[(ex, th, '3TPYMB')], best[(ex, th, 'B3PyMPM')]
            print('Theta %.3f %-34s 3TPYMB %.3f (d %3.0f, p1 %.3f) | B3PyMPM %.3f (d %3.0f, p1 %.3f) | diff %+.1f %%p' % (
                th, ex, a[-1], a[1], a[5], b[-1], b[1], b[5], 100 * (b[-1] - a[-1])))
    bins = [(0, 30), (30, 45), (45, 60), (60, 70), (70, 90)]
    for ex, (BT, BR) in EXT.items():
        print('%-34s first-pass T per bin %s' % (ex, ' '.join('%d-%d: %.2f' % (lo, hi, BT[(TH > lo) & (TH < hi)].mean()) for lo, hi in bins)))

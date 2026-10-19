"""From the preprint Table 1 cell (B3PyMPM 200 / TAPC 180 / ITO 50 / n_sub 1.8 / LightTools hemisphere MLA, Theta 0.8,
Ir(ppy)2acac 400-700 nm) to the be_cone_opt.py device, one change at a time (2 nm spectral grid)."""
import numpy as np, scipy.io as sio, sys
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla']
import cps2, series
import be_cone_opt as B

LT = sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat')['BSDF_MLA']
ITO = B.lib('l_ITO'); C = np.load('/home/user/OLED-/sim/mla/cone_bsdf.npz')
def lt(n): k = int(round((n - 1.30) / 0.05)); b = LT[:, :, k]; return b[:90].sum(0), b[90:180][::-1]
EXT = {'LT1.8': lt(1.8), 'LT1.75': lt(1.75), 'hemi_cp1.8': B.split(C['hemi_cp_n1.8']) if 'hemi_cp_n1.8' in C else None,
       'hemi_cp1.77': B.split(C['hemi_cp_n1.77']), 'cone1.77': B.split(C['cone_n1.77'])}


def pt(a):
    tco, nsub, ext, dE, dH, j = a
    z = ITO[j] if tco == 'ITO' else B.IZOS['SJ_RTP'][j]
    S = cps2.Stack(float(B.LAM[j]), (B.EMo[j], B.EMe[j]), 25.0, 12.5, above=[(B.B3o[j], B.B3e[j], dE), (B.AG[j], B.AG[j], 100.0)],
                   below=[(B.TPo[j], B.TPe[j], dH), (z, z, 50.0)], n_sub=nsub, h=0.8)
    r = cps2.solve_pol(S, npts=3000); R = np.real(cps2.stack_reflectance(S, B.TH)); P = np.real(cps2.sub_angular(S, B.TH))
    es = float(np.real(r['air'] + r['sub'])); BT, BR = EXT[ext]
    return es, es * series.eta_ext(BT, BR, R, P, 150)[0]


STEPS = [('preprint Table 1 cell (ITO, n 1.8, LT hemisphere, 200/180)', 'ITO', 1.8, 'LT1.8', 200, 180),
         ('+ own-tracer close-packed hemisphere (n 1.8)', 'ITO', 1.8, 'hemi_cp1.8', 200, 180),
         ('+ substrate/lens n 1.77 (own hemisphere 1.77)', 'ITO', 1.77, 'hemi_cp1.77', 200, 180),
         ('+ IZO 50 (SJ_RTP) instead of ITO 50', 'IZO', 1.77, 'hemi_cp1.77', 200, 180),
         ('+ HTL 200 (= be_cone hemi_cp optimum under 200)', 'IZO', 1.77, 'hemi_cp1.77', 200, 200),
         ('+ cone array instead of hemisphere', 'IZO', 1.77, 'cone1.77', 200, 200)]
if __name__ == '__main__':
    if EXT['hemi_cp1.8'] is None:
        STEPS = [s for s in STEPS if s[3] != 'hemi_cp1.8']
    js = list(range(0, 301, 2)); w = B.WP[js] / B.WP[js].sum()
    with Pool(4) as p:
        for lab, tco, nsub, ext, dE, dH in STEPS:
            res = np.array(p.map(pt, [(tco, nsub, ext, float(dE), float(dH), j) for j in js]))
            print('%-62s eta_sub %.3f  EQE %.3f' % (lab, w @ res[:, 0], w @ res[:, 1]), flush=True)

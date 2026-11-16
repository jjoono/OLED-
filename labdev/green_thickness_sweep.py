"""Sensitivity of the green bare device to the emission-zone position z0 (from the HTL side of the EML) and the
horizontal dipole ratio; 20 nm spectral grid. Compares EQE, angular factor and I/I0 with the measured bare s2d2."""
import sys, json, numpy as np
sys.path.insert(0, '../selectivity'); import sel
from sel import cps2
import device as DV, simulate as SM, layout as LY
MEAS = np.array([1.000, 0.971, 0.924, 0.853, 0.771, 0.701, 0.575, 0.335, 0.114])
LAMS = np.arange(440, 661, 20).astype(float); J = (LAMS - 400).astype(int)
s = sel.D['s_Ir(ppy)2acac'][J] * LAMS; s = s / s.sum(); PL = SM.planar_kernel(); st0 = None
E = np.arange(91); dOm = 2 * np.pi * (np.cos(np.radians(E[:-1])) - np.cos(np.radians(E[1:])))
res = []
import itertools
for (etl, htl) in itertools.product((35.0, 42.0, 50.0, 60.0, 70.0), (90.0, 105.0, 120.0)):
    z0, th = 12.5, 0.77
    if True:
        A = np.zeros(90); w_sum = 0; Fs = []; fs = []; air_all = []
        for j, lam in zip(J, LAMS):
            org = lambda m, d: (sel.D['n_' + m + '_o'][j].real, sel.D['n_' + m + '_e'][j].real, d)
            met = lambda m, d: (sel.D['n_' + m][j], sel.D['n_' + m][j], d)
            eml = (sel.D['n_EML_o'][j].real, sel.D['n_EML_e'][j].real)
            below = [org('TCTA', 10.0), org('TAPC', htl), met('ITO', 150.0)]
            S = cps2.Stack(lam, eml, 25.0, z0, above=[org('B3', etl), met('LiF', 1.0), met('Al', 100.0)], below=below, n_sub=1.51, h=th)
            r = cps2.solve_pol(S, npts=3000); ref = cps2.solve_pol(cps2.Stack(lam, eml, 25.0, z0, above=[], below=[], n_sub=eml[0], n_top=eml[0], h=th), npts=3000)
            th_f = np.linspace(0, 90, 1801)[:-1] + 0.025; pf = np.clip(np.real(cps2.sub_angular(S, th_f)), 0, None).reshape(90, 20).sum(1); P = pf / pf.sum()
            R = np.clip(np.real(cps2.stack_reflectance(S, DV.TH)), 0, 1)
            Rd = {'pix': R, 'ITO_Al': R, 'noITO_Al': R, 'ITO_noAl': np.ones(90) * 0, 'noITO_noAl': np.ones(90) * 0}
            gst = np.load('layout_stacks_green.npz'); k = list(DV.LAMS).index(lam); Rd.update({n: gst[n][k] for n in gst.files})
            eS, eP, sd = LY.mc(PL, PL, Rd, P, 0.7, 0.0, True, N=60000, seed=j)
            Fs.append(r['P_tot'] / ref['P_tot']); fs.append(r['air'] + r['sub']); air_all.append(eP @ (PL['BT'] / np.maximum(PL['p'][:, None], 1e-15)))
        Fs = np.array(Fs); fs = np.array(fs); wq = 0.98 * s * Fs * fs / (0.98 * np.sum(s * Fs) + 0.02)
        Aa = (wq[:, None] * np.array(air_all)); prof = Aa.sum(0) / dOm; i0 = (prof[:5] * dOm[:5]).sum() / dOm[:5].sum()
        rel = np.array([1.0] + [float((prof[a - 3:a + 3] * dOm[a - 3:a + 3]).sum() / dOm[a - 3:a + 3].sum() / i0) for a in range(10, 90, 10)])
        eqe = Aa.sum(); pe = eqe / (np.pi * i0); rm = float(np.sqrt(np.mean((rel[1:] - MEAS[1:]) ** 2)))
        res.append(dict(etl=etl, htl=htl, z0=z0, theta_h=th, EQE=float(eqe), PE=float(pe), rmse=rm, rel=rel.tolist()))
        print(f'ETL {etl:4} HTL {htl:5}: EQE {eqe*100:.1f}%  PE {pe:.2f}  RMSE {rm:.3f}  I/I0 ' + ' '.join(f'{v:.2f}' for v in rel), flush=True)
json.dump(res, open('green_thickness_sweep.json', 'w'))

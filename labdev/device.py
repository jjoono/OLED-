"""Optics of the lab device ITO 150 / TAPC 50 / TCTA 10 / TCTA:B3PyMPM 1:1 + Ir(dmppy-ph)2tmd 4 wt% (25, dipole at
the centre) / B3PyMPM 45 / LiF 1 / Al 100 on glass (n = 1.51), horizontal dipole ratio 0.865 (CPS + TMM, group library).
Per wavelength (400-700 nm, 10 nm): fraction of dipole power entering the substrate f_sub, its angular distribution
P(theta_sub) in 1-degree bins, the Purcell factor F (total dissipated power relative to the dipole in a homogeneous
EML-like medium), and the back-face reflectance seen from the glass under the pixel (full stack) and outside it
(ITO / organics / air, no LiF/Al). Writes device.npz."""
import sys, numpy as np
sys.path.insert(0, '../selectivity'); import sel
from sel import cps2
NSUB = 1.51; LAMS = np.arange(400, 701, 10).astype(float); J = (LAMS - 400).astype(int)
TH = np.arange(90) + 0.5


def run(theta_h=0.865):
    out = {k: [] for k in ('f_sub', 'P', 'F', 'R_pix', 'R_noAl', 'wg', 'spp', 'abs')}
    for j, lam in zip(J, LAMS):
        org = lambda m, d: (sel.D['n_' + m + '_o'][j].real, sel.D['n_' + m + '_e'][j].real, d)
        met = lambda m, d: (sel.D['n_' + m][j], sel.D['n_' + m][j], d)
        eml = (sel.D['n_EML_o'][j].real, sel.D['n_EML_e'][j].real)
        below = [org('TCTA', 10.0), org('TAPC', 50.0), met('ITO', 150.0)]
        S = cps2.Stack(lam, eml, 25.0, 12.5, above=[org('B3', 45.0), met('LiF', 1.0), met('Al', 100.0)], below=below, n_sub=NSUB, h=theta_h)
        r = cps2.solve_pol(S, npts=4000)
        ref = cps2.solve_pol(cps2.Stack(lam, eml, 25.0, 12.5, above=[], below=[], n_sub=eml[0], n_top=eml[0], h=theta_h), npts=4000)
        th_f = np.linspace(0, 90, 1801)[:-1] + 0.025
        pf = np.clip(np.real(cps2.sub_angular(S, th_f)), 0, None).reshape(90, 20).sum(1)
        S0 = cps2.Stack(lam, eml, 25.0, 12.5, above=[org('B3', 45.0)], below=below, n_sub=NSUB, n_top=1.0, h=theta_h)
        out['f_sub'].append(r['air'] + r['sub']); out['P'].append(pf / pf.sum()); out['F'].append(r['P_tot'] / ref['P_tot'])
        out['R_pix'].append(np.clip(np.real(cps2.stack_reflectance(S, TH)), 0, 1))
        out['R_noAl'].append(np.clip(np.real(cps2.stack_reflectance(S0, TH)), 0, 1))
        for k in ('wg', 'spp', 'abs'): out[k].append(r[k])
        print(f"{lam:.0f} nm: f_sub {out['f_sub'][-1]:.3f}  F {out['F'][-1]:.3f}  R(10°) pix {out['R_pix'][-1][10]:.3f} noAl {out['R_noAl'][-1][10]:.3f}", flush=True)
    s = sel.D['s_Ir(dmppy-ph)2tmd'][J] * LAMS; s = s / s.sum()          # photon spectrum (intrinsic)
    np.savez('device.npz', lam=LAMS, s=s, **{k: np.array(v) for k, v in out.items()})


if __name__ == '__main__':
    run()

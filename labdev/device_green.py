"""Green lab device: ITO 150 / TAPC 105 / TCTA 10 / TCTA:B3PyMPM (1:1):Ir(ppy)2acac 4 wt% (25, dipole at the centre) /
B3PyMPM 42 / LiF 1 / Al 100 on glass (n = 1.51). Horizontal dipole ratio 0.77 (Ir(ppy)2acac in TCTA:B3PyMPM, literature
value), library Ir(ppy)2acac PL spectrum. Writes device_green.npz (same fields as device.npz) and layout_stacks_green.npz."""
import sys, numpy as np
sys.path.insert(0, '../selectivity'); import sel
from sel import cps2
import device as DV
TH_H = 0.77; HTL = 105.0; ETL = 42.0
out = {k: [] for k in ('f_sub', 'P', 'F', 'R_pix', 'R_noAl', 'wg', 'spp', 'abs')}; reg = {k: [] for k in ('ITO_Al', 'noITO_Al', 'ITO_noAl', 'noITO_noAl')}
for j, lam in zip(DV.J, DV.LAMS):
    org = lambda m, d: (sel.D['n_' + m + '_o'][j].real, sel.D['n_' + m + '_e'][j].real, d)
    met = lambda m, d: (sel.D['n_' + m][j], sel.D['n_' + m][j], d)
    eml = (sel.D['n_EML_o'][j].real, sel.D['n_EML_e'][j].real)
    below = [org('TCTA', 10.0), org('TAPC', HTL), met('ITO', 150.0)]
    cath = [org('B3', ETL), met('LiF', 1.0), met('Al', 100.0)]
    S = cps2.Stack(lam, eml, 25.0, 12.5, above=cath, below=below, n_sub=DV.NSUB, h=TH_H)
    r = cps2.solve_pol(S, npts=4000)
    ref = cps2.solve_pol(cps2.Stack(lam, eml, 25.0, 12.5, above=[], below=[], n_sub=eml[0], n_top=eml[0], h=TH_H), npts=4000)
    th_f = np.linspace(0, 90, 1801)[:-1] + 0.025
    pf = np.clip(np.real(cps2.sub_angular(S, th_f)), 0, None).reshape(90, 20).sum(1)
    out['f_sub'].append(r['air'] + r['sub']); out['P'].append(pf / pf.sum()); out['F'].append(r['P_tot'] / ref['P_tot'])
    out['R_pix'].append(np.clip(np.real(cps2.stack_reflectance(S, DV.TH)), 0, 1))
    for k in ('wg', 'spp', 'abs'): out[k].append(r[k])
    for ito in (True, False):
        bl = [org('TCTA', 10.0), org('TAPC', HTL)] + ([met('ITO', 150.0)] if ito else [])
        for al in (True, False):
            ab = cath if al else [org('B3', ETL)]
            S2 = cps2.Stack(lam, eml, 25.0, 12.5, above=ab, below=bl, n_sub=DV.NSUB, n_top=1.0, h=TH_H)
            reg[('ITO' if ito else 'noITO') + '_' + ('Al' if al else 'noAl')].append(np.clip(np.real(cps2.stack_reflectance(S2, DV.TH)), 0, 1))
    out['R_noAl'].append(reg['ITO_noAl'][-1])
s = sel.D['s_Ir(ppy)2acac'][DV.J] * DV.LAMS; s = s / s.sum()
np.savez('device_green.npz', lam=DV.LAMS, s=s, **{k: np.array(v) for k, v in out.items()})
np.savez('layout_stacks_green.npz', **{k: np.array(v) for k, v in reg.items()})
k = list(DV.LAMS).index(520.0)
print('520 nm: f_sub %.3f wg %.3f spp %.3f abs %.3f F %.3f | s-weighted f_sub %.3f F %.3f' % (out['f_sub'][k], out['wg'][k], out['spp'][k], out['abs'][k], out['F'][k], (s * np.array(out['f_sub'])).sum(), (s * np.array(out['F'])).sum()))

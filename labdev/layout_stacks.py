"""Back-face reflectance (seen from the glass) of the four regions of the lab substrate layout:
ITO/organics/LiF/Al, organics/LiF/Al (no ITO), ITO/organics/air, organics/air. PDL is neglected (thin insulator).
Writes layout_stacks.npz (R[lam, bin] for each region)."""
import sys, numpy as np
sys.path.insert(0, '../selectivity'); import sel
from sel import cps2
import device as DV
out = {k: [] for k in ('ITO_Al', 'noITO_Al', 'ITO_noAl', 'noITO_noAl')}
for j, lam in zip(DV.J, DV.LAMS):
    org = lambda m, d: (sel.D['n_' + m + '_o'][j].real, sel.D['n_' + m + '_e'][j].real, d)
    met = lambda m, d: (sel.D['n_' + m][j], sel.D['n_' + m][j], d)
    eml = (sel.D['n_EML_o'][j].real, sel.D['n_EML_e'][j].real)
    for ito in (True, False):
        below = [org('TCTA', 10.0), org('TAPC', 50.0)] + ([met('ITO', 150.0)] if ito else [])
        for al in (True, False):
            above = [org('B3', 45.0)] + ([met('LiF', 1.0), met('Al', 100.0)] if al else [])
            S = cps2.Stack(lam, eml, 25.0, 12.5, above=above, below=below, n_sub=DV.NSUB, n_top=1.0, h=0.865)
            out[('ITO' if ito else 'noITO') + '_' + ('Al' if al else 'noAl')].append(np.clip(np.real(cps2.stack_reflectance(S, DV.TH)), 0, 1))
np.savez('layout_stacks.npz', **{k: np.array(v) for k, v in out.items()})
d = np.load('layout_stacks.npz'); k = 15
print('550 nm, R at 10/50 deg:', {n: (round(float(d[n][k][10]), 3), round(float(d[n][k][50]), 3)) for n in d.files})

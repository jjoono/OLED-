"""Sensitivity of the reference (conventional) device to the anisotropy of the 60 nm B3PyMPM:Cs2CO3 layer next to the Ag
cathode: measured birefringent B3PyMPM (n_o, n_e) vs isotropic (n_e = n_o) vs n_o, n_e both shifted by -0.02 (dilution).
Planar air fraction (n 1.52), substrate delivery at n 1.77, SPP fraction; Ir(ppy)2acac photon weights, Theta 0.714."""
import sys, numpy as np
sys.path[:0] = ['/home/user/OLED-/sim/diag', '/home/user/OLED-/sim/fig3c']
import fit_ref_angles as F, cps2
D = F.D; f, g = 172.5 / 185, 160 / 210
def st(j, mode, ns):
    no, ne = D.B3o[j], D.B3e[j]
    if mode == 'iso': dn = (no, no)
    elif mode == 'shift': dn = (no - 0.02, ne - 0.02)
    elif mode == 'half': dn = (no, (no + ne) / 2)
    else: dn = (no, ne)
    return cps2.Stack(float(D.LAM[j]), (D.EMo[j], D.EMe[j]), 25.0, 1.0,
                      above=[(no, ne, 120 * f), dn + (60 * f,), (1.7, 1.7, 5.0), (D.AG[j], D.AG[j], 100.0)],
                      below=[(D.TCo[j], D.TCe[j], 10.0), (D.TPo[j], D.TPe[j], 70 * g), (D.HAT[j], D.HAT[j], 5.0), (D.TPo[j], D.TPe[j], 70 * g),
                             (D.HAT[j], D.HAT[j], 5.0), (D.TPo[j], D.TPe[j], 70 * g), (D.HAT[j], D.HAT[j], 5.0), (F.dev3.IZO[j], F.dev3.IZO[j], 12.7)],
                      n_sub=ns, h=0.714)
js = list(range(0, 301, 10)); w = (D.LAM * np.clip(np.asarray(D.spec.l_I_Irppy2acac, float)[:301], 0, None))[js]; w = w / w.sum()
print('B3PyMPM at 530 nm: n_o %.3f  n_e %.3f' % (D.B3o[130].real, D.B3e[130].real))
for mode in ('meas', 'half', 'iso', 'shift'):
    a = np.array([[np.real(cps2.solve_pol(st(j, mode, 1.52), npts=2000)['air']), *(np.real([cps2.solve_pol(st(j, mode, 1.77), npts=2000)[k] for k in ('air', 'sub', 'spp')]))] for j in js])
    print('%-6s planar air %.4f | n1.77: eta_sub %.4f  SPP %.4f' % (mode, w @ a[:, 0], w @ (a[:, 1] + a[:, 2]), w @ a[:, 3]))

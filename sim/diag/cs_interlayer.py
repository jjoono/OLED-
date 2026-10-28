"""Thin, nearly pure Cs2CO3 (or reacted Cs-Ag) interlayer between B3PyMPM:Cs2CO3 and the Ag cathode.
Stack as fit_ref_angles (ETL 172.5 nm, TAPC x0.762, IZO 12.7 nm, Theta 0.714, dipole 1 nm from the ETL side), the 5 nm
50 wt% layer replaced by d_cs nm of n 1.7 + i k.  Outputs at 400-700 nm (Ir(ppy)2acac weighted): eta_sub on n 1.77,
the substrate-side round-trip loss A' (cos-sin weighted over the substrate angles, s/p averaged) and R_s(30 deg in air, n 1.52)."""
import sys, numpy as np
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/diag', '/home/user/OLED-/sim/fig3c']
import fit_ref_angles as F, cps2
D = F.D; dev3 = F.dev3
f, g, z, di, th = 172.5 / 185, 160 / 210, 1.0, 12.7, 0.714
LJ = list(range(0, 301, 10)); TH = np.arange(90) + 0.5
I = np.clip(np.asarray(D.spec.l_I_Irppy2acac, float)[:301], 0, None); W = (D.LAM * I)[LJ]; W /= W.sum()
def stack(j, dcs, kcs, ns):
    B3o, B3e, EMo, EMe, TCo, TCe, TPo, TPe, HAT, AG = (D.B3o[j], D.B3e[j], D.EMo[j], D.EMe[j], D.TCo[j], D.TCe[j], D.TPo[j], D.TPe[j], D.HAT[j], D.AG[j])
    above = [(B3o, B3e, 120 * f), (B3o, B3e, 60 * f + 5 - dcs)] + ([(1.7 + 1j * kcs,) * 2 + (dcs,)] if dcs > 0 else []) + [(AG, AG, 100.0)]
    return cps2.Stack(float(D.LAM[j]), (EMo, EMe), 25.0, z, above=above,
                      below=[(TCo, TCe, 10.0), (TPo, TPe, 70 * g), (HAT, HAT, 5.0), (TPo, TPe, 70 * g), (HAT, HAT, 5.0), (TPo, TPe, 70 * g),
                             (HAT, HAT, 5.0), (dev3.IZO[j], dev3.IZO[j], di)], n_sub=ns, h=th)
def one(a):
    j, dcs, kcs = a
    S = stack(j, dcs, kcs, 1.77); r = cps2.solve_pol(S, npts=3000); es = float(np.real(r['air'] + r['sub']))
    R = np.real(cps2.stack_reflectance(S, TH)); w = np.cos(np.radians(TH)) * np.sin(np.radians(TH)); Ap = float(w @ (1 - R) / w.sum())
    S2 = stack(j, dcs, kcs, 1.52); k0 = 2 * np.pi / S2.lam; kx = np.sin(np.radians(30.0)) * k0
    rs, _ = cps2.stack_rt(S2.all_no, S2.all_ne, S2.all_d, np.array([kx]), k0, 's')
    return es, Ap, float(abs(rs[0]) ** 2), float(np.real(r['spp']))
if __name__ == '__main__':
    with Pool(4) as p:
        for dcs, kcs in ((0, 0), (1.5, 0.0), (1.5, 0.3), (1.5, 1.0), (1.5, 2.0), (3.0, 1.0), (5.0, 0.3)):
            o = np.array(p.map(one, [(j, dcs, kcs) for j in LJ]))
            print('d_cs %.1f nm  k %.1f : eta_sub(1.77) %.3f  A_round %.3f  R_s30 %.3f  SPP %.3f' % (dcs, kcs, *(W @ o)), flush=True)

"""EQE of the three devices (flat-glass reference, MLA without / with auxiliary back reflector) for the thickness sets
fitted to the reference angular EL (fit_ref_angles.py), compared with the nominal stack.  PLQY 1, no extra absorption."""
import sys, numpy as np
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/diag']
import dev3
IZO0 = dev3.IZO.copy()
# argv: "('name', f_ETL, Theta, (g_TAPC, dipole distance from ETL side [nm], IZO [nm], EML/TCTA/HATCN scale), [IZO n scale])"
if __name__ == '__main__':
    sets = [('nominal', 1.0, 0.76, (1.0, 12.5, 50.0, 1.0))] + [eval(a) for a in sys.argv[1:]]
    meas = np.array([0.16, 0.51, 0.59])
    for name, f, th, geo, *kz in sets:
        dev3.IZO = IZO0 * (kz[0] if kz else 1.0)        # IZO index scale (set before the pool forks)
        with Pool(4) as p:
            o = dev3.three(p, f=f, theta=th, geo=geo, N=150000, nrun=8)
        m = np.array([o['ref'], o['cathode'], o['full']])
        print('%-28s ref %.3f  MLA %.3f  MLA+refl %.3f | meas/model %s | MLA/ref %.2f refl/MLA %.2f (meas 3.19, 1.16)'
              % (name, *m, np.round(meas / m, 3), m[1] / m[0], m[2] / m[1]), flush=True)

"""Diagnostic: normal-direction EL spectrum of the flat-glass reference vs ETL tooling factor f (cavity fingerprint)."""
import numpy as np, sys
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/diag']
import dev3, cps2
D = dev3.D
JS = list(range(0, 301, 2))
def one(a):
    j, f = a
    B3o, B3e, EMo, EMe, TCo, TCe, TPo, TPe, HAT, AG = (D.B3o[j], D.B3e[j], D.EMo[j], D.EMe[j], D.TCo[j], D.TCe[j], D.TPo[j], D.TPe[j], D.HAT[j], D.AG[j])
    S = cps2.Stack(float(D.LAM[j]), (EMo, EMe), 25.0, 12.5,
                   above=[(B3o, B3e, 120.0 * f), (B3o, B3e, 60.0 * f), (1.7, 1.7, 5.0), (AG, AG, 100.0)],
                   below=[(TCo, TCe, 10.0), (TPo, TPe, 70.0), (HAT, HAT, 5.0), (TPo, TPe, 70.0), (HAT, HAT, 5.0), (TPo, TPe, 70.0),
                          (HAT, HAT, 5.0), (dev3.IZO[j], dev3.IZO[j], 50.0)], n_sub=1.52, h=0.76)
    P = np.real(cps2.sub_angular(S, np.array([0.5, 1.5, 2.5]))).mean()
    return P / np.real(cps2.solve(S, npts=3000)['P_tot'])     # photons per exciton into the substrate near normal
if __name__ == '__main__':
    I = np.clip(np.asarray(D.spec.l_I_Irppy2acac, float)[:301], 0, None)
    lam = D.LAM[JS]; pl = I[JS] / I[JS].max()
    print('PL (input) peak %d nm' % lam[pl.argmax()])
    with Pool(4) as p:
        for f in (0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.3, 1.35):
            cav = np.array(p.map(one, [(j, f) for j in JS]))
            el = pl * cav; el /= el.max(); half = lam[el >= 0.5]
            print('f %.2f (ETL %3.0f nm): normal EL peak %d nm, FWHM %d nm, I(560)/I(peak) %.2f, I(500)/I(peak) %.2f'
                  % (f, 185 * f, lam[el.argmax()], half.max() - half.min(), np.interp(560, lam, el), np.interp(500, lam, el)), flush=True)

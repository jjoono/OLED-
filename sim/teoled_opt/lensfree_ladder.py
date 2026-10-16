"""Ladder from the idealised lens-free limit to the real lens-free device (Theta 0.82, Ir(mphmq)2tmd spectrum, emission into air).
Each step re-optimises the free thicknesses (Nelder-Mead from several starts, 21-point spectrum; reported on 41 points)."""
import sys, numpy as np
from scipy.optimize import minimize
sys.path[:0] = ['/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/teoled_opt']
import cps2, tmodel as T, lensfree_limit as L
TH = L.TH; c0 = dict(theta=0.82, fwhm='Ir', lam0=620, metal='bulk', n_t=1.8, n_e=1.8, n_c=2.36, order='free')


def ideal(x, l, w, metal, dAg):
    return L.eqe([x[0], x[1], dAg if dAg else x[3], x[2]], dict(c0, metal=metal), l, w)


def real(x, l, w, agtop, al):
    dE, dH, z = x; s = 0
    for lam, wi in zip(l, w):
        above, below, i = T.stack(lam, dE, dH, 10., 7.87, agtop=agtop, al=al, tfe=False, zns=z)
        e = (np.real(T.BEB[i]),) * 2
        St = cps2.Stack(lam, e, 25.0, 12.5, above=above, below=below, n_sub=1.0, n_top=1.0, h=0.82)
        s += wi * np.sum(np.real(cps2.sub_angular(St, TH))) * 0.1 / np.real(cps2.solve(St, npts=3000)['P_tot'])
    return s


def opt(f, starts, bnds):
    best = (-1, None)
    for s in starts:
        r = minimize(lambda x: -f(x), s, method='Nelder-Mead', bounds=bnds, options=dict(maxiter=250, xatol=0.3, fatol=1e-5))
        if -r.fun > best[0]: best = (-r.fun, r.x)
    return best


def main():
    lc, wc = L.spec(c0, fine=False); lf, wf = L.spec(c0)
    steps = []
    # A: ideal, bulk Ag thickness free
    f = lambda x, l=lc, w=wc: ideal(x, l, w, 'bulk', None)
    e, x = opt(f, [[65, 70, 45, 28], [65, 70, 70, 20]], [(15, 400), (15, 400), (0, 250), (5, 60)])
    steps.append(('A ideal: bulk Ag, top Ag thickness free, n 1.8 organics, ZnS cap', ideal(x, lf, wf, 'bulk', None), x))
    # B: bulk Ag fixed 7.87 nm
    f = lambda x, l=lc, w=wc: ideal(x, l, w, 'bulk', 7.87)
    e, x = opt(f, [[65, 70, 45, 0], [65, 70, 225, 0], [65, 250, 45, 0]], [(15, 400), (15, 400), (0, 250), (0, 1)])
    steps.append(('B  + top Ag fixed at 7.87 nm (bulk nk)', ideal(x, lf, wf, 'bulk', 7.87), x))
    # C: measured joint nk 7.87 nm
    f = lambda x, l=lc, w=wc: ideal(x, l, w, 'joint8', 7.87)
    e, x = opt(f, [[65, 70, 45, 0], [65, 70, 225, 0], [65, 250, 45, 0]], [(15, 400), (15, 400), (0, 250), (0, 1)])
    steps.append(('C  + measured HATCN/Ag 8 nm nk (joint fit)', ideal(x, lf, wf, 'joint8', 7.87), x))
    # D: real materials (Ag/Al0.2/Liq/Bebq2/EML 25/TAPC o,e/HATCN 10/Ag 8 joint/ZnS), no TFE
    f = lambda x, l=lc, w=wc: real(x, l, w, L.AGJ, 0.2)
    e, x = opt(f, [[60, 60, 45], [80, 60, 225], [60, 250, 45], [100, 80, 100]], [(15, 400), (15, 400), (0, 250)])
    steps.append(('D  + real stack (Bebq2, EML 25, TAPC birefringent, HATCN 10, Al0.2/Liq)', real(x, lf, wf, L.AGJ, 0.2), x))
    for s, v, x in steps:
        print('%-80s %.3f  x=%s' % (s, v, np.round(x, 1)), flush=True)
    print('E  x eta_rad 0.96: %.3f' % (0.96 * steps[-1][1]))


if __name__ == '__main__':
    main()

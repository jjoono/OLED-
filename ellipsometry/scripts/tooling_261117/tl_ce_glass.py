"""Glass pieces with the CE-frame Gen-Osc dispersions: one glass, one backside
fraction, one thickness per piece (as tl_glass_global, normalised N,C,S), and the
glass Cauchy re-expressed as a Gen-Osc substrate (Einf + UV pole) so it can go
into a .mod."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_genosc as GO, tl_glass as GL, tl_glass_global as GG
import jnk_model as M
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
HC = 1239.841984
PIECES = {'1_31': 'NPB', '1_32': 'NPB', '1_33': 'NPB', '1_34': 'TAPC', '2_31': 'Bphen', '2_32': 'Bphen'}

def glass_genosc(A, B):
    """The glass Cauchy as a Gen-Osc substrate: Einf = 1 + one UV pole, i.e. a
    one-term Sellmeier.  Einf is held at 1 for the same reason as for the films -
    a free Einf drifted to 0.5, which no glass has."""
    wl = np.linspace(300, 1000, 200); E = HC / wl; n2 = (A + B / (wl / 1000)**2)**2
    q = least_squares(lambda q: 1.0 + q[0] / (q[1]**2 - E**2) - n2, [100.0, 10.0],
                      bounds=([0.0, 6.0], [2000.0, 30.0])).x
    nfit = np.sqrt(1.0 + q[0] / (q[1]**2 - E**2))
    return dict(einf=1.0, pole_amp=float(q[0]), pole_en=float(q[1])), float(np.max(np.abs(nfit - np.sqrt(n2))))


def main():
    CE = json.load(open(W + 'ce_genosc.json'))
    data = []
    for sh, mat in PIECES.items():
        wl, P, D, ang = L.se(sh, 300.0, 1000.0, 3)
        data.append((sh, wl, L.ncs(P, D), ang, GO.N_of(np.array(CE[mat]['pick']['p']), wl)))
    def piece(x, j):
        A, B, f = x[:3]; sh, wl, meas, ang, Nf = data[j]
        N, C_, S_ = GG.norm(*GL.ncs_glass(wl, Nf, x[3 + j], A + B / (wl / 1000)**2 + 0j, ang, f))
        return np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    res = lambda x: np.concatenate([piece(x, j) for j in range(len(data))])
    x0 = [1.51, 0.005, 0.5] + [{'1_31': 29, '1_32': 49, '1_33': 9.5, '1_34': 49, '2_31': 45.5, '2_32': 46.4}[s] for s in PIECES]
    lo = [1.40, -0.01, 0.0] + [0.5] * len(PIECES); hi = [1.70, 0.03, 1.0] + [150.0] * len(PIECES)
    r = least_squares(res, x0, bounds=(lo, hi), x_scale='jac', max_nfev=600)
    x = r.x; A, B, f = x[:3]
    per = {sh: float(M.mse_ncs(piece(x, j), 1 + 3 / len(PIECES))) for j, sh in enumerate(PIECES)}
    g, dn = glass_genosc(A, B)
    out = dict(cauchy=[float(A), float(B)], f=float(f), d={sh: float(x[3 + j]) for j, sh in enumerate(PIECES)},
               mse=per, mat=PIECES, glass_genosc=g, glass_genosc_max_dn=dn)
    json.dump(out, open(W + 'ce_glass.json', 'w'), indent=1)
    print('glass n(550) %.4f  backside fraction f = %.3f' % (A + B / 0.55**2, f))
    print('glass as Gen-Osc: Einf %.4f + UV pole %.2f eV^2 at %.2f eV  (max |dn| vs Cauchy %.1e)'
          % (g[0], g[1], g[2], out['glass_genosc_max_dn']))
    for sh in PIECES:
        print('  %s %-6s d = %6.2f nm  MSE %5.2f' % (sh, PIECES[sh], out['d'][sh], per[sh]))

if __name__ == '__main__':
    main()

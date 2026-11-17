"""Glass pieces, second pass: the collected back-surface fraction per piece.

The lab's tooling template (Bphen_S5) carries '% 1st Reflection' per sample in
its MultiSamp block, and freeing it per piece here drops 2_31 from MSE 8.2 to
3.5 and 2_34 from 7.1 to 4.7.  One glass (2-term Cauchy) is still shared by all
eight glass pieces.  The lab's no5 glass (n550 1.508 vs 1.521 here) fits the
Bphen glass pieces nearly as well (tl_s5.py: MSE 4.5 / 4.4 / 9.0 vs 3.8 / 4.5 /
8.8, thickness within 1 nm) but not 2_34 (8.3 vs 4.7), so this set's own glass
is kept."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_genosc as GO, tl_glass as GL, tl_glass_global as GG, tl_ce_glass as CG
import jnk_model as M
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
PIECES = {'1_31': 'NPB', '1_32': 'NPB', '1_33': 'NPB', '1_34': 'TAPC',
          '2_31': 'Bphen', '2_32': 'Bphen', '2_33': 'Bphen', '2_34': 'Bphen'}
D0 = {'1_31': 29, '1_32': 49, '1_33': 9, '1_34': 50, '2_31': 46, '2_32': 47, '2_33': 0.7, '2_34': 4.3}

def main():
    CE = json.load(open(W + 'ce_genosc.json'))
    data = []
    for sh, mat in PIECES.items():
        wl, P, D, ang = L.se(sh, 300.0, 1000.0, 3)
        data.append((wl, L.ncs(P, D), ang, GO.N_of(np.array(CE[mat]['pick']['p']), wl)))
    n = len(PIECES)
    def piece(x, j):
        A, B = x[:2]; wl, meas, ang, Nf = data[j]
        N, C_, S_ = GG.norm(*GL.ncs_glass(wl, Nf, x[2 + j], A + B / (wl / 1000)**2 + 0j, ang, x[2 + n + j]))
        return np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    res = lambda x: np.concatenate([piece(x, j) for j in range(n)])
    x0 = [1.51, 0.005] + [D0[s] for s in PIECES] + [0.7] * n
    lo = [1.40, -0.01] + [0.0] * n + [0.0] * n; hi = [1.70, 0.03] + [150.0] * n + [1.0] * n
    r = least_squares(res, x0, bounds=(lo, hi), x_scale='jac', max_nfev=800)
    x = r.x; A, B = x[:2]
    per = {sh: float(M.mse_ncs(piece(x, j), 2 + 2 / n)) for j, sh in enumerate(PIECES)}
    g, dn = CG.glass_genosc(A, B)
    out = dict(cauchy=[float(A), float(B)], d={sh: float(x[2 + j]) for j, sh in enumerate(PIECES)},
               f={sh: float(x[2 + n + j]) for j, sh in enumerate(PIECES)}, mse=per, mat=PIECES,
               glass_genosc=g, glass_genosc_max_dn=dn)
    json.dump(out, open(W + 'ce_glass2.json', 'w'), indent=1)
    print('shared glass n(550) %.4f;  as Gen-Osc: Einf 1 + UV pole %.2f eV^2 @ %.2f eV (max |dn| %.1e)'
          % (A + B / 0.55**2, g['pole_amp'], g['pole_en'], dn))
    for sh in PIECES:
        print('  %s %-6s d = %6.2f nm   backside %% 1st reflection = %5.1f   MSE %5.2f'
              % (sh, PIECES[sh], out['d'][sh], 100 * out['f'][sh], per[sh]))

if __name__ == '__main__':
    main()

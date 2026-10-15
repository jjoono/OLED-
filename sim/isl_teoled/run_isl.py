"""Sweep: planar TEOLED vs ISL TEOLED, top = (TFE + lens 1.7) / lens 1.7 / lens 1.5; ISL matrix 1.85 or 1.6;
S = 0..30, g 0.6, albedo 1 or 0.995.  Emitter-to-electrode distances optimised per case (PLQY = 1, Theta 0.76)."""
import sys, csv, itertools, numpy as np
sys.path[:0] = ['/home/user/OLED-/sim/scatter']
import isl_model as I
from mla_plus_scatter import slab

bs = np.load('mla_hemi_bsdf.npz')
MLA = {n: (bs[f'BT_{n}'], bs[f'BR_{n}']) for n in (1.5, 1.6, 1.7)}
TOPS = {'TFE+lens1.7': (True, 1.7), 'lens1.7': (False, 1.7), 'lens1.5': (False, 1.5)}
DE = (40., 60., 80., 120., 160., 200.); DH = (30., 50., 80., 120.)
NM = (1.85, 1.6); SS = (0, 1, 3, 10, 30); G = 0.6; ALB = (1.0, 0.995)


def top_layers(tfe, dH):
    L = [(I.TAPC[0], I.TAPC[1], dH), I.iso(I.HAT, 10.), I.iso(I.AG15, 15.)]
    return L + ([I.iso(n, d) for n, d in I.TFE] if tfe else [])


def coh(L3):            # (n_o,n_e,d) -> (n,d) isotropic n_o for the incoherent-light transfer
    return [(l[0], l[2]) for l in L3]


def main():
    slabs = {}
    for S, alb in itertools.product(SS, ALB):
        if S == 0:
            slabs[(S, alb)] = (np.eye(90), np.zeros((90, 90)))
        else:
            slabs[(S, alb)] = slab(S, G, N=6000, seed=S * 7 + int(alb * 1000), albedo=alb)

    rows = []
    for top, (tfe, nL) in TOPS.items():
        BT, BR = MLA[nL]
        best_p = None
        for dE, dH in itertools.product(DE, DH):
            tl = top_layers(tfe, dH)
            # planar: Ag 100 behind the ETL
            bel = [I.iso(I.BEB, dE), I.iso(I.AG, 100.)]
            _, sU, _, _ = I.emission(bel, tl, 1.0, nL)
            Ru, _ = I.rt_unpol(coh(tl[::-1]) + [(I.EML[0], 25.)] + coh(bel), nL, 1.0, I.TH)
            e = I.solve_planar(sU, Ru, BT, BR)
            rows.append([top, 'planar', '', '', '', dE, dH, sU.sum(), 0, e])
            if best_p is None or e > best_p[-1]: best_p = (dE, dH, sU.sum(), e)
            for nm in NM:
                bel = [I.iso(I.BEB, dE), I.iso(I.IZO, 50.)]
                sD, sU, _, _ = I.emission(bel, tl, nm, nL)
                stackD = coh(bel[::-1]) + [(I.EML[0], 25.)] + coh(tl)
                Rd, T_DU = I.transfer(stackD, nm, nL)
                Ru, T_UD = I.transfer(stackD[::-1], nL, nm)
                Mg = I.mirror_R(nm)
                for (S, alb), (Ts, Rs) in slabs.items():
                    Reff = Rs + Ts @ (Mg[:, None] * np.linalg.solve(np.eye(90) - Rs * Mg[None, :], Ts))
                    e = I.solve_isl(sD, sU, Reff, Rd, T_DU, Ru, T_UD, BT, BR)
                    rows.append([top, 'ISL', nm, S, alb, dE, dH, sU.sum(), sD.sum(), e])
        print(top, 'planar best', best_p, flush=True)

    with open('isl_sweep.csv', 'w', newline='') as f:
        csv.writer(f).writerows([['top', 'device', 'n_m', 'S', 'albedo', 'd_ETL', 'd_HTL', 'emit_up', 'emit_down', 'EQE']] + rows)

    # best per case
    print('\n%-12s %-8s %5s %3s %6s | %6s %5s %5s | %5s %5s | %s' % ('top', 'device', 'n_m', 'S', 'alb', 'EQE', 'dE', 'dH', 'up', 'down', 'lost at birth'))
    key = lambda r: (r[0], r[1], r[2], r[3], r[4])
    best = {}
    for r in rows:
        k = key(r)
        if k not in best or r[-1] > best[k][-1]: best[k] = r
    for k, r in best.items():
        print('%-12s %-8s %5s %3s %6s | %6.3f %5.0f %5.0f | %5.3f %5.3f | %5.3f' % (r[0], r[1], r[2], r[3], r[4], r[9], r[5], r[6], r[7], r[8], 1 - r[7] - r[8]))


if __name__ == "__main__":
    main()

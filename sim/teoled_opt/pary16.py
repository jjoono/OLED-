"""Practical TEOLED: ... HATCN 10 / Ag 8 (joint nk, 7.87 nm) / ZnS z / Al2O3 30 / ParyC 500 / Al2O3 30 / ParyC 500 / Al2O3 30 /
n=1.6 hemisphere MLA. Ag/Al0.2/Liq bottom, Bebq2 ETL 60-140, TAPC 40-160, ZnS 100-250. 27 wavelengths 560-716 nm.
Also: (a) before TFE (Ag/ZnS/air, planar), (b) after TFE without lens (planar, into air)."""
import csv, itertools, numpy as np, tmodel as T, cps2
from multiprocessing import Pool
T.NL = 1.6
bs = np.load('/home/user/OLED-/sim/isl_teoled/mla_hemi_bsdf.npz'); BT, BR = bs['BT_1.6'], bs['BR_1.6']
a = np.genfromtxt('nk_joint/Ag_HATCN_2-6_8nm_nk.csv', delimiter=',', names=True)
AG8 = np.interp(T.LAMS, a['wavelength_nm'], a['n_joint']) + 1j * np.interp(T.LAMS, a['wavelength_nm'], a['k_joint'])
LF = tuple(np.arange(560, 717, 6)); W = np.array([T.S[T.idx(l)] for l in LF]); W /= W.sum()
KW = dict(agtop=AG8, al=0.2)
def eqe_mla(dE, dH, z):
    eo = [T.eta_out(l, dE, dH, 10, 7.87, BT, BR, tfe=True, zns=z, **KW)[0] for l in LF]
    return T.ETA_RAD * float(W @ np.array(eo))
def eqe_air(dE, dH, z, tfe):
    out = []
    for l in LF:
        above, below, i = T.stack(l, dE, dH, 10, 7.87, tfe=tfe, zns=z, **KW)
        e = (np.real(T.BEB[i]),) * 2
        St = cps2.Stack(l, e, 25.0, 12.5, above=above, below=below, n_sub=1.0, n_top=1.0, h=T.THETA)
        P = np.real(cps2.solve(St)['P_tot'])
        out.append(np.nan_to_num(np.real(cps2.sub_angular(St, T.TH))).sum() / P)
    return T.ETA_RAD * float(W @ np.array(out))
DE = range(60, 141, 10); DH = range(40, 161, 10); ZS = (100, 125, 150, 175, 200, 225, 250)
def g(a): z, dE, dH = a; return (z, dE, dH, eqe_mla(dE, dH, z))
def h(a): z, dE, dH, e = a; return (z, dE, dH, e, eqe_air(dE, dH, z, False), eqe_air(dE, dH, z, True))
if __name__ == '__main__':
    with Pool(4) as p:
        G = p.map(g, list(itertools.product(ZS, DE, DH)), chunksize=6)
        with open('pary16_grid.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['ZnS', 'd_ETL', 'd_HTL', 'EQE_MLA']] + G)
        best = {}
        for r in G:
            if r[0] not in best or r[3] > best[r[0]][3]: best[r[0]] = r
        top = sorted(G, key=lambda r: -r[3])[:3]
        sel = list(best.values()) + [r for r in top if r not in best.values()]
        H = p.map(h, sel)
    with open('pary16_best.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['ZnS', 'd_ETL', 'd_HTL', 'EQE_MLA', 'EQE_air_beforeTFE', 'EQE_air_afterTFE_noLens']] + H)
    for r in H: print('ZnS %3d ETL %3d HTL %3d | MLA %.3f | before TFE %.3f | TFE no lens %.3f' % r)

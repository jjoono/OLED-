"""Top electrode IZO instead of HATCN/Ag/ZnS: ... TAPC / HATCN 10 / IZO d / Al2O3 30 / ParyC 500 / Al2O3 30 / ParyC 500 / Al2O3 30 /
n1.6 MLA; Ag/Al0.2/Liq bottom; l_IZO (2.044+0.0053i @550). Grid ETL 60-140, HTL 40-160, IZO 30-130; 27 wavelengths."""
import csv, itertools, numpy as np, tmodel as T, pary16 as P
from multiprocessing import Pool
LF, W = P.LF, P.W
def eqe_mla(dE, dH, dI, nk=None):
    eo = [T.eta_out(l, dE, dH, 10, dI, P.BT, P.BR, tfe=True, topm='IZO', izo_nk=nk, al=0.2)[0] for l in LF]
    return T.ETA_RAD * float(W @ np.array(eo))
def eqe_air(dE, dH, dI, tfe, nk=None):
    import cps2
    out = []
    for l in LF:
        above, below, i = T.stack(l, dE, dH, 10, dI, tfe=tfe, topm='IZO', izo_nk=nk, al=0.2)
        e = (np.real(T.BEB[i]),) * 2
        St = cps2.Stack(l, e, 25.0, 12.5, above=above, below=below, n_sub=1.0, n_top=1.0, h=T.THETA)
        out.append(np.nan_to_num(np.real(cps2.sub_angular(St, T.TH))).sum() / np.real(cps2.solve(St)['P_tot']))
    return T.ETA_RAD * float(W @ np.array(out))
DE = range(60, 141, 10); DH = range(40, 161, 10); DI = (30, 50, 70, 100, 130)
def g(a): dI, dE, dH = a; return (dI, dE, dH, eqe_mla(dE, dH, dI))
def h(a):
    dI, dE, dH, e = a
    lowk = np.asarray(T.nk('l_IZO_SJ_RTP'))
    return (dI, dE, dH, e, eqe_air(dE, dH, dI, False), eqe_air(dE, dH, dI, True), eqe_mla(dE, dH, dI, lowk))
if __name__ == '__main__':
    with Pool(4) as p:
        G = p.map(g, list(itertools.product(DI, DE, DH)), chunksize=6)
        with open('izo_top_grid.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['IZO', 'd_ETL', 'd_HTL', 'EQE_MLA']] + G)
        best = {}
        for r in G:
            if r[0] not in best or r[3] > best[r[0]][3]: best[r[0]] = r
        H = p.map(h, list(best.values()))
    with open('izo_top_best.csv', 'w', newline='') as fh: csv.writer(fh).writerows([['IZO', 'd_ETL', 'd_HTL', 'EQE_MLA', 'EQE_air_beforeTFE', 'EQE_air_afterTFE_noLens', 'EQE_MLA_IZO_SJ_RTP']] + H)
    for r in H: print('IZO %3d ETL %3d HTL %3d | MLA %.3f | before TFE %.3f | TFE no lens %.3f | MLA low-k IZO %.3f' % r)

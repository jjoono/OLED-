"""Bottom-emitting OLED on an n 1.77 substrate with an index-matched hexagonal cone array
(r 1.4 um, h 1.8 um, 0.2 um gap; sim/mla/cone_bsdf.npz 'cone_n1.77'):
  Ag 100 (McPeak) | B3PyMPM d_E (o/e) | TCTA:B3PyMPM:Ir(ppy)2acac 25 (dipole at centre, Theta 0.8) | TAPC d_H (o/e)
  | IZO 50 | glass n 1.77 (incoherent) | cone array | air
EQE = sum_lambda w(lambda) eta_sub(lambda) eta_ext(lambda), PLQY 1, w = photon-number weighted Ir(ppy)2acac spectrum
(400-700 nm); eta_ext from the full recycling series with the ray-traced BSDF.
Grid d_E, d_H = 20..500 nm step 10 on a 10 nm spectral grid; the best point under each bound (d_E, d_H <= 100/200/300/500)
is then re-evaluated on a 2 nm spectral grid (with a +-10 nm local search on that grid).
IZO: l_IZO_SJ_RTP (main) and l_IZO.  Hemisphere arrays (close-packed and 0.2 um gap, same tracer) and a flat
substrate are evaluated on the same grid for reference.  -> be_cone_grid.csv, be_cone_best.csv
ETL from the environment: BE_ETL=3TPYMB uses the isotropic library 3TPYMB (TPYMB_3) -> be_cone_grid_3TPYMB.csv, ..."""
import numpy as np, scipy.io as sio, sys, csv, itertools, os
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla']
import cps2, series

M = sio.loadmat('/home/user/OLED-/sim/table1_author/nk_JH_total.mat', squeeze_me=True, struct_as_record=False)
mat, spec = M['material'], M['spectrum']
LAM = np.arange(400, 701)
def lib(f): return np.asarray(getattr(mat, f))[:301]
ETL = os.environ.get('BE_ETL', 'B3PyMPM'); SUF = '' if ETL == 'B3PyMPM' else '_' + ETL
AG = lib('l_Ag_McPeak'); B3o, B3e = (lib('l_B3_o_JO'), lib('l_B3_e_JO')) if ETL == 'B3PyMPM' else (lib('TPYMB_3'),) * 2; EMo, EMe = lib('l_TCTA_B3_o_JO'), lib('l_TCTA_B3_e_JO')
TPo, TPe = lib('l_TAPC_o_JO'), lib('l_TAPC_e_JO'); IZOS = {'SJ_RTP': lib('l_IZO_SJ_RTP'), 'l_IZO': lib('l_IZO')}
I = np.clip(np.asarray(spec.I_Irppy2acac, float)[:301], 0, None); WP = LAM * I
NSUB, THETA = 1.77, 0.8
TH = np.arange(90) + 0.5
def split(B): return B[:90].sum(0), B[90:180][::-1]
C = np.load('/home/user/OLED-/sim/mla/cone_bsdf.npz')
c = np.cos(np.radians(TH)); s2 = (NSUB * np.sin(np.radians(TH))) ** 2; ok = s2 < 1; ct = np.sqrt(np.clip(1 - s2, 0, 1))
Rf = np.where(ok, 0.5 * (((NSUB * c - ct) / (NSUB * c + ct)) ** 2 + ((NSUB * ct - c) / (NSUB * ct + c)) ** 2), 1.0)
EXT = {'cone': split(C['cone_n1.77']), 'hemi_cp': split(C['hemi_cp_n1.77']), 'hemi_gap': split(C['hemi_gap_n1.77']),
       'flat': (1 - Rf, np.diag(Rf))}


def point(a):
    dE, dH, izo, j = a
    z = IZOS[izo][j]; lam = float(LAM[j])
    S = cps2.Stack(lam, (EMo[j], EMe[j]), 25.0, 12.5, above=[(B3o[j], B3e[j], dE), (AG[j], AG[j], 100.0)],
                   below=[(TPo[j], TPe[j], dH), (z, z, 50.0)], n_sub=NSUB, h=THETA)
    r = cps2.solve_pol(S, npts=3000); R = np.real(cps2.stack_reflectance(S, TH)); P = np.real(cps2.sub_angular(S, TH))
    es = float(np.real(r['air'] + r['sub']))
    return a, es, [series.eta_ext(BT, BR, R, P, 150)[0] for BT, BR in EXT.values()], float(np.real(r['spp'])), float(np.real(r['abs']))


def evaluate(pts, step, pool):
    js = list(range(0, 301, step)); w = WP[js] / WP[js].sum()
    res = pool.map(point, [(dE, dH, izo, j) for (dE, dH, izo) in pts for j in js], chunksize=16)
    out = {}
    for k, (dE, dH, izo) in enumerate(pts):
        blk = res[k * len(js):(k + 1) * len(js)]
        es = np.array([b[1] for b in blk]); ee = np.array([b[2] for b in blk])
        out[(dE, dH, izo)] = dict(eta_sub=float(w @ es), **{f'EQE_{x}': float(w @ (es * ee[:, i])) for i, x in enumerate(EXT)},
                                  spp=float(w @ np.array([b[3] for b in blk])))
    return out


if __name__ == '__main__':
    D = list(range(20, 501, 10))
    with Pool(4) as pool:
        grid = evaluate([(dE, dH, izo) for izo in IZOS for dE in D for dH in D], 10, pool)
        with open('be_cone_grid%s.csv' % SUF, 'w', newline='') as f:
            wr = csv.writer(f); keys = list(next(iter(grid.values())).keys()); wr.writerow(['IZO', 'd_ETL', 'd_HTL'] + keys)
            for (dE, dH, izo), v in grid.items(): wr.writerow([izo, dE, dH] + [v[k] for k in keys])
        rows = []
        for izo, ext, bound in itertools.product(IZOS, ('cone', 'hemi_cp', 'hemi_gap', 'flat'), (100, 200, 300, 500)):
            cand = {k: v for k, v in grid.items() if k[2] == izo and k[0] <= bound and k[1] <= bound}
            (dE, dH, _), _ = max(cand.items(), key=lambda kv: kv[1][f'EQE_{ext}'])
            loc = [(e, h, izo) for e in range(max(20, dE - 10), min(bound, dE + 10) + 1, 5) for h in range(max(20, dH - 10), min(bound, dH + 10) + 1, 5)]
            fine = evaluate(loc, 2, pool)
            (bE, bH, _), v = max(fine.items(), key=lambda kv: kv[1][f'EQE_{ext}'])
            row = dict(IZO=izo, extractor=ext, bound=bound, d_ETL=bE, d_HTL=bH, EQE=v[f'EQE_{ext}'], eta_sub=v['eta_sub'],
                       eta_ext=v[f'EQE_{ext}'] / v['eta_sub'], spp=v['spp'], EQE_grid10=grid[(dE, dH, izo)][f'EQE_{ext}'])
            rows.append(row); print({k: (round(x, 4) if isinstance(x, float) else x) for k, x in row.items()}, flush=True)
    with open('be_cone_best%s.csv' % SUF, 'w', newline='') as f:
        wr = csv.DictWriter(f, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)

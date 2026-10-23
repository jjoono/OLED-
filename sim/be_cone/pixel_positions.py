"""EQE of the 8 pixels of the mask set on the 24 x 24 x 0.5 mm, n 1.77 cone-array substrate, back fully covered by a
reflector, side walls absorbing.  Pixel centres (21 x 23 pattern centred on the substrate): x = +-2.48 mm,
y = +-1.5 mm ("centre" pixels) and +-4.5 mm ("edge" pixels); pixel 2.0 x 1.96 mm (1.47 x 1.96 mm at x +-2.21 if the strip
was trimmed).  Back reflector outside the pixels: (B1) same reflectance as the OLED stack, (B2) Ag directly on the glass.
By symmetry one pixel of each kind is traced."""
import numpy as np, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/be_cone']
import finite_substrate as F

AG = F.D.AG; NS = F.NS; TH = F.TH
def r_ag(j):                                     # glass (n 1.77) / Ag, unpolarised, per 1-deg bin
    n2 = AG[j]; ci = np.cos(np.radians(TH)); si = NS * np.sin(np.radians(TH)); ct = np.sqrt(1 - (si / n2) ** 2 + 0j)
    rs = (NS * ci - n2 * ct) / (NS * ci + n2 * ct); rp = (n2 * ci - NS * ct) / (n2 * ci + NS * ct)
    return 0.5 * (abs(rs) ** 2 + abs(rp) ** 2)


def run(args):
    seed, N, cx, cy, hx, hy, mode = args
    rng = np.random.default_rng(seed)
    j = rng.choice(len(WL), N, p=WL); Pc = np.cumsum(PS / PS.sum(1, keepdims=True), 1)
    ib = (Pc[j] < rng.random(N)[:, None]).sum(1).clip(0, 89)
    x = cx + (rng.random(N) - 0.5) * 2 * hx; y = cy + (rng.random(N) - 0.5) * 2 * hy
    ph = rng.random(N) * 2 * np.pi; w = np.ones(N); top = edge = 0.0; alive = np.ones(N, bool)
    for _ in range(600):
        a = np.flatnonzero(alive)
        if a.size == 0: break
        for half in (0, 1):
            th = np.radians(TH[ib[a]] + rng.random(a.size) - 0.5); dd = F.T * np.tan(th)
            x[a] += dd * np.cos(ph[a]); y[a] += dd * np.sin(ph[a])
            out = (np.abs(x[a]) > F.HALF) | (np.abs(y[a]) > F.HALF); edge += w[a[out]].sum(); w[a[out]] = 0
            if half == 0:
                bt = F.BT[ib[a]]; top += (w[a] * bt).sum(); w[a] *= 1 - bt; u = rng.random(a.size)
                ib[a] = (F.BRc[:, ib[a]].T < u[:, None]).sum(1).clip(0, 89); ph[a] = rng.random(a.size) * 2 * np.pi
        Rd = RS[j[a], ib[a]]
        if mode == 'B2':                         # OLED stack on all 8 pixel areas, Ag on glass elsewhere
            onpix = np.zeros(a.size, bool)
            for px in PIXELS:
                onpix |= (np.abs(x[a] - px[0]) <= px[2]) & (np.abs(y[a] - px[1]) <= px[3])
            Rd = np.where(onpix, Rd, RAG[j[a], ib[a]])
        w[a] *= Rd; alive &= w > 1e-5
    return top / N, edge / N


if __name__ == '__main__':
    with Pool(4) as p: res = p.map(F.stack, F.LAMJ)
    ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
    RAG = np.array([r_ag(j) for j in F.LAMJ])
    I = np.clip(np.asarray(F.D.spec.l_I_Irppy2acac, float)[:301], 0, None); wp = (F.D.LAM * I)[F.LAMJ]; wp /= wp.sum()
    es = float(wp @ ES); WL = wp * ES / (wp @ ES)
    print('glass/Ag reflectance at 530 nm, 0/40/70 deg: %s ; OLED stack: %s'
          % (np.round(RAG[13][[0, 40, 70]], 3), np.round(RS[13][[0, 40, 70]], 3)))
    rows = []
    for geo, hx, xc in (('as designed (2.0 x 1.96 mm)', 1.0, 2.48), ('strip trimmed (1.47 x 1.96 mm)', 0.735, 2.21)):
        PIXELS = [(sx * xc, sy * yc, hx, 0.98) for sx in (-1, 1) for sy in (-1, 1) for yc in (1.5, 4.5)]
        for mode in ('B1', 'B2'):
            for lab, yc in (('centre pixel (y +-1.5)', 1.5), ('edge pixel (y +-4.5)', 4.5)):
                with Pool(4) as p:
                    out = np.mean(p.map(run, [(100 * k + 5, 200000, xc, yc, hx, 0.98, mode) for k in range(8)]), 0)
                rows.append(dict(geometry=geo, back=mode, pixel=lab, eta_ext=out[0], EQE=es * out[0], edge_loss=out[1]))
                print('%-32s %s %-24s eta_ext %.4f  EQE %.4f  edge loss %.4f' % (geo, mode, lab, out[0], es * out[0], out[1]), flush=True)
    with open('pixel_positions.csv', 'w', newline='') as fh:
        wr = csv.DictWriter(fh, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)

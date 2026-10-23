"""No auxiliary reflector: only the cathode metal of the mask set reflects at the back.  Back-face map (0.05 mm raster):
inside the cathode pattern (fingers, C-leads, middle leads, anode-contact pads; geometry of figures/mask_check/align_etch.py,
21 x 23 pattern centred on the 24 x 24 substrate) -> OLED-stack reflectance; elsewhere bare glass/air (Fresnel/TIR, the
transmitted part leaves through the back).  Compared with: pixel-only Ag (A) and a full back reflector (B1)."""
import numpy as np, sys, csv
from multiprocessing import Pool
from matplotlib.path import Path
sys.path[:0] = ['/home/user/OLED-/sim/be_cone']
import finite_substrate as F

CXP, CYP = 10.5, 11.5                                     # pattern centre -> substrate origin
R = lambda x0, y0, x1, y1: [(x0, y0), (x1, y0), (x1, y1), (x0, y1)]
mir = lambda p: [(2 * 10.45 - x, y) for x, y in p]
TL = [(3.96, 0), (6.91, 0), (6.91, 2.5), (5.97, 2.5), (5.97, 6.0), (9.98, 6.0), (9.98, 8.0), (3.96, 8.0)]
CATP = [TL, R(1.2, 8.99, 9.98, 10.95), R(1.2, 12.05, 9.98, 14.01), [(x, 23 - y) for x, y in TL]]
CATP = CATP + [mir(p) for p in CATP] + [R(8.0, 0, 12.9, 2.5), R(8.0, 20.5, 12.9, 23)]
DX = 0.05; GX = np.arange(-12, 12, DX) + DX / 2
XX, YY = np.meshgrid(GX, GX, indexing='ij'); PTS = np.c_[XX.ravel(), YY.ravel()]
MET = np.zeros(len(PTS), bool)
for p in CATP: MET |= Path([(x - CXP, y - CYP) for x, y in p]).contains_points(PTS)
MET = MET.reshape(XX.shape)


def run(args):
    seed, N, cx, cy, mode = args
    rng = np.random.default_rng(seed)
    j = rng.choice(len(WL), N, p=WL); Pc = np.cumsum(PS / PS.sum(1, keepdims=True), 1)
    ib = (Pc[j] < rng.random(N)[:, None]).sum(1).clip(0, 89)
    x = cx + (rng.random(N) - 0.5) * 2.0; y = cy + (rng.random(N) - 0.5) * 1.96
    ph = rng.random(N) * 2 * np.pi; w = np.ones(N); top = back = edge = 0.0; alive = np.ones(N, bool)
    for _ in range(600):
        a = np.flatnonzero(alive)
        if a.size == 0: break
        for half in (0, 1):
            th = np.radians(F.TH[ib[a]] + rng.random(a.size) - 0.5); dd = F.T * np.tan(th)
            x[a] += dd * np.cos(ph[a]); y[a] += dd * np.sin(ph[a])
            out = (np.abs(x[a]) > F.HALF) | (np.abs(y[a]) > F.HALF); edge += w[a[out]].sum(); w[a[out]] = 0
            if half == 0:
                bt = F.BT[ib[a]]; top += (w[a] * bt).sum(); w[a] *= 1 - bt; u = rng.random(a.size)
                ib[a] = (F.BRc[:, ib[a]].T < u[:, None]).sum(1).clip(0, 89); ph[a] = rng.random(a.size) * 2 * np.pi
        Rd = RS[j[a], ib[a]]
        if mode == 'B1':
            Rb = Rd
        else:
            if mode == 'cathode':
                ix = np.clip(((x[a] + 12) / DX).astype(int), 0, len(GX) - 1); iy = np.clip(((y[a] + 12) / DX).astype(int), 0, len(GX) - 1)
                metal = MET[ix, iy]
            else:                                                  # 'A': metal only on the lit pixel
                metal = (np.abs(x[a] - cx) <= 1.0) & (np.abs(y[a] - cy) <= 0.98)
            Rb = np.where(metal, Rd, F.RF[ib[a]]); back += (w[a] * np.where(metal, 0, 1 - F.RF[ib[a]])).sum()
        w[a] *= Rb; alive &= w > 1e-5
    return top / N, back / N, edge / N


if __name__ == '__main__':
    with Pool(4) as p: res = p.map(F.stack, F.LAMJ)
    ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
    I = np.clip(np.asarray(F.D.spec.l_I_Irppy2acac, float)[:301], 0, None); wp = (F.D.LAM * I)[F.LAMJ]; wp /= wp.sum()
    es = float(wp @ ES); WL = wp * ES / (wp @ ES)
    print('back covered by cathode metal: %.1f %% of 24 x 24 mm (%.0f mm^2)' % (100 * MET.mean(), MET.sum() * DX * DX))
    for rr in (3, 5, 8):
        sel = (np.hypot(XX - (-2.48), YY - 1.5) < rr)
        print('   metal fraction within %d mm of the centre pixel: %.0f %%' % (rr, 100 * MET[sel].mean()))
    rows = []
    for mode, lab in (('A', 'Ag only on the lit pixel'), ('cathode', 'cathode pattern only (no auxiliary reflector)'), ('B1', 'full back reflector')):
        for plab, cy in (('centre pixel', 1.5), ('edge pixel', 4.5)):
            with Pool(4) as p: out = np.mean(p.map(run, [(100 * k + 3, 200000, -2.48, cy, mode) for k in range(8)]), 0)
            rows.append(dict(back=lab, pixel=plab, eta_ext=out[0], EQE=es * out[0], back_loss=out[1], edge_loss=out[2]))
            print('%-48s %-13s eta_ext %.4f  EQE %.4f | back loss %.4f  edge loss %.4f' % (lab, plab, out[0], es * out[0], out[1], out[2]), flush=True)
    with open('cathode_as_reflector.csv', 'w', newline='') as fh:
        wr = csv.DictWriter(fh, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)

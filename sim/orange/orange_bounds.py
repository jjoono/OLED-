"""Orange device: reflector-limited bound (every return lands on the auxiliary reflector; series with R_aux) and the
finite-geometry dependence on pixel size and substrate thickness (Al cathode)."""
import numpy as np, sys
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/orange']
import orange_reflectors as O


def mcg(args):
    seed, N, half, T, PS, RS, RB, WL = args
    rng = np.random.default_rng(seed)
    j = rng.choice(len(WL), N, p=WL); Pc = np.cumsum(PS / PS.sum(1, keepdims=True), 1)
    ib = (Pc[j] < rng.random(N)[:, None]).sum(1).clip(0, 89)
    x = (rng.random(N) - 0.5) * 2 * half; y = (rng.random(N) - 0.5) * 2 * half; ph = rng.random(N) * 2 * np.pi; w = np.ones(N)
    top = 0.0; alive = np.ones(N, bool)
    for _ in range(800):
        a = np.flatnonzero(alive)
        if a.size == 0: break
        for h in (0, 1):
            th = np.radians(O.TH[ib[a]] + rng.random(a.size) - 0.5); dd = T * np.tan(th)
            x[a] += dd * np.cos(ph[a]); y[a] += dd * np.sin(ph[a])
            out = (np.abs(x[a]) > 12) | (np.abs(y[a]) > 12); w[a[out]] = 0
            if h == 0:
                bt = O.BT[ib[a]]; top += (w[a] * bt).sum(); w[a] *= 1 - bt; u = rng.random(a.size)
                ib[a] = (O.BRc[:, ib[a]].T < u[:, None]).sum(1).clip(0, 89); ph[a] = rng.random(a.size) * 2 * np.pi
        inpix = (np.abs(x[a]) <= half) & (np.abs(y[a]) <= half)
        w[a] *= np.where(inpix, RS[j[a], ib[a]], RB[j[a], ib[a]]); alive &= w > 1e-5
    return top / N


if __name__ == '__main__':
    wp = (O.LAM * O.SPEC)[O.LAMJ]; wp = wp / wp.sum()
    for cat in ('Al', 'Ag'):
        with Pool(4) as p: res = p.map(O.device, [(cat, j) for j in O.LAMJ])
        ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
        es = float(wp @ ES); WL = wp * ES / (wp @ ES)
        for kind in ('DBR', 'Ag'):
            RB = np.array([O.back_R(j, kind) for j in O.LAMJ])
            ext = float(sum(WL[k] * O.series.eta_ext(O.BT, O.BR, RB[k], PS[k], 400)[0] for k in range(len(O.LAMJ))))
            print('%s cathode, every return on %s: eta_ext %.3f  EQE %.3f (PLQY 1: %.3f)' % (cat, kind, ext, O.PLQY * es * ext, es * ext), flush=True)
        if cat == 'Al':
            RB = np.array([O.back_R(j, 'DBR') for j in O.LAMJ])
            for half, T in ((1.0, 0.5), (1.0, 0.7), (1.0, 1.1), (0.5, 0.5), (0.5, 1.1), (1.5, 0.7)):
                with Pool(4) as p: e = np.mean(p.map(mcg, [(10 * k + 1, 120000, half, T, PS, RS, RB, WL) for k in range(8)]))
                print('   Al + DBR, pixel %.0f x %.0f mm, substrate %.1f mm: eta_ext %.3f  EQE %.3f' % (2 * half, 2 * half, T, e, O.PLQY * es * e), flush=True)

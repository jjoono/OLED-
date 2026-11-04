"""Monte Carlo transfer matrices of a non-absorbing scattering slab index-matched on both sides (n_slab = n_sub).
Henyey-Greenstein phase function g, optical thickness S = mu_s d.  Photons enter the bottom at polar angle bin j
(1-deg bins, centres 0.5..89.5 deg, uniform in azimuth), leave the top (T) or the bottom (R); exit polar angle binned.
Output npz: T[k, j], R[k, j] = fraction of the power entering in bin j that leaves in bin k (top / bottom).
Because the slab is symmetric the same matrices describe light entering from the top.
Usage: python3 slab_mc.py S g N_per_bin out.npz"""
import sys, numpy as np
def hg_cos(g, u):
    if abs(g) < 1e-6: return 2 * u - 1
    return (1 + g * g - ((1 - g * g) / (1 - g + 2 * g * u)) ** 2) / (2 * g)
def run(S, g, N, seed=1):
    rng = np.random.default_rng(seed); T = np.zeros((90, 90)); R = np.zeros((90, 90))
    for a in range(90):
        th = np.radians(a + 0.5); ph = rng.uniform(0, 2 * np.pi, N); st = np.sin(th)
        d = np.stack([st * np.cos(ph), st * np.sin(ph), np.full(N, np.cos(th))], 1); z = np.zeros(N); alive = np.ones(N, bool)
        while alive.any():
            idx = np.flatnonzero(alive)
            step = -np.log(rng.random(idx.size)) / S if S > 0 else np.full(idx.size, 1e9)
            zn = z[idx] + step * d[idx, 2]
            top = zn >= 1; bot = zn <= 0; mid = ~(top | bot)
            if top.any():
                ii = idx[top]; k = np.minimum(np.degrees(np.arccos(np.clip(d[ii, 2], 0, 1))).astype(int), 89)
                np.add.at(T[:, a], k, 1); alive[ii] = False
            if bot.any():
                ii = idx[bot]; k = np.minimum(np.degrees(np.arccos(np.clip(-d[ii, 2], 0, 1))).astype(int), 89)
                np.add.at(R[:, a], k, 1); alive[ii] = False
            if mid.any():
                ii = idx[mid]; z[ii] = zn[mid]
                ct = hg_cos(g, rng.random(ii.size)); stt = np.sqrt(np.clip(1 - ct ** 2, 0, 1)); p2 = rng.uniform(0, 2 * np.pi, ii.size)
                ux, uy, uz = d[ii, 0], d[ii, 1], d[ii, 2]; den = np.sqrt(np.clip(1 - uz ** 2, 1e-12, None))
                nx = stt * (ux * uz * np.cos(p2) - uy * np.sin(p2)) / den + ux * ct
                ny = stt * (uy * uz * np.cos(p2) + ux * np.sin(p2)) / den + uy * ct
                nz = -stt * np.cos(p2) * den + uz * ct
                v = den < 1e-6; nx[v] = stt[v] * np.cos(p2[v]); ny[v] = stt[v] * np.sin(p2[v]); nz[v] = np.sign(uz[v]) * ct[v]
                nn = np.sqrt(nx ** 2 + ny ** 2 + nz ** 2); d[ii] = np.stack([nx / nn, ny / nn, nz / nn], 1)
    return T / N, R / N
if __name__ == '__main__':
    S, g, N, out = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    T, R = run(S, g, N); np.savez(out, T=T, R=R, S=S, g=g)
    print(out, 'col sums %.4f..%.4f' % ((T + R).sum(0).min(), (T + R).sum(0).max()), 'T(normal) %.3f' % T[:, 0].sum(), flush=True)

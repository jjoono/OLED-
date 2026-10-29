"""Monte Carlo redistribution matrix of a scattering layer on an Ag reflector (the 'back diffuser').
Photons enter the layer (index n_s, index-matched top) at polar angle theta_i, Henyey-Greenstein scattering
(anisotropy g, optical thickness S = mu_s d), absorption mu_a d (Beer-Lambert along the path), specular reflection
at the Ag (R_Ag(theta), s/p average, McPeak Ag at 530 nm seen from n_s); they leave through the top at theta_o.
B[o, i] = fraction leaving into bin o (2-deg bins of the polar angle in the layer) per photon entering in bin i.
-> scat_ag_B.npz with keys S{S}_g{g}_a{mua}"""
import numpy as np, sys
sys.path[:0] = ['/home/user/OLED-/sim/diag']
NB = 45; TB = np.arange(NB) * 2 + 1.0
def hg_cos(g, u):
    if abs(g) < 1e-6: return 2 * u - 1
    return (1 + g * g - ((1 - g * g) / (1 - g + 2 * g * u)) ** 2) / (2 * g)
def R_ag(ns=1.8, lam=530.0, dag=0.0):
    import dev3
    from front_lum import tmm
    n = dev3.D.AG[int(lam - 400)] + dag; out = []
    for th in np.linspace(0, 90, 181):
        kx = 2 * np.pi / lam * ns * np.sin(np.radians(min(th, 89.9)))
        out.append(0.5 * (tmm([n], [100.0], ns, 1.0, lam, kx, 's')[0] + tmm([n], [100.0], ns, 1.0, lam, kx, 'p')[0]))
    return np.array(out)
def redist(S, g, mua, Rag, N=40000, seed=1):
    rng = np.random.default_rng(seed); B = np.zeros((NB, NB))
    for i in range(NB):
        th = np.radians(TB[i] + rng.uniform(-1, 1, N)); ph = rng.uniform(0, 2 * np.pi, N); st = np.sin(th)
        d = np.stack([st * np.cos(ph), st * np.sin(ph), -np.cos(th)], 1)        # downwards
        z = np.ones(N); w = np.ones(N); alive = np.ones(N, bool); acc = np.zeros(NB)
        for it in range(100000):
            idx = np.flatnonzero(alive)
            if idx.size == 0: break
            dz = d[idx, 2]
            step = -np.log(rng.random(idx.size)) / S if S > 0 else np.full(idx.size, 1e9)
            # distance to the boundary along the direction
            tb = np.where(dz < 0, z[idx] / np.maximum(-dz, 1e-12), (1 - z[idx]) / np.maximum(dz, 1e-12))
            hit = step >= tb; s = np.where(hit, tb, step)
            w[idx] *= np.exp(-mua * s); z[idx] += s * dz
            hb = idx[hit & (dz < 0)]; ht = idx[hit & (dz >= 0)]; sc = idx[~hit]
            if hb.size:                                   # Ag: specular, weight x R(theta)
                ang = np.degrees(np.arccos(np.clip(-d[hb, 2], 0, 1)))
                w[hb] *= np.interp(ang, np.linspace(0, 90, 181), Rag); d[hb, 2] *= -1; z[hb] = 0.0
            if ht.size:                                   # leaves through the top
                ang = np.degrees(np.arccos(np.clip(d[ht, 2], 0, 1)))
                np.add.at(acc, np.minimum((ang / 2).astype(int), NB - 1), w[ht]); alive[ht] = False
            if sc.size:
                ct = hg_cos(g, rng.random(sc.size)); stt = np.sqrt(np.clip(1 - ct ** 2, 0, 1)); p2 = rng.uniform(0, 2 * np.pi, sc.size)
                dd = d[sc]; uz = dd[:, 2]; ux, uy = dd[:, 0], dd[:, 1]; den = np.sqrt(np.clip(1 - uz ** 2, 1e-12, None))
                nx = stt * (ux * uz * np.cos(p2) - uy * np.sin(p2)) / den + ux * ct
                ny = stt * (uy * uz * np.cos(p2) + ux * np.sin(p2)) / den + uy * ct
                nz = -stt * np.cos(p2) * den + uz * ct
                v = den < 1e-6; nx[v] = stt[v] * np.cos(p2[v]); ny[v] = stt[v] * np.sin(p2[v]); nz[v] = np.sign(uz[v]) * ct[v]
                d[sc] = np.stack([nx, ny, nz], 1)
            dead = idx[w[idx] < 1e-4]; alive[dead] = False
        B[:, i] = acc / N
    return B
if __name__ == '__main__':
    from multiprocessing import Pool
    Rag = R_ag()
    grid = [(S, g, a) for a in (0.0, 0.01, 0.03) for S in (0.0, 0.3, 1.0, 3.0, 10.0, 30.0) for g in (0.0, 0.6, 0.9) if not (S == 0 and g > 0)]
    with Pool(4) as p:
        Bs = p.starmap(redist, [(S, g, a, Rag) for S, g, a in grid])
    np.savez('scat_ag_B.npz', **{'S%g_g%g_a%g' % k: B for k, B in zip(grid, Bs)}, Rag=Rag)
    for k, B in zip(grid, Bs): print(k, 'mean reflectance (cos-sin weighted) %.3f' % (B.sum(0) @ (np.cos(np.radians(TB)) * np.sin(np.radians(TB))) / (np.cos(np.radians(TB)) * np.sin(np.radians(TB))).sum()))

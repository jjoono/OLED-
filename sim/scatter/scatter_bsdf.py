"""Monte Carlo BSDF of a non-absorbing scattering layer index-matched to the substrate (n), planar top to air.
Henyey-Greenstein phase function (g), scattering optical thickness S = mu_s d. Photons enter from the substrate at
polar angle theta_i (1-deg bins, 0.5..89.5); they leave either through the top (Fresnel, unpolarised; recorded as air
angle -> B_T) or back through the bottom into the substrate (index-matched, recorded as substrate angle -> B_R).
Output: BSDF(180, 90) in the author's convention: rows 0-89 transmitted (air angle), rows 90-179 reflected
stored in reverse (row 179 - k = substrate angle k)."""
import numpy as np, sys
def fresnel_R(n, ci):
    si2 = (n**2) * (1 - ci**2); out = np.ones_like(ci); ok = si2 < 1
    ct = np.sqrt(1 - si2[ok]); c = ci[ok]
    rs = ((n*c - ct)/(n*c + ct))**2; rp = ((n*ct - c)/(n*ct + c))**2
    out[ok] = 0.5*(rs + rp); return out
def hg_cos(g, u):
    if abs(g) < 1e-6: return 2*u - 1
    return (1 + g*g - ((1 - g*g)/(1 - g + 2*g*u))**2)/(2*g)
def bsdf(n, S, g, N=20000, seed=1):
    rng = np.random.default_rng(seed); B = np.zeros((180, 90))
    for a in range(90):
        th = np.radians(a + 0.5); mu = np.full(N, np.cos(th)); ph = rng.uniform(0, 2*np.pi, N)
        st = np.sin(th); d = np.stack([st*np.cos(ph), st*np.sin(ph), mu], 1)
        z = np.zeros(N); alive = np.ones(N, bool); w_T = np.zeros(90); w_R = np.zeros(90)
        for it in range(20000):
            idx = np.where(alive)[0]
            if idx.size == 0: break
            dz = d[idx, 2]; step = -np.log(rng.random(idx.size))/S if S > 0 else np.full(idx.size, 1e9)
            znew = z[idx] + step*dz
            top = znew >= 1; bot = znew <= 0; mid = ~(top | bot)
            # bottom exit -> reflected into substrate
            if bot.any():
                ii = idx[bot]; ang = np.degrees(np.arccos(np.clip(-d[ii, 2], 0, 1)))
                np.add.at(w_R, np.minimum(ang.astype(int), 89), 1); alive[ii] = False
            if top.any():
                ii = idx[top]; ci = d[ii, 2]; R = fresnel_R(n, ci); refl = rng.random(ii.size) < R
                out = ii[~refl]; s_air = np.clip(n*np.sqrt(1 - d[out, 2]**2), 0, 1)
                np.add.at(w_T, np.minimum(np.degrees(np.arcsin(s_air)).astype(int), 89), 1); alive[out] = False
                back = ii[refl]; d[back, 2] *= -1; z[back] = 1.0
            if mid.any():
                ii = idx[mid]; z[ii] = znew[mid]
                ct = hg_cos(g, rng.random(ii.size)); stt = np.sqrt(np.clip(1 - ct**2, 0, 1)); p2 = rng.uniform(0, 2*np.pi, ii.size)
                dd = d[ii]; uz = dd[:, 2]; ux, uy = dd[:, 0], dd[:, 1]
                den = np.sqrt(np.clip(1 - uz**2, 1e-12, None))
                nx = stt*(ux*uz*np.cos(p2) - uy*np.sin(p2))/den + ux*ct
                ny = stt*(uy*uz*np.cos(p2) + ux*np.sin(p2))/den + uy*ct
                nz = -stt*np.cos(p2)*den + uz*ct
                vert = den < 1e-6
                nx[vert] = stt[vert]*np.cos(p2[vert]); ny[vert] = stt[vert]*np.sin(p2[vert]); nz[vert] = np.sign(uz[vert])*ct[vert]
                d[ii] = np.stack([nx, ny, nz], 1)
        B[:90, a] = w_T/N; B[90:, a] = (w_R/N)[::-1]
    return B
if __name__ == '__main__':
    n = float(sys.argv[1]); out = {}
    for S in (1, 3, 10):
        for g in (0.5, 0.9):
            out[f'S{S}_g{g}'] = bsdf(n, S, g); print(n, S, g, 'col sums', out[f'S{S}_g{g}'].sum(0).min().round(4), flush=True)
    np.savez(f'/home/user/OLED-/sim/scatter/scatter_bsdf_n{n}.npz', **out)

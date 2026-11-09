"""Free-thickness angular filters.  With a fully randomizing layer, every encounter with the filter
sees a Lambertian distribution, so the multi-pass series has an exact per-wavelength closed form:

    eta_cone(lam) = p_cone / (1 - r * rho),   p_cone = sum_bins WL Tf F_cone,
    r = sum_bins WL Rf,  rho = sum_bins WL R_LED            (Rs = 0, Ts = Lambertian projector)

which makes a layer-by-layer optimization (differential evolution over all thicknesses) cheap.
The result is then re-checked with the full matrix model, including the physical scattering slab."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sys, json, time, numpy as np
from scipy.optimize import differential_evolution
import sel, common


def tmm_grid(nH, nL, ds, lam, n_in, th):
    """Unpolarized R, T of alternating H/L layers (thicknesses ds, first layer H) on (lam, theta).
    nH, nL: complex arrays over lam."""
    kx = (n_in * np.sin(np.radians(th)))[None, :].astype(complex)          # (1, A)
    lam = lam[:, None]
    R = 0.0; T = 0.0
    for pol in ('s', 'p'):
        def eta(n):
            c = np.sqrt(n ** 2 - kx ** 2 + 0j); c = np.where(c.imag < 0, -c, c)
            return (c if pol == 's' else n ** 2 / c), c
        m00 = np.ones(kx.shape, complex) * np.ones(lam.shape); m01 = np.zeros_like(m00); m10 = np.zeros_like(m00); m11 = np.ones_like(m00)
        for k, d in enumerate(ds):
            n = (nH if k % 2 == 0 else nL)[:, None]
            e, c = eta(n); dl = 2 * np.pi / lam * c * d
            cs, sn = np.cos(dl), np.sin(dl)
            a00, a01, a10, a11 = cs, 1j * sn / e, 1j * e * sn, cs
            m00, m01, m10, m11 = m00 * a00 + m01 * a10, m00 * a01 + m01 * a11, m10 * a00 + m11 * a10, m10 * a01 + m11 * a11
        e0, _ = eta(np.complex128(n_in)); es, cs_ = eta(np.complex128(1.0))
        B = m00 + m01 * es; C = m10 + m11 * es
        r = (e0 * B - C) / (e0 * B + C); t = 2 * e0 / (e0 * B + C)
        Tp = np.real(es) / np.real(e0) * np.abs(t) ** 2
        Tp = np.where(np.abs(cs_.real) < 1e-9, 0.0, Tp)
        R = R + 0.5 * np.abs(r) ** 2; T = T + 0.5 * np.nan_to_num(Tp)
    return np.clip(R, 0, 1), np.clip(T, 0, 1)


class Problem:
    def __init__(self, stack, spec_name, beta, mat, key='cone', wmin=1e-3, coarse=False):
        DV = np.load('data/devices.npz')
        self.ns = float(DV[f'{stack}_n_sub']); self.spec = common.SPECTRA[spec_name]
        self.es = DV[f'{stack}_eta_sub']; self.rho = (DV[f'{stack}_R'] * sel.WL).sum(1)
        keep = self.spec > wmin * self.spec.max()                         # wavelengths that matter
        self.k = np.where(keep)[0]; self.lam = sel.LAM[self.k]
        H, L = sel.DBR_MATS[mat]
        self.nH = sel.nk(H)[self.k]; self.nL = sel.nk(L)[self.k]
        F = sel.cone_frac(self.ns)[common.BI[5.0 if key == 'fwd' else beta]]
        self.key, self.beta = key, beta
        if coarse:                                                        # 1-deg bins for the search
            self.th = sel.TH.reshape(-1, 2).mean(1); self.w = sel.WL.reshape(-1, 2).sum(1)
            self.F = (F * sel.WL).reshape(-1, 2).sum(1) / self.w
        else:
            self.th, self.w, self.F = sel.TH, sel.WL, F

    def eqe(self, ds):
        R, T = tmm_grid(self.nH, self.nL, ds, self.lam, self.ns, self.th)
        p = (T * self.w * self.F).sum(1); r = (R * self.w).sum(1)
        eta = p / (1 - r * self.rho[self.k])
        e = float((self.spec[self.k] * self.es[self.k] * eta).sum())
        return e / np.sin(np.radians(5.0)) ** 2 if self.key == 'fwd' else e


def edge_seq(mat, N, lb):
    """Quarter-wave edge filter (L/2 H L/2)^N written as an H-first alternating sequence."""
    lay = sel.edge_dbr(mat, N, lb); dL2, dH = lay[0][1], lay[1][1]
    return np.array([0.0, dL2] + [dH, 2 * dL2] * (N - 1) + [dH, dL2])


def optimize(stack, spec_name, beta, mat, N, key='cone', seed=1, maxiter=120, hi=400.0):
    """DE over all 2N+2 thicknesses, seeded with the best quarter-wave edge filter of the sweep."""
    SW = common.Sweep(); spec = common.SPECTRA[spec_name]
    d0, _ = SW.best_edge(stack, 'full', spec, beta, mat=mat, N=N, key=key)
    x0 = np.clip(edge_seq(mat, N, d0[3]), 0, hi); M = x0.size
    P = Problem(stack, spec_name, beta, mat, key, wmin=1e-2, coarse=True)
    rng = np.random.default_rng(seed); pop = 8 * M
    init = np.vstack([x0, np.clip(x0 * rng.uniform(0.8, 1.2, (pop // 2, M)), 0, hi), rng.uniform(0, hi, (pop - 1 - pop // 2, M))])
    t0 = time.time()
    res = differential_evolution(lambda x: -P.eqe(x), [(0.0, hi)] * M, seed=seed, maxiter=maxiter, init=init,
                                 tol=1e-8, mutation=(0.5, 1.0), recombination=0.7, polish=True)
    Pf = Problem(stack, spec_name, beta, mat, key)
    return dict(stack=stack, spec=spec_name, beta=beta, mat=mat, M=M, key=key, eqe=Pf.eqe(res.x), eqe_edge=Pf.eqe(x0),
                edge=str(d0), ds=[round(float(d), 2) for d in res.x], nfev=int(res.nfev), sec=round(time.time() - t0, 1))


if __name__ == '__main__':
    cases = json.loads(sys.argv[1])
    from multiprocessing import Pool
    with Pool(4) as p:
        out = p.starmap(optimize, [tuple(c) for c in cases])
    old = json.load(open('data/optfilter.json')) if os.path.exists('data/optfilter.json') else []
    json.dump(old + out, open('data/optfilter.json', 'w'), indent=1)
    for o in out:
        print(o['stack'], o['spec'], o['beta'], o['mat'], o['M'], o['key'], round(o['eqe'], 4), o['sec'], 's', flush=True)

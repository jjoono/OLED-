"""Finite square pixel (opening w) of the recycling OLED at 530 nm: lateral Monte Carlo of the recycling passes
(as psf_recycle.py) with emission uniform over the opening and the pixel edge either open (no wall) or a reflective
sidewall (reflectance R_w, specular, spanning the thin-film stack).  Records on-axis light (air angle <~ 11 deg) leaving
inside the own opening vs outside (crosstalk), relative to an infinitely large pixel.
Usage: python3 pixel_mc.py <structure json> <case F1|F2> <N>"""
import sys, json, numpy as np
import ce_fab as F, ce_opt as C
D = C.D; j = 130; lam = 530.0
def setup(fn, case, N, S=10.0, g=0.0, mua=0.01, d_s=1.0):
    s = json.load(open(fn))['structure']; x = list(s.values())
    P, Ld, Pair, (Ps, Tf, Rf, ns) = F.point((case, j, tuple(x), N, 3000))
    fr, back, lamc, _ = F.layers(case, j, x); jc = int(round(lamc - 400))
    lay = [(l[0].real if hasattr(l[0], 'real') else l[0], l[2]) for l in back] + [(D.EMo[j].real, 25.0)] + [(np.real(l[0]), l[2]) for l in fr]
    lay += [(C.ZNS[j].real, lamc / 4 / C.ZNS[jc].real), (C.LIF[j].real, lamc / 4 / C.LIF[jc].real)] * N
    lay = [(np.real(n), d) for n, d in lay if np.real(n) > 1.0]           # metal layers (thin Ag) do not carry propagating light
    TB = np.radians(C.TB)
    def up(th):
        sn = ns * np.sin(th); return sum(d * min(sn / n, 0.999) / np.sqrt(1 - min(sn / n, 0.999) ** 2) for n, d in lay)
    UP = np.array([up(t) for t in TB]) / 1000.0
    return dict(Ps=Ps, Tf=Tf, Rf=Rf, B=C.BZ['S%g_g%g_a%g' % (S, g, mua)], UP=UP, TS=np.tan(np.minimum(TB, np.radians(85))), d_s=d_s, Ld=Ld)
def run(M, w, Rw, n_ph=300000, seed=1):
    rng = np.random.default_rng(seed); Ps = M['Ps']; B = M['B']
    Bc = B / np.maximum(B.sum(0), 1e-12); Bcum = np.cumsum(Bc, 0); refl = B.sum(0)
    ib = rng.choice(45, n_ph, p=Ps / Ps.sum()); wgt = np.full(n_ph, Ps.sum())
    xu = (rng.random((n_ph, 2)) - 0.5) * (w if np.isfinite(w) else 0.0)            # unfolded position
    alive = np.ones(n_ph, bool); own = 0.0; out = 0.0
    def fold(v):                                                               # crossings of the walls at +-w/2
        return np.floor((v + w / 2) / w) if np.isfinite(w) else np.zeros_like(v)
    k0 = fold(xu)
    for it in range(400):
        a = np.flatnonzero(alive)
        if a.size == 0: break
        th = ib[a]; wgt[a] *= refl[th]
        ob = (Bcum[:, th].T < rng.random(a.size)[:, None]).sum(1).clip(0, 44)
        ph = rng.uniform(0, 2 * np.pi, (a.size, 3))
        xu[a] += M['d_s'] * (M['TS'][th][:, None] * np.c_[np.cos(ph[:, 0]), np.sin(ph[:, 0])] + M['TS'][ob][:, None] * np.c_[np.cos(ph[:, 1]), np.sin(ph[:, 1])])
        u = np.c_[np.cos(ph[:, 2]), np.sin(ph[:, 2])]; xu[a] += M['UP'][ob][:, None] * u
        if Rw is not None and np.isfinite(w):
            k = fold(xu[a]); wgt[a] *= Rw ** np.abs(k - k0[a]).sum(1); k0[a] = k
        on = ob < 3; e = wgt[a] * M['Tf'][ob] * on
        if Rw is None and np.isfinite(w):
            inside = (np.abs(xu[a]) <= w / 2).all(1); own += e[inside].sum(); out += e[~inside].sum()
        else: own += e.sum()
        wgt[a] *= M['Rf'][ob]; xu[a] += M['UP'][ob][:, None] * u
        if Rw is not None and np.isfinite(w):
            k = fold(xu[a]); wgt[a] *= Rw ** np.abs(k - k0[a]).sum(1); k0[a] = k
        ib[a] = ob; alive[a[wgt[a] < 1e-4]] = False
    return own / n_ph, out / n_ph
if __name__ == '__main__':
    fn, case, N = sys.argv[1], sys.argv[2], int(sys.argv[3]); M = setup(fn, case, N)
    inf, _ = run(M, np.inf, None)
    print('infinite pixel: recycled on-axis (arb.) %.4g' % inf)
    for w in (10, 20, 40, 80, 160):
        o_open, x_open = run(M, float(w), None)
        line = 'w %3d um | open edge: own %.2f of infinite, crosstalk %.0f %% of its on-axis light' % (w, o_open / inf, 100 * x_open / max(o_open + x_open, 1e-12))
        for Rw in (0.90, 0.97):
            o, _ = run(M, float(w), Rw); line += ' | wall R %.2f: own %.2f' % (Rw, o / inf)
        print(line, flush=True)

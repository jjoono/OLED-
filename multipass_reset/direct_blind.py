"""(1) Direct multi-pass ray tracing in (theta, phi) with a specular, angle-resolved device reflection, as an
independent check of the Markov-on-polar-angle model (which averages azimuth and position at every encounter).
(2) Blind test: structure families that were not used to develop the reduced model.
Output: blind/<name>.npz (BT, BR), blind/results.json."""
import os, sys, json, time, inspect, numpy as np
from numba import njit, prange
import hftrace as T, hftrace_gen as G

OUT = 'blind'; os.makedirs(OUT, exist_ok=True)
DV = np.load('../selectivity/data/devices.npz')


def make_dir_trace(src_fn, sig_old, sig_new, surf_name, surf):
    src = inspect.getsource(src_fn.py_func); src = src[src.index('def _trace'):]
    src = src.replace(sig_old, sig_new)
    src = src.replace('return 1, np.arccos(min(1.0, -dz))', 'return 1, dx, dy, dz')
    src = src.replace('return 0, np.arccos(min(1.0, dz))', 'return 0, dx, dy, dz')
    src = src.replace('return 2, 0.0', 'return 2, 0.0, 0.0, 0.0')
    ns = {'np': np, surf_name: surf}; exec(src, ns); return njit(ns['_trace'])


_dir_hex = make_dir_trace(T._trace, 'X', 'X', '_surf', T._surf)
_g_src_sig = 'def _trace(th, phi, sx, sy, n, htab, dtab, hmax, lip, max_ev):'
_g_new_sig = 'def _trace(th, phi, sx, sy, n, P0, P1, P2, P3, P4, P5, hmax, lip, max_ev):'


def _gen_dir(surf):
    src = inspect.getsource(T._trace.py_func); src = src[src.index('def _trace'):]
    src = src.replace(_g_src_sig, _g_new_sig)
    for a in ('px, py', 'qx, qy', 'px + mid * dx, py + mid * dy'):
        src = src.replace(f'_surf({a}, htab, dtab)', f'_SURF({a}, P0, P1, P2, P3, P4, P5)')
    src = src.replace('return 1, np.arccos(min(1.0, -dz))', 'return 1, dx, dy, dz')
    src = src.replace('return 0, np.arccos(min(1.0, dz))', 'return 0, dx, dy, dz')
    src = src.replace('return 2, 0.0', 'return 2, 0.0, 0.0, 0.0')
    ns = {'np': np, '_SURF': surf}; exec(src, ns); return njit(ns['_trace'])


_dir_lens = _gen_dir(G._surf_lens); _dir_wave = _gen_dir(G._surf_wave)


def _make_mp(trace):
    @njit(parallel=True)
    def mp(N, n, P0, P1, P2, P3, P4, P5, hmax, lip, max_ev, seed, bx1, by1, bx2, by2, cdf, R, randphi):
        esc = np.zeros(N); air = np.zeros((N,), np.int64); npass = np.zeros(N, np.int64)
        for i in prange(N):
            np.random.seed(seed + i)
            u = np.random.random(); k = np.searchsorted(cdf, u)
            s2a = np.sin(np.radians(k)) ** 2; s2b = np.sin(np.radians(k + 1)) ** 2
            th = np.arcsin(np.sqrt(s2a + (s2b - s2a) * np.random.random()))
            phi = 2 * np.pi * np.random.random()
            air[i] = -1
            for ps in range(400):
                a = np.random.random(); b = np.random.random()
                sx = a * bx1 + b * bx2; sy = a * by1 + b * by2
                oc, dx, dy, dz = trace(th, phi, sx, sy, n, P0, P1, P2, P3, P4, P5, hmax, lip, max_ev)
                npass[i] = ps + 1
                if oc == 0:
                    esc[i] = 1.0; air[i] = min(int(np.degrees(np.arccos(min(1.0, dz)))), 89); break
                if oc == 2:
                    break
                th = np.arccos(min(1.0, -dz)); kb = min(int(np.degrees(th)), 89)
                if np.random.random() >= R[kb]:
                    break
                phi = np.arctan2(dy, dx) if randphi == 0 else 2 * np.pi * np.random.random()
        return esc, air, npass
    return mp


_mp_lens = _make_mp(_dir_lens); _mp_wave = _make_mp(_dir_wave)


@njit
def _hexwrap(th, phi, sx, sy, n, H, D, P2, P3, P4, P5, hmax, lip, max_ev):
    return _dir_hex(th, phi, sx, sy, n, H, D, hmax, lip, max_ev)


_mp_hex = _make_mp(_hexwrap)


def direct(kind, P, box, hmax, lip, dev='Al', i=75, N=1_000_000, randphi=0, seed=3):
    Pv = DV[f'{dev}_P'][i].reshape(90, 2).sum(1); R = DV[f'{dev}_R'][i].reshape(90, 2).mean(1)
    cdf = np.cumsum(Pv) / Pv.sum()
    fn = {'hex': _mp_hex, 'lens': _mp_lens, 'wave': _mp_wave}[kind]
    esc, air, npass = fn(N, 1.51, *P, hmax, lip, 80, seed, *box, cdf, R, randphi)
    eta = esc.mean() * Pv.sum()
    S = np.bincount(air[air >= 0], minlength=90)
    return float(eta), float(esc.std() / np.sqrt(N) * Pv.sum()), float(S[40:60].sum() / max(S.sum(), 1)), float(npass.mean())


# ------------------------------------------------------------------ geometries
NT = T.NT; r = np.linspace(0, 1, NT)


def hex_params(h):
    htab, dtab = T.tables(r, h); P = (htab, dtab, np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1))
    return P, (2.0, 0.0, 1.0, np.sqrt(3.0)), float(htab.max()) + 1e-6, float(np.abs(dtab).max())


def lens_params(H, S):
    H = np.ascontiguousarray(np.atleast_2d(H), float); D = np.ascontiguousarray(np.gradient(H, r, axis=1))
    S = np.ascontiguousarray(S, float); K = S.shape[0]
    hmax = float((H.max(1)[S[..., 3].astype(int)] * S[..., 2]).max()) + 1e-6
    reach = float((S[..., 2] + np.hypot(S[..., 0], S[..., 1])).max()); ring = 1 if reach <= 1.6 else 2
    P = (H, D, S, np.zeros(K), np.array([float(ring)]), np.zeros(1))
    return P, (K * G.A, 0.0, K * G.A / 2.0, K * G.A * G.SQ3 / 2.0), hmax, float(np.abs(D).max())


def wave_params(amp, gx, gy, ph, L):
    h0 = float(np.abs(amp).sum()) + 1e-6
    xs = np.linspace(0, L, 400); X, Y = np.meshgrid(xs, xs)
    a = gx[:, None, None] * X + gy[:, None, None] * Y + ph[:, None, None]
    hx = -(amp[:, None, None] * gx[:, None, None] * np.sin(a)).sum(0); hy = -(amp[:, None, None] * gy[:, None, None] * np.sin(a)).sum(0)
    return (amp, gx, gy, ph, np.array([h0]), np.zeros(1)), (L, 0.0, 0.0, L), 2 * h0, 1.3 * float(np.hypot(hx, hy).max())


def ring_wave(rms_slope, k0=2.0, dk=0.25, L=12.0, seed=0):
    """Quasi-periodic surface: spectrum concentrated on a ring |k| = k0 (not Gaussian-correlated)."""
    rng = np.random.default_rng(seed); g = []
    M = int(np.ceil((k0 + 3 * dk) * L / (2 * np.pi)))
    for m in range(0, M + 1):
        for k in range(-M, M + 1):
            if m == 0 and k <= 0:
                continue
            gx, gy = 2 * np.pi * m / L, 2 * np.pi * k / L
            if abs(np.hypot(gx, gy) - k0) < 3 * dk:
                g.append((gx, gy))
    g = np.array(g); kk = np.hypot(g[:, 0], g[:, 1])
    amp = np.exp(-((kk - k0) / dk) ** 2 / 2) * rng.rayleigh(1.0, len(g)); ph = rng.uniform(0, 2 * np.pi, len(g))
    amp *= rms_slope / np.sqrt(0.5 * np.sum((amp * kk) ** 2))
    return amp, g[:, 0], g[:, 1], ph, L


def cap(alpha_deg):
    a = np.radians(alpha_deg); Rs = 1 / np.sin(a)
    return np.sqrt(np.clip(Rs ** 2 - r ** 2, 0, None)) - Rs * np.cos(a)


def bsdf(kind, P, box, hmax, lip, N):
    if kind == 'hex':
        BT = np.zeros((90, 90)); BR = np.zeros((90, 90)); lost = np.zeros(90); e = np.linspace(0, 90, 91)
        for i in range(90):
            oc, to = T._batch(np.radians(i + .5), N, 1.51, P[0], P[1], hmax, lip, 80, 7919 + i * 100003)
            d = np.degrees(to); BT[i] = np.histogram(d[oc == 0], e)[0] / N; BR[i] = np.histogram(d[oc == 1], e)[0] / N; lost[i] = np.mean(oc == 2)
        return BT, BR, lost
    batch = G._batch_lens if kind == 'lens' else G._batch_wave
    return G._run(batch, P, box, hmax, lip, 1.51, N, 80, 1)


def structures():
    rng = np.random.default_rng(2026); hemi = np.sqrt(np.clip(1 - r ** 2, 0, None)); out = {}
    MB = np.load('bsdfs.npz')
    # direct-trace checks on library geometries
    from scipy.interpolate import PchipInterpolator
    def prof(xv):
        px = np.concatenate(([0.0], np.sort(xv[0:5]), [1.0])); py = np.concatenate(([1.0], xv[5:10], [0.0])); px, u = np.unique(px, return_index=True)
        return PchipInterpolator(px, py[u])(r) * xv[12]
    out['hemi_total'] = ('hex', *hex_params(prof(MB['hemi_total_x'])), 'library')
    out['free_40-60'] = ('hex', *hex_params(prof(MB['free_40-60_x'])), 'library')
    out['W_slope0.2'] = ('wave', *wave_params(*G.wave_surface(0.2, seed=0)), 'library')
    out['W_slope1.0'] = ('wave', *wave_params(*G.wave_surface(1.0, seed=0)), 'library')
    # blind families
    out['cap60'] = ('hex', *hex_params(cap(60.0)), 'blind')
    out['cap35'] = ('hex', *hex_params(cap(35.0)), 'blind')
    out['paraboloid'] = ('hex', *hex_params(1.0 - r ** 2), 'blind')
    out['cone45'] = ('hex', *hex_params(1.0 - r), 'blind')
    out['cone60'] = ('hex', *hex_params(np.sqrt(3.0) * (1.0 - r)), 'blind')
    g = np.exp(-r ** 2 / 0.18); out['gauss_bump'] = ('hex', *hex_params(1.5 * (g - g[-1]) / (1 - g[-1])), 'blind')
    K = 6; S = np.zeros((K, K, 4)); S[..., 0:2] = rng.uniform(-1, 1, (K, K, 2)) * 0.15 * G.A
    S[..., 2] = 1 + 0.2 * rng.uniform(-1, 1, (K, K)); S[..., 3] = rng.integers(0, 2, (K, K))
    out['mixed_cone_hemi'] = ('lens', *lens_params(np.stack([hemi, np.sqrt(3.0) * (1 - r)]), S), 'blind')
    out['ring_wave0.6'] = ('wave', *wave_params(*ring_wave(0.6, seed=4)), 'blind')
    out['ring_wave1.2'] = ('wave', *wave_params(*ring_wave(1.2, seed=5)), 'blind')
    return out


if __name__ == '__main__':
    sys.path.insert(0, '../selectivity'); os.chdir('../selectivity'); import mixing; os.chdir('../multipass_reset')
    W1 = np.diff(np.sin(np.radians(np.linspace(0, 90, 91))) ** 2)
    res = json.load(open(f'{OUT}/results.json')) if os.path.exists(f'{OUT}/results.json') else {}
    for name, (kind, P, box, hmax, lip, tag) in structures().items():
        if name in res and ('direct' in res[name] or name in ('cap60', 'paraboloid', 'cone45', 'gauss_bump', 'ring_wave1.2')):
            continue
        t0 = time.time(); f = f'{OUT}/{name}.npz'
        if tag == 'library':
            src = {'hemi_total': 'bsdfs.npz', 'free_40-60': 'bsdfs.npz'}.get(name)
            if src:
                z = np.load(src); BT, BR = z[name + '_BT'], z[name + '_BR']
            else:
                z = np.load(f'random_sweep/{name}.npz'); BT, BR = z['BT'], z['BR']
        elif os.path.exists(f):
            z = np.load(f); BT, BR = z['BT'], z['BR']
        else:
            BT, BR, lost = bsdf(kind, P, box, hmax, lip, 20000 if kind != 'wave' else 10000)
            np.savez(f, BT=BT, BR=BR, lost=lost)
        BRc = BR.T; pT = BT.sum(1); lam2 = mixing.lam2_cut(BRc); pinf, pt = mixing.pinf_cut(BRc)
        r_ = dict(tag=tag, kind=kind, p=float(pT @ (W1 / W1.sum())), lam2=lam2, N_mix=(1 / (1 - lam2) if lam2 < 1 - 1e-9 else float('inf')), p_inf=pinf, xi=pinf / pt)
        for dev in ('Al', 'Ag', 'LL15'):
            for i, lam in ((50, 500), (75, 550), (100, 600)):
                Pv = DV[f'{dev}_P'][i].reshape(90, 2).sum(1); R = DV[f'{dev}_R'][i].reshape(90, 2).mean(1)
                r_[f'{dev}_{lam}'] = dict(full=mixing.full(pT, BRc, R, Pv), red=mixing.reduced(pT, W1, R, Pv, max(1 - lam2, 0.0)),
                                          red_w=mixing.reduced(pT, W1, R, Pv, max(1 - lam2, 0.0), recip=False),
                                          memoryless=mixing.closed(pT, W1, R), m1=mixing.reduced(pT, W1, R, Pv, 1.0))
        if kind in ('hex', 'wave', 'lens') and name in ('hemi_total', 'free_40-60', 'W_slope0.2', 'W_slope1.0', 'cap35', 'cone60', 'mixed_cone_hemi', 'ring_wave0.6'):
            d0 = direct(kind, P, box, hmax, lip); d1 = direct(kind, P, box, hmax, lip, randphi=1, seed=5)
            x = np.linalg.solve(np.eye(90) - (DV['Al_R'][75].reshape(90, 2).mean(1))[:, None] * BRc, DV['Al_P'][75].reshape(90, 2).sum(1))
            air = x @ BT
            r_['direct'] = dict(eta=d0[0], se=d0[1], S4060=d0[2], passes=d0[3], eta_randphi=d1[0], se_randphi=d1[1],
                                S4060_model=float(air[40:60].sum() / air.sum()))
        res[name] = r_; json.dump(res, open(f'{OUT}/results.json', 'w'), indent=1)
        a = r_['Al_550']
        print(f"{name:16s} [{tag}] p {r_['p']:.3f} N_mix {r_['N_mix']:.2f} xi {r_['xi']:.2f} | Al550 full {a['full']:.4f} red {a['red']:.4f} q~w {a['red_w']:.4f} memless {a['memoryless']:.4f}"
              + (f" | direct {r_['direct']['eta']:.4f}±{r_['direct']['se']:.4f} randphi {r_['direct']['eta_randphi']:.4f}" if 'direct' in r_ else '') + f"  [{time.time() - t0:.0f}s]", flush=True)

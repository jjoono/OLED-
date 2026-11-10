"""Randomness sweep of passive surfaces: p_Lamb, lambda2, N_mix and multi-pass extraction (Al, Ag, low-loss, 550 nm).
  A  shape spread s: spline control points of every lenslet = quarter-circle points + s * noise (6 x 6 supercell)
  B  aspect spread : hemispheres with heights scaled by 1 + s * U(-1, 1)
  C  position jitter d (fraction of the pitch) with hemispheres; C2 adds radius spread 20 %
  W  Gaussian-correlated wavy surface, rms slope sigma
Each structure is saved to random_sweep/<name>.npz; a summary goes to random_sweep/summary.json."""
import os, sys, json, time, numpy as np
import hftrace_gen as G, hftrace_rand as T

OUT = 'random_sweep'; os.makedirs(OUT, exist_ok=True)
K = 6; N = 20000
th = np.arange(90) + .5; W = np.cos(np.radians(th)) * np.sin(np.radians(th)); W /= W.sum()
HX = np.array([0.3, 0.5, 0.7, 0.85, 0.95]); HY = np.sqrt(1 - HX ** 2)
DV = np.load('../selectivity/data/devices.npz')


def metrics(BT, BR):
    M = BR.T / np.maximum(BR.sum(1), 1e-12)[None, :]
    ev = np.linalg.eigvals(M); ev = ev[np.argsort(-np.abs(ev))]
    lam2 = float(np.abs(ev[1])); out = dict(p=float(BT.sum(1) @ W), lam2=lam2, N_mix=1 / max(1 - lam2, 1e-9))
    for dev in ('Al', 'Ag', 'LL15'):
        P = DV[f'{dev}_P'][75].reshape(90, 2).sum(1); R = DV[f'{dev}_R'][75].reshape(90, 2).mean(1)
        v = P.copy(); air = np.zeros(90)
        for _ in range(600):
            air += v @ BT; v = (v @ BR) * R
            if v.sum() < 1e-9:
                break
        out[f'eta_{dev}'] = float(air.sum()); out[f'S4060_{dev}'] = float(air[40:60].sum() / air.sum())
    return out


def hemi_sites(rng, jitter=0.0, rspread=0.0, nprof=1):
    S = np.zeros((K, K, 4))
    S[..., 0:2] = rng.uniform(-1, 1, (K, K, 2)) * jitter * G.A
    S[..., 2] = 1.0 + rspread * rng.uniform(-1, 1, (K, K))
    S[..., 3] = rng.integers(0, nprof, (K, K)) if nprof > 1 else 0
    return S


def run(name, fn):
    f = f'{OUT}/{name}.npz'
    if os.path.exists(f):
        z = np.load(f); return metrics(z['BT'], z['BR'])
    t = time.time(); BT, BR, lost = fn()
    np.savez(f, BT=BT, BR=BR, lost=lost)
    m = metrics(BT, BR); m['lost_max'] = float(lost.max()); m['sec'] = round(time.time() - t)
    return m


if __name__ == '__main__':
    rng = np.random.default_rng(11); res = {}
    r = np.linspace(0, 1, T.NT)
    jobs = []
    for s in (0.0, 0.1, 0.2, 0.4, 0.7, 1.0):                                    # A: shape spread
        prof = []
        while len(prof) < K * K:
            xv = np.concatenate([np.clip(HX + s * rng.normal(0, 0.25, 5), 0.02, 0.98), np.clip(HY + s * rng.normal(0, 0.4, 5), 0, 1.5), [1.0]])
            h = np.clip(T.spline_profile(xv), 0, None)
            if np.abs(np.gradient(h, r)).max() < 60:
                prof.append(h)
        Hs = np.array(prof); S = hemi_sites(rng, nprof=K * K); S[..., 3] = rng.permutation(K * K).reshape(K, K)
        jobs.append((f'A_shape{s}', lambda Hs=Hs, S=S: G.bsdf_lens(Hs, S, N=N)))
    hemi = np.sqrt(np.clip(1 - r ** 2, 0, None))
    for s in (0.2, 0.4, 0.6):                                                   # B: aspect spread
        sc = 1 + s * rng.uniform(-1, 1, K * K); Hs = hemi[None, :] * sc[:, None]
        S = hemi_sites(rng, nprof=K * K); S[..., 3] = np.arange(K * K).reshape(K, K)
        jobs.append((f'B_aspect{s}', lambda Hs=Hs, S=S: G.bsdf_lens(Hs, S, N=N)))
    for d in (0.05, 0.1, 0.2, 0.3):                                             # C: position jitter
        S = hemi_sites(rng, jitter=d)
        jobs.append((f'C_jitter{d}', lambda S=S: G.bsdf_lens(hemi, S, N=N)))
        S2 = hemi_sites(rng, jitter=d, rspread=0.2)
        jobs.append((f'C2_jitter{d}_r20', lambda S2=S2: G.bsdf_lens(hemi, S2, N=N)))
    for sg in (0.2, 0.4, 0.7, 1.0, 1.5):                                         # W: wavy surface
        jobs.append((f'W_slope{sg}', lambda sg=sg: G.bsdf_wave(sg, N=N // 2)))
    for name, fn in jobs:
        res[name] = run(name, fn)
        json.dump(res, open(f'{OUT}/summary.json', 'w'), indent=1)
        m = res[name]
        print(f"{name:20s} p {m['p']:.4f}  lam2 {m['lam2']:.3f}  N_mix {m['N_mix']:5.1f} | eta Al {m['eta_Al']:.3f} Ag {m['eta_Ag']:.3f} LL15 {m['eta_LL15']:.3f} | "
              f"S40-60 {m['S4060_Al']:.3f}  [{m.get('sec', 0)}s, lost {m.get('lost_max', 0):.4f}]", flush=True)

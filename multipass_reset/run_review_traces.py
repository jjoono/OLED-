"""Extra traces for the technical review (items 9 and 16).
  size : random MLAs on a 10x10 supercell (s = 1.0 shape spread; 30% position jitter) and wavy surfaces on a 24x24 box
  rays : ray-count convergence for A_shape0.0 (N = 5k, 10k, 40k) and W_slope0.2 (N = 2.5k, 5k, 20k)
Output: review_traces/<name>.npz"""
import os, numpy as np
import hftrace_gen as G, hftrace_rand as T
OUT = 'review_traces'; os.makedirs(OUT, exist_ok=True)
HX = np.array([0.3, 0.5, 0.7, 0.85, 0.95]); HY = np.sqrt(1 - HX ** 2); r = np.linspace(0, 1, T.NT)


def profiles(s, n, rng):
    prof = []
    while len(prof) < n:
        xv = np.concatenate([np.clip(HX + s * rng.normal(0, 0.25, 5), 0.02, 0.98), np.clip(HY + s * rng.normal(0, 0.4, 5), 0, 1.5), [1.0]])
        h = np.clip(T.spline_profile(xv), 0, None)
        if np.abs(np.gradient(h, r)).max() < 60:
            prof.append(h)
    return np.array(prof)


def sites(K, rng, jitter=0.0, nprof=1, perm=False):
    S = np.zeros((K, K, 4)); S[..., 0:2] = rng.uniform(-1, 1, (K, K, 2)) * jitter * G.A; S[..., 2] = 1.0
    S[..., 3] = rng.permutation(K * K).reshape(K, K) if perm else 0
    return S


def save(name, fn):
    f = f'{OUT}/{name}.npz'
    if os.path.exists(f): return
    BT, BR, lost = fn(); np.savez(f, BT=BT, BR=BR, lost=lost); print(name, 'done', flush=True)


rng = np.random.default_rng(31)
hemi = np.sqrt(np.clip(1 - r ** 2, 0, None))
H0 = profiles(0.0, 1, rng)
save('A_shape0.0_N5000', lambda: G.bsdf_lens(H0, sites(6, rng), N=5000))
save('A_shape0.0_N10000', lambda: G.bsdf_lens(H0, sites(6, rng), N=10000))
save('A_shape0.0_N40000', lambda: G.bsdf_lens(H0, sites(6, rng), N=40000))
save('W_slope0.2_N2500', lambda: G.bsdf_wave(0.2, N=2500))
save('W_slope0.2_N5000', lambda: G.bsdf_wave(0.2, N=5000))
save('W_slope0.2_N20000', lambda: G.bsdf_wave(0.2, N=20000))
K = 10
H1 = profiles(1.0, K * K, rng)
save('A_shape1.0_K10', lambda: G.bsdf_lens(H1, sites(K, rng, nprof=K * K, perm=True), N=20000))
save('C_jitter0.3_K10', lambda: G.bsdf_lens(hemi, sites(K, rng, jitter=0.3), N=20000))
save('W_slope0.2_L24', lambda: G.bsdf_wave(0.2, N=10000, L=24.0))
save('W_slope1.0_L24', lambda: G.bsdf_wave(1.0, N=10000, L=24.0))

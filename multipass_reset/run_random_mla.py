"""Demo: BSDF of a hexagonal MLA whose lenslets carry independent random spline profiles (K x K supercell),
with escape probability p and angular mixing number N_mix = 1/(1 - lambda2)."""
import sys, time, numpy as np
import hftrace_rand as T

K = int(sys.argv[1]) if len(sys.argv) > 1 else 6
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 3)
prof = []
while len(prof) < K * K:
    xv = np.concatenate([np.sort(rng.uniform(0, 1, 5)), rng.uniform(0, 1.5, 5), [rng.uniform(0.3, 1.5)]])
    h = np.clip(T.spline_profile(xv), 0, None)
    if np.abs(np.gradient(h, np.linspace(0, 1, T.NT))).max() < 60:       # drop near-vertical walls (not a height field)
        prof.append(h)
H = np.array(prof); assign = rng.permutation(K * K).reshape(K, K)
t = time.time(); BT, BR, lost = T.bsdf(H, assign, n=1.51, N=20000)
th = np.arange(90) + .5; w = np.cos(np.radians(th)) * np.sin(np.radians(th)); w /= w.sum()
M = BR.T / np.maximum(BR.sum(1), 1e-12)[None, :]
ev = np.sort(np.real(np.linalg.eigvals(M)))[::-1]
print(f'K={K}: p_Lamb {BT.sum(1) @ w:.4f}  lambda2 {ev[1]:.3f}  N_mix {1 / (1 - ev[1]):.1f}  lost max {lost.max():.4f}  [{time.time() - t:.0f}s]')
np.savez(f'bsdf_random_K{K}.npz', BT=BT, BR=BR, lost=lost, H=H, assign=assign)

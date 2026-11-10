"""Angular mixing of passive outcoupling structures: definitions and validation.

Return operator of a structure: B_R[out, in] (substrate angle bins), transmission per encounter p(theta) = 1 - sum_out B_R.
M = B_R normalised per input column is the angular Markov kernel of the returned light.  For a passive reciprocal
structure, w_in B_R(in->out) = w_out B_R(out->in) with the etendue weight w = d(sin^2 theta) (detailed balance), so M is
reversible w.r.t. its stationary weight and its spectrum is real.

Definitions compared
  lam2      : second eigenvalue of M (relaxation time t_rel = 1/(1 - lam2), Levin-Peres-Wilmer, Ch. 12)
  lam_eff   : Rayleigh quotient of M for the deviation of the device's own emission from the stationary distribution,
              in the etendue-weighted inner product  <f, M f>_w / <f, f>_w,  f = v0/w - <v0/w>  (memory as experienced
              by the light that actually circulates)
  m_fit     : mixing fraction that makes the reduced model reproduce the full series
Reduced model: the structure is replaced by its angle-resolved transmission p(theta) and a single mixing fraction m;
returned light = (1 - p(theta)) [ (1 - m) same angle + m Lambertian ].  m = 1, scalar p: eq. (2) of the round-trip-loss paper.
"""
import json, numpy as np
from scipy.optimize import brentq
import sel

N_SUB = 1.51


def full(BT, BR, R, P0):
    x = np.linalg.solve(np.eye(len(P0)) - R[:, None] * BR, P0)
    return float(BT @ x)


def reduced(pT, w, R, P0, m, recip=True):
    """Reduced model: returned light = (1 - p_in) [ (1 - m) same angle + m q_out ].
    recip=True : q ∝ w (1 - p), the only rank-one mixing term compatible with reciprocity (detailed balance
                 w_i B(i->j) = w_j B(j->i)); light is mixed preferentially into the angles that also return light.
    recip=False: q ∝ w (Lambertian redistribution, violates reciprocity when p depends on angle; kept for comparison)."""
    q = w * (1 - pT) if recip else w.astype(float)
    q = q / q.sum(); n = len(P0)
    BR = (1 - pT)[None, :] * ((1 - m) * np.eye(n) + m * q[:, None])
    return full(pT, BR, R, P0)


def closed(pT, w, R):
    wn = w / w.sum(); p = float(pT @ wn); A = 1 - float(R @ wn)
    return p / (p + (1 - p) * A)


def lam2_cut(BR, cut_deg=85.0):
    """|lambda2| of the column-normalised return kernel restricted to theta < cut_deg (default: excludes 0.76 % of the
    etendue at grazing angles, where a few nearly unreachable states can return light at their own angle)."""
    nb = BR.shape[0]; th = (np.arange(nb) + 0.5) * 90.0 / nb; k = th < cut_deg
    B = BR[np.ix_(k, k)]; M = B / np.maximum(B.sum(0, keepdims=True), 1e-12)
    ev = np.linalg.eigvals(M); ev = ev[np.argsort(-np.abs(ev))]
    return float(np.abs(ev[1]))


def pinf_cut(BR, cut_deg=85.0):
    """Late-pass escape probability p_inf = 1 - rho(B_R) on theta < cut_deg: the per-encounter escape probability of
    the quasi-stationary distribution of light that has circulated many times (Perron root of the sub-stochastic
    return operator).  By reciprocity B_R is self-adjoint in the etendue-weighted inner product, so the Rayleigh
    quotient with the Lambertian vector gives rho >= 1 - p~ and hence  min p(theta) <= p_inf <= p~, where p~ is the
    uniform-radiance probability of leaving the truncated chain (escape plus leakage beyond cut_deg).  Returns (p_inf, p~)."""
    nb = BR.shape[0]; th = (np.arange(nb) + 0.5) * 90.0 / nb; k = th < cut_deg
    B = BR[np.ix_(k, k)]; rho = float(np.max(np.abs(np.linalg.eigvals(B))))
    w = np.diff(np.sin(np.radians(np.linspace(0, 90, nb + 1))) ** 2)[k]
    p = 1 - float((B.sum(0) * w).sum() / w.sum())          # escape plus leakage beyond the cutoff (truncated chain)
    return 1 - rho, p


def mixing_defs(BR, w, P0):
    cs = BR.sum(0); ok = cs > 1e-6
    M = np.where(ok[None, :], BR / np.maximum(cs, 1e-12)[None, :], 0)
    ev = np.sort(np.real(np.linalg.eigvals(M[np.ix_(ok, ok)])))[::-1]
    lam2 = float(ev[1])
    # stationary distribution and the deviation of the emitted light from it
    vals, vecs = np.linalg.eig(M[np.ix_(ok, ok)]); pi = np.abs(np.real(vecs[:, np.argmax(np.real(vals))])); pi /= pi.sum()
    v0 = P0[ok] / P0[ok].sum(); f = v0 / pi - 1.0                      # density ratio w.r.t. stationary
    Mf = M[np.ix_(ok, ok)].T @ f                                       # backward operator acting on functions
    lam_eff = float((pi * f * Mf).sum() / (pi * f * f).sum())
    return lam2, lam_eff


def m_fit(pT, w, R, P0, target):
    g = lambda m: reduced(pT, w, R, P0, m) - target
    a, b = g(0.0), g(1.0)
    if a * b > 0:
        return 0.0 if abs(a) < abs(b) else 1.0
    return brentq(g, 0.0, 1.0)


def structures():
    """Return {name: (BT_per_input, BR[out,in], weights, grid)} on their native grids."""
    SL = np.load('data/slabs.npz'); MB = np.load('../multipass_reset/bsdfs.npz'); out = {}
    Rf, Tf = sel.tmm_RT([], 550.0, sel.TH, N_SUB); I = np.eye(sel.NB)
    def slab(Ts, Rs):
        U = np.linalg.solve(I - Rs * Rf[None, :], Ts)
        return (Tf[:, None] * U).sum(0), Rs + (Ts * Rf[None, :]) @ U
    out['planar'] = slab(I, np.zeros_like(I)) + (sel.WL, '0.5')
    for k in ('g0.0_t0.25_w1', 'g0.0_t1.0_w1', 'g0.8_t0.25_w1', 'g0.8_t1.0_w1', 'g0.9_t1.0_w1', 'g0.8_t4.0_w1'):
        out['scatter ' + k.replace('_w1', '')] = slab(*SL[k]) + (sel.WL, '0.5')
    W1 = sel.WL.reshape(90, 2).sum(1)
    for d in sorted({k.rsplit('_', 1)[0] for k in MB.files if k.endswith('_BT')} - {'free_0-20'}):
        out['MLA ' + d] = (MB[d + '_BT'].sum(1), MB[d + '_BR'].T, W1, '1')
    return out


if __name__ == '__main__':
    DV = np.load('data/devices.npz'); ST = structures(); res = []
    for dev in ('Al', 'Ag', 'LL15'):
        for i, lam in ((50, 500.0), (75, 550.0), (100, 600.0)):
            P = DV[f'{dev}_P'][i]; R = DV[f'{dev}_R'][i]
            for name, (BT, BR, w, g) in ST.items():
                if g == '1':
                    P0 = P.reshape(90, 2).sum(1); Rg = R.reshape(90, 2).mean(1)
                else:
                    P0, Rg = P, R
                pT = BT if BT.ndim == 1 else BT
                eta = full(pT, BR, Rg, P0)
                lam2_full, lam_eff = mixing_defs(BR, w, P0); lam2 = lam2_cut(BR)
                r = dict(device=dev, lam=lam, structure=name, eta_full=eta, eta_closed=closed(pT, w, Rg),
                         eta_m0=reduced(pT, w, Rg, P0, 0.0), eta_m1=reduced(pT, w, Rg, P0, 1.0),
                         lam2=lam2, lam2_full=lam2_full, lam_eff=lam_eff, eta_red_lam2=reduced(pT, w, Rg, P0, 1 - lam2),
                         eta_red_eff=reduced(pT, w, Rg, P0, float(np.clip(1 - lam_eff, 0, 1))),
                         m_fit=m_fit(pT, w, Rg, P0, eta), p_lamb=float(pT @ (w / w.sum())))
                res.append(r)
    json.dump(res, open('data/mixing_validation.json', 'w'), indent=1)
    print(f"{'device':5s} {'lam':4s} {'structure':22s} {'full':>6s} | {'eq2':>6s} {'m=0':>6s} {'m=1':>6s} | {'lam2':>5s} {'red':>6s} | {'lam_eff':>7s} {'red':>6s} | m_fit")
    for r in res:
        if r['lam'] != 550.0: continue
        print(f"{r['device']:5s} {r['lam']:.0f} {r['structure']:22s} {r['eta_full']:.3f} | {r['eta_closed']:.3f} {r['eta_m0']:.3f} {r['eta_m1']:.3f} | "
              f"{r['lam2']:5.2f} {r['eta_red_lam2']:.3f} | {r['lam_eff']:7.2f} {r['eta_red_eff']:.3f} | {r['m_fit']:.2f}")
    err = lambda k: np.array([abs(r[k] - r['eta_full']) for r in res])
    for k in ('eta_closed', 'eta_m1', 'eta_red_lam2', 'eta_red_eff'):
        e = err(k); print(f'{k:14s} mean |err| {e.mean() * 100:.2f} pp, max {e.max() * 100:.2f} pp  (all devices, 3 wavelengths)')

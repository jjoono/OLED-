"""Polarization check of the scalar (polarization-averaged) model.
Scalar model: R = (Rs+Rp)/2, v0 = v0s + v0p, unpolarized Fresnel in the tracer.
Polarization-resolved limit without s/p conversion: each polarization keeps its label at every encounter
(exact for a planar interface; for curved surfaces the true result lies between this limit and full depolarization).
eta_pol = sum_pol p_pol (I - diag R_pol B_R,pol)^-1 v0_pol.  Al stack (61.3/97.5 nm), 500/550/600 nm.
Output: polcheck/results.json"""
import os, sys, json, inspect, numpy as np
from numba import njit
import hftrace as T, hftrace_gen as G
sys.path.insert(0, '../selectivity'); sys.path.insert(0, '../collimation/cps')
os.chdir('../selectivity'); import sel, cps2, mixing; os.chdir('../multipass_reset')
lam2_cut = mixing.lam2_cut
OUT = 'polcheck'; os.makedirs(OUT, exist_ok=True)
TH = np.arange(90) + 0.5; W1 = np.diff(np.sin(np.radians(np.linspace(0, 90, 91))) ** 2)

def pol_trace(src_fn, pol, gen=False, surf=None):
    src = inspect.getsource(src_fn.py_func); src = src[src.index('def _trace'):]
    a = 'Rf = 0.5 * (rs * rs + rp * rp)'; assert a in src
    src = src.replace(a, 'Rf = rs * rs' if pol == 's' else 'Rf = rp * rp')
    if gen:
        src = src.replace('def _trace(th, phi, sx, sy, n, htab, dtab, hmax, lip, max_ev):', 'def _trace(th, phi, sx, sy, n, P0, P1, P2, P3, P4, P5, hmax, lip, max_ev):')
        for a_ in ('px, py', 'qx, qy', 'px + mid * dx, py + mid * dy'):
            src = src.replace(f'_surf({a_}, htab, dtab)', f'_SURF({a_}, P0, P1, P2, P3, P4, P5)')
        ns = {'np': np, '_SURF': surf}
    else:
        ns = {'np': np, '_surf': T._surf}
    exec(src, ns); return njit(ns['_trace'])

def hex_bsdf(h, pol, N=20000):
    tr = pol_trace(T._trace, pol)
    from numba import prange
    @njit(parallel=True)
    def batch(th, N, n, htab, dtab, hmax, lip, max_ev, seed):
        oc = np.empty(N, np.int64); to = np.empty(N)
        for i in prange(N):
            np.random.seed(seed + i); u = np.random.random(); w = np.random.random()
            sx = u * 2.0 + w * 1.0; sy = w * np.sqrt(3.0); phi = 2 * np.pi * np.random.random()
            a, b = tr(th, phi, sx, sy, n, htab, dtab, hmax, lip, max_ev); oc[i] = a; to[i] = b
        return oc, to
    r = np.linspace(0, 1, T.NT); htab, dtab = T.tables(r, h); hmax = float(htab.max()) + 1e-6; lip = float(np.abs(dtab).max())
    e = np.linspace(0, 90, 91); BT = np.zeros((90, 90)); BR = np.zeros((90, 90))
    for i in range(90):
        oc, to = batch(np.radians(i + .5), N, 1.51, htab, dtab, hmax, lip, 60, 7919 + i * 100003)
        d = np.degrees(to); BT[i] = np.histogram(d[oc == 0], e)[0] / N; BR[i] = np.histogram(d[oc == 1], e)[0] / N
    return BT, BR

def wave_bsdf(slope, pol, N=10000):
    tr = pol_trace(T._trace, pol, gen=True, surf=G._surf_wave); batch = G._make_batch(tr)
    amp, gx, gy, ph, L = G.wave_surface(slope, seed=0); h0 = float(np.abs(amp).sum()) + 1e-6
    xs = np.linspace(0, L, 400); X, Y = np.meshgrid(xs, xs); a = gx[:, None, None] * X + gy[:, None, None] * Y + ph[:, None, None]
    hx = -(amp[:, None, None] * gx[:, None, None] * np.sin(a)).sum(0); hy = -(amp[:, None, None] * gy[:, None, None] * np.sin(a)).sum(0)
    lip = 1.3 * float(np.hypot(hx, hy).max())
    BT, BR, lost = G._run(batch, (amp, gx, gy, ph, np.array([h0]), np.zeros(1)), (L, 0.0, 0.0, L), 2 * h0, lip, 1.51, N, 80, 1)
    return BT, BR

def planar_bsdf(pol):
    n = 1.51; t = np.radians(TH); st = n * np.sin(t); BT = np.zeros((90, 90)); BR = np.zeros((90, 90))
    for i in range(90):
        if st[i] >= 1:
            BR[i, i] = 1.0; continue
        ct = np.sqrt(1 - st[i] ** 2); ci = np.cos(t[i])
        rs = (n * ci - ct) / (n * ci + ct); rp = (ct - n * ci) / (ct + n * ci); Rf = rs ** 2 if pol == 's' else rp ** 2
        BR[i, i] = Rf; BT[i, min(int(np.degrees(np.arcsin(st[i]))), 89)] = 1 - Rf
    return BT, BR

def device_pol(j):
    above, below, ns = sel.stack('Al', j); lam = float(sel.LAM[list(sel.IL).index(j)])
    S = cps2.Stack(lam, (sel.nk('EML_o', j).real, sel.nk('EML_e', j).real), 25.0, 12.5, above=above, below=below, n_sub=ns, h=0.865)
    th_f = np.linspace(0, 90, 1801)[:-1] + 0.025; out = {}
    k0 = 2 * np.pi / S.lam; kx = S.n_sub.real * np.sin(np.radians(TH)) * k0
    for pol in ('s', 'p'):
        pf = np.clip(np.real(cps2.sub_angular(S, th_f, pols=(pol,))), 0, None); out['P' + pol] = pf.reshape(90, 20).sum(1)
        r, _ = cps2.stack_rt(S.all_no, S.all_ne, S.all_d, kx, k0, pol); out['R' + pol] = np.clip(np.abs(r) ** 2, 0, 1)
    tot = out['Ps'].sum() + out['Pp'].sum(); out['Ps'] /= tot; out['Pp'] /= tot
    return out

def eta(BT, BR, R, P):
    return float((np.linalg.solve(np.eye(90) - R[:, None] * BR.T, P) @ BT).sum())

if __name__ == '__main__':
    MB = np.load('bsdfs.npz'); from scipy.interpolate import PchipInterpolator
    x = MB['hemi_total_x']; px = np.concatenate(([0.0], np.sort(x[0:5]), [1.0])); py = np.concatenate(([1.0], x[5:10], [0.0]))
    hprof = PchipInterpolator(px, py)(np.linspace(0, 1, T.NT)) * x[12]
    structs = {'planar': lambda pol: planar_bsdf(pol), 'hemi_total': lambda pol: hex_bsdf(hprof, pol), 'W_slope1.0': lambda pol: wave_bsdf(1.0, pol)}
    J = {500: 100, 550: 150, 600: 200}
    devs = {lam: device_pol(j) for lam, j in J.items()}
    res = {}
    for nm, fn in structs.items():
        f = f'{OUT}/{nm}.npz'
        if os.path.exists(f):
            z = np.load(f); B = {pol: (z['BT' + pol], z['BR' + pol]) for pol in ('s', 'p')}
        else:
            B = {pol: fn(pol) for pol in ('s', 'p')}; np.savez(f, BTs=B['s'][0], BRs=B['s'][1], BTp=B['p'][0], BRp=B['p'][1])
        BTa = 0.5 * (B['s'][0] + B['p'][0]); BRa = 0.5 * (B['s'][1] + B['p'][1])
        res[nm] = {}
        for lam, d in devs.items():
            e_sc = eta(BTa, BRa, 0.5 * (d['Rs'] + d['Rp']), d['Ps'] + d['Pp'])
            e_pol = eta(*B['s'], d['Rs'], d['Ps']) + eta(*B['p'], d['Rp'], d['Pp'])
            # complete reciprocal mixing and reduced model, scalar vs per polarization channel
            def mix_red(BT, BR, R, P):
                pT = BT.sum(1); l2 = lam2_cut(BR.T)
                return mixing.reduced(pT, W1, R, P, 1.0), mixing.reduced(pT, W1, R, P, max(1 - l2, 0.0))
            ms, rs_ = mix_red(BTa, BRa, 0.5 * (d['Rs'] + d['Rp']), d['Ps'] + d['Pp'])
            mp = [mix_red(*B[pl], d['R' + pl], d['P' + pl]) for pl in ('s', 'p')]
            e_mix_pol = sum(x[0] for x in mp); e_red_pol = sum(x[1] for x in mp)
            res[nm][lam] = dict(scalar=e_sc, pol_no_conversion=e_pol, diff_pp=(e_pol - e_sc) * 100, s_fraction=float(d['Ps'].sum()),
                                gain_mix_scalar_pp=(ms - e_sc) * 100, gain_mix_pol_pp=(e_mix_pol - e_pol) * 100,
                                err_red_scalar_pp=(rs_ - e_sc) * 100, err_red_pol_pp=(e_red_pol - e_pol) * 100)
            print(f'   mixing gain scalar {(ms - e_sc) * 100:+.2f} pol {(e_mix_pol - e_pol) * 100:+.2f} pp | reduced err scalar {(rs_ - e_sc) * 100:+.2f} pol {(e_red_pol - e_pol) * 100:+.2f} pp', flush=True)
            print(f'{nm:12s} {lam} nm  scalar {e_sc:.4f}  pol-resolved (no s/p conversion) {e_pol:.4f}  diff {100 * (e_pol - e_sc):+.2f} pp  (s share of v0 {d["Ps"].sum():.2f})', flush=True)
        json.dump(res, open(f'{OUT}/results.json', 'w'), indent=1)

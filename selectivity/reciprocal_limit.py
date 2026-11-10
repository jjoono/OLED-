"""Complete reciprocal mixing: closed form and its gap to the memoryless (statistical ray optics) formula.

For B_R(i->r) = (1-p_i) q_r with q = w(1-p)/sum w(1-p) (the only reciprocal rank-one return), the mixed light escapes
with  p_mix = sum q p = p - var_w(p)/(1-p)  (var_w: etendue-weighted variance of p(theta)), and the series has the
closed form  eta = p.v0 + (p.u)(v.v0)/(1 - v.u),  u = R*q, v = 1-p.  Output: data/reciprocal_limit.json"""
import json, numpy as np
src = open('oe_figs.py').read(); exec(src[:src.index('# ------------------------------------------------------------------ Fig. 1')])
w = W1 / W1.sum(); out = []
for nm, (c, BT, BR) in LIB.items():
    pT = BT.sum(1); p = float(pT @ w); var = float(((pT - p) ** 2) @ w)
    q = w * (1 - pT); q /= q.sum(); pmix = float(q @ pT)
    assert abs(pmix - (p - var / (1 - p))) < 1e-12
    for dev in DEVS:
        P = DV[f'{dev}_P'][75].reshape(90, 2).sum(1); R = DV[f'{dev}_R'][75].reshape(90, 2).mean(1)
        u = R * q; v = 1 - pT
        eta_cf = float(pT @ P + (pT @ u) * (v @ P) / (1 - v @ u))
        e = EV[(nm, dev, 550.0)]
        assert abs(eta_cf - e['m1']) < 1e-9, (nm, eta_cf, e['m1'])
        out.append(dict(structure=nm, cls=c, device=dev, p=p, var=var, p_mix=pmix, xi_max=pmix / p,
                        eta_full=e['eta'], eta_mix=eta_cf, eta_m0=e['m0'], eta_memoryless=e['closed'],
                        p_inf=e['pinf'], lam2=e['lam2']))
json.dump(out, open('data/reciprocal_limit.json', 'w'), indent=1)
for r in out:
    if r['device'] == 'Al' and (r['structure'] in ('planar', 'MLA hemi_total', 'MLA free_total', 'scatter g0.8_t1.0', 'scatter g0.8_t4.0', 'A_shape1.0', 'W_slope1.0', 'W_slope0.2', 'C_jitter0.3')):
        print(f"{r['structure']:20s} p {r['p']:.3f} sd {np.sqrt(r['var']):.3f} p_mix {r['p_mix']:.3f} p_inf {r['p_inf']:.3f} | eta full {r['eta_full']:.3f} m0 {r['eta_m0']:.3f} mix {r['eta_mix']:.3f} memless {r['eta_memoryless']:.3f}")
import collections
for dev in DEVS:
    rr = [r for r in out if r['device'] == dev and r['cls'] != 'planar']
    g = np.array([r['eta_mix'] - r['eta_full'] for r in rr]) * 100; o = np.array([r['eta_memoryless'] - r['eta_mix'] for r in rr]) * 100
    print(dev, 'mixing headroom (pp) %.2f..%.2f' % (g.min(), g.max()), ' memoryless over reciprocal-mixed %.2f..%.2f' % (o.min(), o.max()))
    for cl in ('scatter', 'MLA hemi', 'MLA free', 'random lens', 'wavy'):
        s = [r for r in rr if r['cls'] == cl]
        print('   ', cl, 'p_mix/p %.2f-%.2f' % (min(r['xi_max'] for r in s), max(r['xi_max'] for r in s)),
              'headroom %.2f..%.2f' % (min(r['eta_mix'] - r['eta_full'] for r in s) * 100, max(r['eta_mix'] - r['eta_full'] for r in s) * 100),
              'memless-mix %.2f..%.2f' % (min(r['eta_memoryless'] - r['eta_mix'] for r in s) * 100, max(r['eta_memoryless'] - r['eta_mix'] for r in s) * 100))

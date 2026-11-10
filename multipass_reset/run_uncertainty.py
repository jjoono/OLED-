"""Monte Carlo uncertainty of the return matrices of representative random structures.
  ray  : same geometry as random_sweep/, ray seeds 2 and 3 (seed 1 is the sweep itself)
  real : new geometry realizations with the same statistical parameters (geometry rng 12, 13), ray seed 1
Output: random_sweep_unc/<name>__<kind><k>.npz"""
import os, numpy as np
import hftrace_gen as G
OUT = 'random_sweep_unc'; os.makedirs(OUT, exist_ok=True)
TARGETS = ('A_shape0.0', 'A_shape1.0', 'C_jitter0.3', 'B_aspect0.6', 'W_slope0.2', 'W_slope1.0')
src = open('run_random_sweep.py').read()
build = src[src.index("    rng = np.random.default_rng(11)"):src.index("    for name, fn in jobs:")]
build = '\n'.join(l[4:] for l in build.split('\n'))
base = src[:src.index("if __name__ == '__main__':")]


def jobs_for(geo_seed, ray_seed, wave_seed):
    ns = {}
    exec(base, ns)
    g = ns['G']
    lens0, wave0 = g.bsdf_lens, g.bsdf_wave
    ns['G'] = type('G', (), {'bsdf_lens': staticmethod(lambda H, S, **kw: lens0(H, S, seed=ray_seed, **kw)),
                             'bsdf_wave': staticmethod(lambda s, **kw: wave0(s, seed=ray_seed, surf_seed=wave_seed, **kw)), 'A': g.A})
    exec(build.replace('default_rng(11)', f'default_rng({geo_seed})'), ns)
    return dict(ns['jobs'])


if __name__ == '__main__':
    plan = [('ray2', 11, 2, 0), ('ray3', 11, 3, 0), ('real1', 12, 1, 1), ('real2', 13, 1, 2)]
    for tag, gs, rs, ws in plan:
        J = jobs_for(gs, rs, ws)
        for nm in TARGETS:
            f = f'{OUT}/{nm}__{tag}.npz'
            if os.path.exists(f):
                continue
            BT, BR, lost = J[nm](); np.savez(f, BT=BT, BR=BR, lost=lost); print(nm, tag, 'done', flush=True)

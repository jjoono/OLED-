"""Textured exit faces on slabs of higher index (Sec. 4.2, Note S11).
Slab of index n with an isotropic internal source, an index-matched textured exit face (hemispherical MLA or Gaussian
wavy surface, rms slope 0.4 / 1.0; traced by ../multipass_reset/run_highindex.py) and an angle-independent specular
back reflectance R. Compares no mixing, complete reciprocal mixing at the exit face, and complete mixing at the back
reflector (= complete transmissive mixing = fully randomized model for constant R). Writes data/highindex.json."""
import os, json, numpy as np
D = '../multipass_reset/highindex'
E = np.radians(np.arange(91)); w = np.diff(np.sin(E) ** 2); w /= w.sum(); I = np.eye(90)
out = []
for name in ('hemiMLA', 'wavy0.4', 'wavy1.0'):
    for n in (1.51, 2.4, 3.5):
        f = f'{D}/{name}_n{n}.npz'
        if not os.path.exists(f): continue
        z = np.load(f); pT = z['BT'].sum(1); BR = z['BR'].T                    # BR[out, in]
        p = float(w @ pT); q = w * (1 - pT); q /= q.sum(); pm = float(q @ pT); s2 = float(w @ (pT - p) ** 2)
        for R in (0.90, 0.98):
            none = float(pT @ np.linalg.solve(I - R * BR, w))
            exit_mix = p + (1 - p) * R * pm / (1 - R * (1 - pm))
            rand = p / (p + (1 - p) * (1 - R))
            out.append(dict(structure=name, n=n, R=R, p=p, p_mix=pm, var=s2, none=none, exit_mix=exit_mix, rand=rand))
json.dump(out, open('data/highindex.json', 'w'), indent=1)
for r in out:
    print(f"{r['structure']:8s} n {r['n']:4} R {r['R']}: p {r['p']:.3f} p_mix {r['p_mix']:.3f} | none {r['none']:.3f} exit-mix {r['exit_mix']:.3f} back/rand {r['rand']:.3f}")

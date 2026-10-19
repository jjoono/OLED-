"""BSDF of a hexagonal cone array (base radius 1.4 um, height 1.8 um, 0.2 um gap between base circles -> pitch 3.0 um,
flat glass/air in the gaps), index-matched to the substrate, and of a hemisphere array with the same base radius and gap
for comparison.  Geometric optics (feature size >> lambda/n).  Output BSDF(180, 90) in the author's LightTools layout
(rows 0-89 transmitted, rows 90-179 reflected in reverse order), one 1-degree incidence bin per column.
-> cone_bsdf.npz  keys 'cone_n<n>', 'hemi_gap_n<n>'"""
import sys, numpy as np
from multiprocessing import Pool
import raytrace as rt

R_UM, H_UM, GAP_UM = 1.4, 1.8, 0.2
GEOM = {'cone': ('cone', H_UM / R_UM, (2 * R_UM + GAP_UM) / R_UM), 'hemi_gap': ('cap', 1.0, (2 * R_UM + GAP_UM) / R_UM)}
N = 40000
TH = np.arange(90) + 0.5


def col(a):
    key, n, i = a; shape, AR, pitch = GEOM[key]
    esc, rc = rt._trace_chunk(AR, n, N, 80, np.random.default_rng([7, i, int(n * 100)]), float(TH[i]), shape, pitch)
    good = np.isfinite(rc)
    h, _ = np.histogram(np.degrees(np.arccos(np.clip(rc[good], 0, 1))), bins=np.linspace(0, 90, 91))
    return key, n, i, esc.mean(), h / N


if __name__ == '__main__':
    ns = [float(x) for x in sys.argv[1:]] or [1.8]
    jobs = [(k, n, i) for k in GEOM for n in ns for i in range(90)]
    out = {}
    with Pool(4) as p:
        for key, n, i, bt, br in p.imap_unordered(col, jobs):
            B = out.setdefault(f'{key}_n{n:g}', np.zeros((180, 90)))
            B[0, i] = bt                     # transmitted total kept in row 0 (air-angle detail not needed by the series)
            B[90:180, i] = br[::-1]
    try:
        old = dict(np.load('cone_bsdf.npz'))
    except FileNotFoundError:
        old = {}
    old.update(out); np.savez('cone_bsdf.npz', **old)
    G = np.cos(np.radians(TH)) * np.sin(np.radians(TH)); G /= G.sum()
    for k, B in sorted(out.items()):
        BT = B[:90].sum(0); BR = B[90:][::-1]
        print('%-14s p_Lamb %.3f | B_T at 0.5/20/40/50/60/70/80: %s | closure %.4f-%.4f' % (
            k, BT @ G, ' '.join('%.2f' % BT[j] for j in (0, 20, 40, 50, 60, 70, 80)), (BT + BR.sum(0)).min(), (BT + BR.sum(0)).max()))

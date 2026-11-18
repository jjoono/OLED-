"""Prospective test of the applicability conditions of the reduced model (main text Sec. 4.3, Fig. 3c).
Twelve structures that were not used before: hexagonal cone arrays (side angle 30-70 deg), square pyramids (30, 45,
54.7 deg), inverted pyramids (54.7 deg) and one-dimensional prism gratings (30, 45, 60 deg).
Step 1 (prediction): from M and the first-pass distribution (Al, 550 nm) only, each structure is classified as
 'valid' (irreducible, lambda_2 - lambda_eff <= 0.35) or 'outside' (reducible: lambda_2 > 0.999, or gap > 0.35).
 The predictions are written to prospective/predictions.json before any multi-pass series is evaluated.
Step 2 (test): reduced-model error against the full series for three devices and three wavelengths.
Usage: python3 prospective.py trace | predict | test"""
import os, sys, json, time, numpy as np
from numba import njit
import hftrace as T, hftrace_gen as G
OUT = 'prospective'; os.makedirs(OUT, exist_ok=True)
r = np.linspace(0, 1, T.NT) if hasattr(T, 'NT') else np.linspace(0, 1, 2001)


@njit(cache=True)
def _surf_pyr(x, y, P0, P1, P2, P3, P4, P5):
    # P0[0] = tan(side angle); P0[1] = +1 pyramid, -1 inverted pyramid, 0 prism (ridges along y); period 2
    H = P0[0]; mode = P0[1]
    u = (x % 2.0) - 1.0; v = (y % 2.0) - 1.0
    if mode == 0.0:
        return H * (1.0 - abs(u)), H, -np.sign(u), 0.0
    if abs(u) >= abs(v):
        m = abs(u); gx = -np.sign(u); gy = 0.0
    else:
        m = abs(v); gx = 0.0; gy = -np.sign(v)
    if mode > 0:
        return H * (1.0 - m), H, gx, gy
    return H * m, H, -gx, -gy


_trace_pyr = G._make_trace(_surf_pyr); _batch_pyr = G._make_batch(_trace_pyr)


def hexcone(alpha):
    h = np.tan(np.radians(alpha)) * (1 - r); htab, dtab = T.tables(r, h)
    return ('hex', (htab, dtab), float(htab.max()) + 1e-6, float(np.abs(dtab).max()))


def pyr(alpha, mode):
    H = np.tan(np.radians(alpha)); P = (np.array([H, float(mode)]), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1))
    return ('pyr', P, H + 1e-6, H)


STRUCT = {**{f'cone{a}': hexcone(a) for a in (30, 40, 50, 55, 70)},
          **{f'pyramid{a:g}': pyr(a, 1) for a in (30, 45, 54.7)}, 'inv_pyramid54.7': pyr(54.7, -1),
          **{f'prism{a}': pyr(a, 0) for a in (30, 45, 60)}}


def trace(name, N=20000):
    kind, P, hmax, lip = STRUCT[name]
    if kind == 'hex':
        BT = np.zeros((90, 90)); BR = np.zeros((90, 90)); lost = np.zeros(90); e = np.linspace(0, 90, 91)
        for i in range(90):
            oc, to = T._batch(np.radians(i + .5), N, 1.51, P[0], P[1], hmax, lip, 80, 7919 + i * 100003)
            d = np.degrees(to); BT[i] = np.histogram(d[oc == 0], e)[0] / N; BR[i] = np.histogram(d[oc == 1], e)[0] / N; lost[i] = np.mean(oc == 2)
    else:
        BT, BR, lost = G._run(_batch_pyr, P, (2.0, 0.0, 0.0, 2.0), hmax, lip, 1.51, N, 80, 1)
    np.savez(f'{OUT}/{name}.npz', BT=BT, BR=BR, lost=lost)
    return BT, BR, lost


if __name__ == '__main__':
    step = sys.argv[1]
    if step == 'trace':
        for nm in (sys.argv[2:] or STRUCT):
            if os.path.exists(f'{OUT}/{nm}.npz'): continue
            t = time.time(); BT, BR, lost = trace(nm)
            W = np.diff(np.sin(np.radians(np.linspace(0, 90, 91))) ** 2)
            print(f'{nm}: {time.time()-t:.0f}s p {BT.sum(1) @ W / W.sum():.4f} lost max {lost.max():.4f}', flush=True)

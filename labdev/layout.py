"""Lab devices in the actual substrate layout (25.4 x 25.4 mm, five Al fingers, ITO on the upper half, PDL-defined
pixel r = 1 mm at the centre of the upper end-cap of the middle finger; black sides). Coordinates in mm, origin at the
pixel centre, +y towards the ITO edge of the substrate. Back-face reflectance by region: pixel (full stack), Al-covered
(fingers, optionally the upper rectangle) with/without ITO, uncovered with/without ITO. Lens patch: square, centred on
the pixel. Output counts light leaving the front face; spectral and Purcell weighting as in simulate.py."""
import os, json, numpy as np
import simulate as SM
TAN = SM.TAN
X0, X1, Y0, Y1 = -12.7, 12.7, -14.4882, 10.9118          # substrate edges relative to the pixel centre
FX = np.array([-7.9796, -3.9898, 0.0, 3.9898, 7.9796]); FR = 1.3919; FL = 9.9263
RECT = (-10.9999, 10.9999, 4.59, 9.16); Y_ITO = -1.7882


def al_mask(x, y, rect_al):
    m = np.zeros(x.shape, bool)
    for xc in FX:
        yy = np.clip(y, -FL, 0.0)
        m |= np.hypot(x - xc, y - yy) <= FR
    if rect_al:
        m |= (x >= RECT[0]) & (x <= RECT[1]) & (y >= RECT[2]) & (y <= RECT[3])
    return m


def mc(S, PL, R, P, t_mm, half_patch, rect_al, N=100000, seed=0, maxit=5000):
    rng = np.random.default_rng(seed)
    th = np.minimum(np.searchsorted(np.cumsum(P) / P.sum(), rng.random(N)), 89)
    r0 = np.sqrt(rng.random(N)); f0 = 2 * np.pi * rng.random(N)
    x = r0 * np.cos(f0); y = r0 * np.sin(f0); phi = 2 * np.pi * rng.random(N); g = np.ones(N)
    x += t_mm * TAN[th] * np.cos(phi); y += t_mm * TAN[th] * np.sin(phi)
    eS = np.zeros(90); eP = np.zeros(90); side = 0.0
    for it in range(maxit):
        out = (x < X0) | (x > X1) | (y < Y0) | (y > Y1); side += g[out].sum()
        k = ~out; th, x, y, phi, g = th[k], x[k], y[k], phi[k], g[k]
        if g.size == 0: break
        ins = (np.abs(x) <= half_patch) & (np.abs(y) <= half_patch)
        pp = np.where(ins, S['p'][th], PL['p'][th])
        np.add.at(eS, th[ins], (g * pp)[ins]); np.add.at(eP, th[~ins], (g * pp)[~ins]); g = g * (1 - pp)
        u = rng.random(g.size); new = th.copy()
        if ins.any() and S.get('cond') is not None:
            new[ins] = np.minimum((S['cdf'][th[ins]] < u[ins, None]).sum(1), 89)
            c = S['cond'][th[ins], new[ins]]; kk = np.minimum((c < rng.random(ins.sum())[:, None]).sum(1), S['NA'] - 1)
            phi[ins] = phi[ins] + (kk + rng.random(ins.sum())) * 2 * np.pi / S['NA']
        th = new
        x += t_mm * TAN[th] * np.cos(phi); y += t_mm * TAN[th] * np.sin(phi)
        out = (x < X0) | (x > X1) | (y < Y0) | (y > Y1); side += g[out].sum(); g[out] = 0.0
        pix = np.hypot(x, y) <= 1.0; al = al_mask(x, y, rect_al); ito = y >= Y_ITO
        Rv = np.where(pix, R['pix'][th], np.where(al, np.where(ito, R['ITO_Al'][th], R['noITO_Al'][th]),
                                                  np.where(ito, R['ITO_noAl'][th], R['noITO_noAl'][th])))
        g = g * Rv
        x += t_mm * TAN[th] * np.cos(phi); y += t_mm * TAN[th] * np.sin(phi)
        k = g > 1e-7; th, x, y, phi, g = th[k], x[k], y[k], phi[k], g[k]
    return eS / N, eP / N, side / N


if __name__ == '__main__':
    dv = np.load('device.npz'); st = np.load('layout_stacks.npz'); L = dv['lam']; s = dv['s']; F = dv['F']; fsub = dv['f_sub']
    Q = SM.Q; wq = Q * s * F * fsub / (Q * np.sum(s * F) + 1 - Q)
    PL = SM.planar_kernel(); KER = {'bare': PL, 'hemisphere': SM.lens_kernel('hemisphere'), 'freeform': SM.lens_kernel('freeform')}
    E = np.arange(91); dOm = 2 * np.pi * (np.cos(np.radians(E[:-1])) - np.cos(np.radians(E[1:])))
    BTn = lambda K: K['BT'] / np.maximum(K['p'][:, None], 1e-15)
    res = {}; cases = [('bare', 0.7, 0.0)] + [(nm, 1.4, hp) for nm in ('hemisphere', 'freeform') for hp in (1.5, 5.0)]
    for rect_al in (True, False):
        for nm, t, hp in cases:
            air = np.zeros((len(L), 90)); side = np.zeros(len(L))
            for k in range(len(L)):
                R = {'pix': dv['R_pix'][k], **{n: st[n][k] for n in st.files}}
                eS, eP, sd = mc(KER[nm], PL, R, dv['P'][k], t, hp, rect_al, seed=k)
                air[k] = eS @ BTn(KER[nm]) + eP @ BTn(PL); side[k] = sd
            A = wq[:, None] * air; I0 = A[:, :5].sum(1) / dOm[:5].sum(); eqe = float(A.sum()); fl = float(np.pi * I0.sum())
            key = f"{nm}|patch{2 * hp:g}mm|rect{'Al' if rect_al else 'noAl'}"
            res[key] = dict(EQE=eqe, front_EQE=fl, PE=eqe / fl, side=float((wq * side).sum() / wq.sum()),
                            bands={b: float(A[:, lo:hi].sum()) for b, (lo, hi) in {'0-20': (0, 20), '20-40': (20, 40), '40-60': (40, 60), '60-80': (60, 80), '80-90': (80, 90)}.items()},
                            profile=(A.sum(0) / dOm).tolist())
            print(f"{key:36s} EQE {eqe*100:.1f}%  front {fl*100:.1f}%  PE {eqe/fl:.2f}  side {res[key]['side']:.3f}", flush=True)
    json.dump(res, open('layout_results.json', 'w'))

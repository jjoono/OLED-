"""Scan-plane-resolved goniometer model for the green devices.
Light transport: layout.mc extended to absolute azimuth. Inside the lens patch the lattice-resolved kernels
(revtrace.kernel_abs) are used with the lattice rotated by beta from the substrate x axis (x across the Al fingers,
y along them, +y towards the ITO edge). Output: radiance L(theta, phi) in the substrate frame and the exit-position map.
Goniometer: square PD of side s at D = 59 mm from the rotation axis point C = (0, 0, z_c) (z from the emitter plane),
scan plane at azimuth phi_s. At goniometer angle theta_g a photon leaving P = (x, y, t) is collected when its
in-plane angle lies within theta_g - b/D +- w and its out-of-plane angle within -y_perp/D +- w, w = s/2D,
b = x_s cos(theta_g) - (t - z_c) sin(theta_g). Measured PE = 2 int (I/I0) sin(theta) d(theta) from the 10-degree points
(natural spline, I(90) = 0), as in the lab spreadsheet.
Usage: python3 gonio_az.py ETL HTL out.json"""
import sys, json, numpy as np
from scipy.interpolate import CubicSpline
import simulate as SM, layout as LY, green_joint_fit as GJ
sys.path.insert(0, '../multipass_reset'); from revtrace import reduce_az
TAN = SM.TAN; NA = 72; NB = 6; DPS = 2 * np.pi / NA; XB = np.linspace(-10, 10, 101); XC = 0.5 * (XB[1:] + XB[:-1])
D = 59.0


def load_abs(nm):
    z = np.load(f'kernelabs_{nm}.npz'); KTa = z['KTa'].astype(float); Ka = z['Ka'].astype(float)
    p = KTa.sum((2, 3)); r = Ka.sum((2, 3))
    c = np.cumsum(Ka.reshape(90 * NB, 90 * NA), 1) / np.maximum(r.reshape(-1, 1), 1e-15); c[:, -1] = 1.0
    flat = (c + np.arange(90 * NB)[:, None]).ravel()
    return dict(KTa=KTa, Kf=np.fft.fft(KTa, axis=3), p=p, r=r, flat=flat)


def mc_abs(K, PL, R, P, t, hp, beta, N, seed, maxit=5000):
    rng = np.random.default_rng(seed)
    st = lambda: rng.permutation((np.arange(N) + rng.random(N)) / N)          # stratified uniforms
    th = np.minimum(np.searchsorted(np.cumsum(P) / P.sum(), st()), 89)
    r0 = np.sqrt(st()); f0 = 2 * np.pi * st()
    x = r0 * np.cos(f0); y = r0 * np.sin(f0); phi = 2 * np.pi * st(); g = np.ones(N)
    x += t * TAN[th] * np.cos(phi); y += t * TAN[th] * np.sin(phi)
    Wl = np.zeros((90, NB, 2, NA)); Wp = np.zeros((90, NA)); H = np.zeros((100, 100))
    for it in range(maxit):
        out = (x < LY.X0) | (x > LY.X1) | (y < LY.Y0) | (y > LY.Y1); k = ~out
        th, x, y, phi, g = th[k], x[k], y[k], phi[k], g[k]
        if g.size == 0: break
        ins = (np.abs(x) <= hp) & (np.abs(y) <= hp) if K is not None else np.zeros(g.size, bool)
        pb = ((phi % (2 * np.pi)) / DPS).astype(int) % NA
        if ins.any():
            a, s = reduce_az(phi[ins] - beta); ab = np.minimum((a / 5.0).astype(int), NB - 1); sb = (s < 0).astype(int)
            np.add.at(Wl, (th[ins], ab, sb, pb[ins]), g[ins])
            pl = K['p'][th[ins], ab]; rl = K['r'][th[ins], ab]
        np.add.at(Wp, (th[~ins], pb[~ins]), g[~ins])
        pe = PL['p'][th].copy(); rr = 1 - pe
        if ins.any(): pe[ins] = pl; rr[ins] = rl
        H += np.histogram2d(x, y, (XB, XB), weights=g * pe)[0]
        g = g * rr; new = th.copy()
        if ins.any():
            row = th[ins] * NB + ab; idx = np.searchsorted(K['flat'], rng.random(row.size) + row) - row * 90 * NA
            idx = np.clip(idx, 0, 90 * NA - 1); new[ins] = idx // NA
            phi[ins] = phi[ins] + s * ((idx % NA) + rng.random(row.size)) * DPS
        th = new; x += t * TAN[th] * np.cos(phi); y += t * TAN[th] * np.sin(phi)
        out = (x < LY.X0) | (x > LY.X1) | (y < LY.Y0) | (y > LY.Y1); g[out] = 0.0
        pix = np.hypot(x, y) <= 1.0; al = LY.al_mask(x, y, True); ito = y >= LY.Y_ITO
        g = g * np.where(pix, R['pix'][th], np.where(al, np.where(ito, R['ITO_Al'][th], R['noITO_Al'][th]),
                                                      np.where(ito, R['ITO_noAl'][th], R['noITO_noAl'][th])))
        x += t * TAN[th] * np.cos(phi); y += t * TAN[th] * np.sin(phi)
        k = g > 1e-7; th, x, y, phi, g = th[k], x[k], y[k], phi[k], g[k]
    # radiant power per (theta_out 1 deg, phi grid point k * DPS)
    BTp = PL['BT']; Lp = np.einsum('if,io->of', Wp, BTp); L = 0.5 * (Lp + np.roll(Lp, 1, axis=1))
    if K is not None:
        Wf = np.fft.fft(Wl, axis=3); sh = np.exp(-2j * np.pi * np.arange(NA) / NA)
        Lf = np.einsum('iaf,iaof->of', Wf[:, :, 0], K['Kf']) * sh + np.einsum('iaf,iaof->of', Wf[:, :, 1], np.conj(K['Kf']))
        L = L + np.real(np.fft.ifft(Lf, axis=1))
    return L / N, H / N


def gonio(L, H, t, zc, phis, s=5.0, thg=np.arange(-85, 86, 1.0), nw=5):
    th_c = np.arange(90) + 0.5; dO = (np.cos(np.radians(th_c - 0.5)) - np.cos(np.radians(th_c + 0.5))) * DPS
    Lr = L / dO[:, None]; Lr = np.c_[Lr, Lr[:, :1]]                       # periodic in phi
    m = H > H.max() * 1e-4; X, Y = np.meshgrid(XC, XC, indexing='ij'); xs0, ys0, wt = X[m], Y[m], H[m]
    ps = np.radians(phis); xs = xs0 * np.cos(ps) + ys0 * np.sin(ps); yp = -xs0 * np.sin(ps) + ys0 * np.cos(ps)
    w = s / 2 / D; al = np.linspace(-w, w, nw); I = np.zeros(thg.size)
    for i, tg in enumerate(np.radians(thg)):
        b = xs * np.cos(tg) - (t - zc) * np.sin(tg)
        z = (tg - b / D)[:, None, None] + al[None, :, None]; e = (-yp / D)[:, None, None] + al[None, None, :]
        dx = np.cos(e) * np.sin(z) * np.cos(ps) - np.sin(e) * np.sin(ps); dy = np.cos(e) * np.sin(z) * np.sin(ps) + np.sin(e) * np.cos(ps)
        dz = np.cos(e) * np.cos(z); th = np.degrees(np.arccos(np.clip(dz, -1, 1))); ph = np.arctan2(dy, dx) % (2 * np.pi) / DPS
        ti = np.clip(th - 0.5, 0, 88.999); t0 = ti.astype(int); ft = ti - t0; p0 = ph.astype(int) % NA; fp = ph - np.floor(ph)
        v = ((1 - ft) * ((1 - fp) * Lr[t0, p0] + fp * Lr[t0, p0 + 1]) + ft * ((1 - fp) * Lr[t0 + 1, p0] + fp * Lr[t0 + 1, p0 + 1]))
        I[i] = (wt[:, None, None] * v).sum() / nw ** 2
    return thg, I


def summarize(thg, I):
    I0 = I[thg == 0][0]; pos = np.array([I[thg == a][0] for a in range(0, 81, 10)]) / I0
    neg = np.array([I[thg == -a][0] for a in range(0, 81, 10)]) / I0; sym = 0.5 * (pos + neg)
    a = np.radians(np.r_[np.arange(0, 81, 10), 90]); cs = CubicSpline(a, np.r_[sym, 0.0], bc_type='natural')
    tt = np.linspace(0, np.pi / 2, 2001); PE = 2 * np.trapezoid(cs(tt) * np.sin(tt), tt) if hasattr(np, 'trapezoid') else None
    return dict(sym=sym.tolist(), pos=pos.tolist(), neg=neg.tolist(), PE=float(PE))


if __name__ == '__main__':
    etl, htl = float(sys.argv[1]), float(sys.argv[2]); N = 300000; NCH = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    ST, wq = GJ.stacks(etl, htl); PL = SM.planar_kernel()
    KS = {'hemisphere': load_abs('hemisphere'), 'freeform': load_abs('green_freeform')}
    cases = [('bare', None, 0.7, 0.0, 0.0)] + [(nm, KS[nm], 1.4, 1.5, b) for nm in KS for b in (0.0, 30.0)]
    res = dict(etl=etl, htl=htl, N=N * NCH)
    for nm, K, t, hp, beta in cases:
        L = 0; H = 0
        for c in range(NCH):
            for j, (P, R) in enumerate(ST):
                Lj, Hj = mc_abs(K, PL, R, P, t, hp, np.radians(beta), N, seed=1000 * c + 10 + j); L = L + wq[j] * Lj / NCH; H = H + wq[j] * Hj / NCH
        th_c = np.arange(90) + 0.5; dOm = np.cos(np.radians(th_c - 0.5)) - np.cos(np.radians(th_c + 0.5))
        az = L.sum(1) / (2 * np.pi * dOm); ref = GJ.rel(az)[0]
        key = f'{nm}|beta{beta:g}'; res[key] = dict(EQE=float(L.sum()), azavg_rel=ref.tolist(), L=L.tolist(), H=H.tolist(), scans={})
        print(f'{key}: EQE {L.sum()*100:.1f}%  azimuth-averaged I/I0 ' + ' '.join(f'{v:.3f}' for v in ref[1:7]), flush=True)
        for zc_name, zc in (('device', 0.0), ('middle', t / 2), ('surface', t)):
            for phis in (0.0, 45.0, 90.0):
                for s in (5.0, 10.0):
                    thg, I = gonio(L, H, t, zc, phis, s); sm = summarize(thg, I)
                    res[key]['scans'][f'{zc_name}|phis{phis:g}|s{s:g}'] = sm
            sm = res[key]['scans'][f'{zc_name}|phis0|s5']; sm9 = res[key]['scans'][f'{zc_name}|phis90|s5']
            print(f'   axis {zc_name:7s} scan x : ' + ' '.join(f'{v:.3f}' for v in sm['sym'][1:7]) + f"  PE {sm['PE']:.3f}  +/-asym@40 {sm['pos'][4]/sm['neg'][4]:.3f}", flush=True)
            print(f'   axis {zc_name:7s} scan y : ' + ' '.join(f'{v:.3f}' for v in sm9['sym'][1:7]) + f"  PE {sm9['PE']:.3f}  +/-asym@40 {sm9['pos'][4]/sm9['neg'][4]:.3f}", flush=True)
    json.dump(res, open(sys.argv[3], 'w'))

"""Freeform lenslet profile digitized from the LightTools profile plot (7 control points, curve pixels).
Axes were calibrated on the five interior control-point markers (residual < 0.4 px, ~0.003 lens radii).
Returns h(rho) on rho in [0, 1] (lens-radius units); the vertical sidewall at rho = 1 (h = 0 ... 0.195)
is represented by a linear drop over the last DRHO of rho."""
import numpy as np
from scipy.signal import savgol_filter
CP_Z = np.array([2.2779, 2.4307, 2.4922, 1.8176, 0.98070, 0.19474, 0.0])
CP_Y = np.array([0.0, 0.30919, 0.41117, 0.73240, 0.90834, 1.0, 1.0])


def profile(drho=0.002, n=4001):
    """Least-squares cubic B-spline (r(t), z(t)) through the digitized pixels, parameterized by the projection of each
    pixel on the 7-point chord-length spline; then h(rho) with rho clipped to <= 1 and made single-valued."""
    from scipy.interpolate import CubicSpline, make_lsq_spline
    Zc, Yc = np.load(__file__.replace('profile.py', 'freeform_curve_digitized.npy')).T
    s = np.r_[0, np.cumsum(np.hypot(np.diff(CP_Z), np.diff(CP_Y)))]
    tt = np.linspace(0, s[-1], 4000); cr = CubicSpline(s, CP_Y)(tt); cz = CubicSpline(s, CP_Z)(tt)
    d2 = (Yc[:, None] - cr[None]) ** 2 + (Zc[:, None] - cz[None]) ** 2
    keep = d2.min(1) < 0.03 ** 2; tp = tt[d2.argmin(1)][keep]; Y, Z = Yc[keep], Zc[keep]
    # anchor the end points with heavy weight
    tp = np.r_[tp, [0.0] * 20, [s[-1]] * 20, [s[-2]] * 20]; Y = np.r_[Y, [0.0] * 20, [1.0] * 20, [1.0] * 20]; Z = np.r_[Z, [CP_Z[0]] * 20, [0.0] * 20, [CP_Z[-2]] * 20]
    o = np.argsort(tp, kind='stable'); tp, Y, Z = tp[o], Y[o], Z[o]
    k = 3; inner = np.linspace(0, s[-1], 18)[1:-1]; knots = np.r_[[0.0] * (k + 1), inner, [s[-1]] * (k + 1)]
    sr = make_lsq_spline(tp, Y, knots, k); sz = make_lsq_spline(tp, Z, knots, k)
    t = np.linspace(0, s[-1], 20000); r = np.clip(sr(t), 0, 1); z = sz(t)
    r = np.maximum.accumulate(r)                           # single-valued h(rho)
    rho = np.linspace(0, 1, n); rmax = 1.0 - drho
    ok = r < rmax; ur, ui = np.unique(r[ok], return_index=True); h = np.interp(rho, ur, z[ok][ui])
    h_wall = float(np.interp(rmax, ur, z[ok][ui])); wall = rho > rmax
    h[wall] = h_wall * (1 - (rho[wall] - rmax) / drho)
    return rho, np.clip(h, 0, None)


if __name__ == '__main__':
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from scipy.interpolate import CubicSpline
    rho, h = profile()
    s = np.r_[0, np.cumsum(np.hypot(np.diff(CP_Z), np.diff(CP_Y)))]; t = np.linspace(0, s[-1], 400)
    fig, ax = plt.subplots(1, 2, figsize=(8, 3.4))
    Zc, Yc = np.load('freeform_curve_digitized.npy').T
    for a in ax:
        a.plot(Yc, Zc, '.', color='#bbb', ms=2, label='digitized pixels')
        a.plot(CubicSpline(s, CP_Y)(t), CubicSpline(s, CP_Z)(t), '--', color='C1', lw=1, label='7-point natural spline')
        a.plot(rho, h, 'C0-', lw=1.2, label='h(ρ) used'); a.plot(CP_Y, CP_Z, 'ro', ms=3, label='control points')
    ax[0].set_aspect('equal'); ax[0].set_xlabel('ρ = r / r_lens'); ax[0].set_ylabel('z / r_lens'); ax[0].legend(fontsize=6)
    ax[1].set_xlim(0.85, 1.01); ax[1].set_ylim(-0.05, 1.3); ax[1].set_title('sidewall detail', fontsize=8)
    fig.tight_layout(); fig.savefig('freeform_profile_check.png', dpi=160)
    print('h(0) %.3f, max %.3f at rho %.3f, h(0.99) %.3f' % (h[0], h.max(), rho[h.argmax()], np.interp(0.99, rho, h)))

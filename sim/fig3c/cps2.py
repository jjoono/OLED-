"""Dipole CPS power budget for an arbitrary stack of uniaxial layers.

Generalisation of sim/fig3b/cps.py.  Same converged quadrature (u = sin(th)
below the light line, u = sqrt(1+v^2) above it, composite Simpson, exact
channel boundaries), but every layer may be uniaxial with the optic axis along
the surface normal.

Convention, taken from the author's TMF_birefringence_whole.m so that the two
agree term by term:

    p (TM):  k_z,j = k0 n_o,j sqrt(1 - (k_x/k0/n_e,j)^2),   k_x = n_e,EML u k0
    s (TE):  k_z,j = sqrt(n_o,j^2 k0^2 - k_x^2),            k_x = n_o,EML u k0
    Im(k_z) >= 0 on both branches (decaying evanescent wave)
    r_p = (n_o,j+1^2 k_z,j - n_o,j^2 k_z,j+1)/(...+...)     t_p = 2 n_o,j n_o,j+1 k_z,j/(...)
    r_s = (k_z,j - k_z,j+1)/(k_z,j + k_z,j+1)               t_s = 1 + r_s

i.e. the isotropic formulae with n -> n_o, and n_e entering only through k_z,p.
"""
import numpy as np

EDGE = 1e-5          # how closely the branch point u = 1 is approached


def _fix(kz):
    return np.where(kz.imag < 0, -kz, kz)


def _kz(no, ne, kx, k0, pol):
    no = np.asarray(no, dtype=complex)
    ne = np.asarray(ne, dtype=complex)
    if pol == 'p':
        return _fix(k0 * no * np.sqrt(1.0 - (kx / (k0 * ne)) ** 2 + 0j))
    return _fix(np.sqrt((no * k0) ** 2 - kx ** 2 + 0j))


def stack_rt(no, ne, ds, kx, k0, pol):
    """r, t of a stack referenced to a source plane inside layer 0.
    ds[0] is the distance from that plane to the first interface; ds[-1] unused."""
    L = len(no)
    kz = [_kz(no[j], ne[j], kx, k0, pol) for j in range(L)]
    n2 = [np.asarray(no[j], dtype=complex) ** 2 for j in range(L)]

    def fres(j):
        a, b = j, j + 1
        if pol == 'p':
            d = n2[b] * kz[a] + n2[a] * kz[b]
            return ((n2[b] * kz[a] - n2[a] * kz[b]) / d,
                    2 * no[a] * no[b] * kz[a] / d)
        d = kz[a] + kz[b]
        return (kz[a] - kz[b]) / d, 2 * kz[a] / d

    r, t = fres(L - 2)
    for j in range(L - 3, -1, -1):
        ph = np.exp(1j * kz[j + 1] * ds[j + 1])
        ph2 = ph * ph
        rj, tj = fres(j)
        den = 1.0 + rj * r * ph2
        t = tj * t * ph / den
        r = (rj + r * ph2) / den
    ph0 = np.exp(1j * kz[0] * ds[0])
    return r * ph0 * ph0, t * ph0


class Stack(object):
    """air | above (reversed: nearest the EML first) | EML | below | substrate

    `above` and `below` are lists of (n_o, n_e, thickness) ordered outwards from
    the EML.  The dipole sits z0 above the bottom face of the EML."""

    def __init__(self, lam, eml, d_eml, z0, above, below, n_sub, n_top=1.0):
        self.lam = lam
        self.no_e, self.ne_e = eml
        self.d_eml, self.z0 = d_eml, z0
        self.n_sub = complex(n_sub)
        self.bot_no = [self.no_e] + [l[0] for l in below] + [n_sub]
        self.bot_ne = [self.ne_e] + [l[1] for l in below] + [n_sub]
        self.bot_d = [d_eml - z0] + [l[2] for l in below] + [0.0]
        self.top_no = [self.no_e] + [l[0] for l in above] + [n_top]
        self.top_ne = [self.ne_e] + [l[1] for l in above] + [n_top]
        self.top_d = [z0] + [l[2] for l in above] + [0.0]
        # whole OLED seen from inside the substrate (for the incoherent substrate)
        self.all_no = [n_sub] + self.bot_no[-2:0:-1] + [self.no_e] + \
                      [l[0] for l in above] + [n_top]
        self.all_ne = [n_sub] + self.bot_ne[-2:0:-1] + [self.ne_e] + \
                      [l[1] for l in above] + [n_top]
        self.all_d = [0.0] + [l[2] for l in below][::-1] + [d_eml] + \
                     [l[2] for l in above] + [0.0]


def _prop(u, S):
    k0 = 2 * np.pi / S.lam
    out = {}
    for pol, nref in (('p', S.ne_e), ('s', S.no_e)):
        kx = nref * u * k0
        rb, tb = stack_rt(S.bot_no, S.bot_ne, S.bot_d, kx, k0, pol)
        rt, _ = stack_rt(S.top_no, S.top_ne, S.top_d, kx, k0, pol)
        out[pol] = (rb, tb, rt)
    return out


def _outcoupled(u, S):
    """T_sub-air / (1 - R_oled R_sub-air), per polarisation."""
    k0 = 2 * np.pi / S.lam
    res = []
    for pol, nref in (('p', S.ne_e), ('s', S.no_e)):
        kx = nref * u * k0
        ro, _ = stack_rt(S.all_no, S.all_ne, S.all_d, kx, k0, pol)
        Ro = np.abs(ro) ** 2
        cS = np.sqrt(1.0 - (kx / (k0 * S.n_sub)) ** 2 + 0j)
        cA = np.sqrt(1.0 - (kx / k0) ** 2 + 0j)
        if pol == 'p':
            rf = (1.0 * cS - S.n_sub * cA) / (1.0 * cS + S.n_sub * cA)
        else:
            rf = (S.n_sub * cS - 1.0 * cA) / (S.n_sub * cS + 1.0 * cA)
        Rf = np.abs(rf) ** 2
        res.append((1.0 - Rf) / (1.0 - Ro * Rf))
    return res[0], res[1]


def _weights(S, h=2.0 / 3.0):
    """The two dipole-orientation prefactors of the source script, with the
    1/lambda^4 dropped (everything is used as a ratio)."""
    no, ne = S.no_e, S.ne_e
    return (1.0 - h) * ne, h * no * (3.0 + (ne / no) ** 2) / 4.0


def integrands_theta(th, S, need_air=True):
    """(W, Sb, A) du/dth on the propagating branch u = sin(th)."""
    no, ne = S.no_e, S.ne_e
    u, c = np.sin(th), np.cos(th)
    P = _prop(u, S)
    rb_p, tb_p, rt_p = P['p']
    rb_s, tb_s, rt_s = P['s']
    Xp = (1 + rb_p) * (1 + rt_p) / (1 - rb_p * rt_p)
    Yp = (1 - rb_p) * (1 - rt_p) / (1 - rb_p * rt_p)
    Xs = (1 + rb_s) * (1 + rt_s) / (1 - rb_s * rt_s)
    Apv = (1 + rt_p) * tb_p / (1 - rb_p * rt_p)
    Aph = (1 - rt_p) * tb_p / (1 - rb_p * rt_p)
    Ash = (1 + rt_s) * tb_s / (1 - rb_s * rt_s)

    cv, ch = _weights(S)
    # kernels already multiplied by cos(th), the Jacobian
    Kpv = 0.75 * (ne / no) * u ** 2 * np.real(Xp)
    Kph = 3.0 / (6 * (no / ne) ** 2 + 2) * c ** 2 * np.real(Yp)
    Ksh = 3.0 / (2 * (ne / no) ** 2 + 6) * np.real(Xs)
    W = u * (cv * Kpv + ch * (Kph + Ksh))

    ns = S.n_sub.real
    inside = u < ns / ne
    safe = np.where(c > 0, c, 1.0)
    cSp = np.sqrt(np.clip(1.0 - (ne * u / ns) ** 2, 0.0, None))      # cos th_sub (p)
    cSs = np.sqrt(np.clip((ns / no) ** 2 - u ** 2, 0.0, None))       # (n_sub/n_o) cos th_sub (s)
    Kpv2 = 0.375 * ne * ns / no ** 2 * cSp * u ** 2 * np.abs(Apv) ** 2 / safe
    Kph2 = 3.0 * (ns / no) * cSp * np.abs(Aph) ** 2 / (12 * (no / ne) ** 2 + 4) * c
    Ksh2 = 3.0 * cSs * np.abs(Ash) ** 2 / ((4 * (ne / no) ** 2 + 12) * safe)
    if abs(ns / ne - 1.0) < 1e-12 and abs(ns / no - 1.0) < 1e-12:
        Kpv2 = 0.375 * ne * ns / no ** 2 * u ** 2 * np.abs(Apv) ** 2
        Kph2 = 3.0 * (ns / no) * np.abs(Aph) ** 2 / (12 * (no / ne) ** 2 + 4) * c ** 2
        Ksh2 = 3.0 * np.abs(Ash) ** 2 / (4 * (ne / no) ** 2 + 12)
    Sb = u * (cv * Kpv2 + ch * (Kph2 + Ksh2)) * inside

    if need_air:
        with np.errstate(divide='ignore', invalid='ignore'):
            Fp, Fs = _outcoupled(u, S)
        Fp = np.where((u < 1.0 / ne) & np.isfinite(Fp), Fp, 0.0)
        Fs = np.where((u < 1.0 / no) & np.isfinite(Fs), Fs, 0.0)
        A = u * (cv * Kpv2 * Fp + ch * (Kph2 * Fp + Ksh2 * Fs)) * inside
    else:
        A = np.zeros_like(W)
    return tuple(np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
                 for x in (W, Sb, A))


def integrand_v(v, S):
    """W du/dv on the evanescent branch u = sqrt(1+v^2)."""
    no, ne = S.no_e, S.ne_e
    u = np.sqrt(1.0 + v ** 2)
    P = _prop(u, S)
    rb_p, _, rt_p = P['p']
    rb_s, _, rt_s = P['s']
    Xp = (1 + rb_p) * (1 + rt_p) / (1 - rb_p * rt_p)
    Yp = (1 - rb_p) * (1 - rt_p) / (1 - rb_p * rt_p)
    Xs = (1 + rb_s) * (1 + rt_s) / (1 - rb_s * rt_s)
    cv, ch = _weights(S)
    Kpv = 0.75 * (ne / no) * u ** 2 * np.imag(Xp)
    Kph = -3.0 / (6 * (no / ne) ** 2 + 2) * v ** 2 * np.imag(Yp)
    Ksh = 3.0 / (2 * (ne / no) ** 2 + 6) * np.imag(Xs)
    return np.nan_to_num(cv * Kpv + ch * (Kph + Ksh),
                         nan=0.0, posinf=0.0, neginf=0.0)


def _simpson(f, a, b, n):
    if n % 2:
        n += 1
    x = np.linspace(a, b, n + 1)
    y = f(x)
    w = np.ones(n + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return np.sum(w * y, axis=-1) * (b - a) / (3.0 * n)


def solve(S, u_max=3.0, npts=12000):
    """Five-channel power budget; the fractions sum to 1 by construction."""
    no, ne = S.no_e, S.ne_e
    ns = S.n_sub.real
    half = np.pi / 2 - EDGE
    th_air = np.arcsin(min(1.0 / ne, ns / ne, 1.0))
    rat = ns / ne
    th_c = np.arcsin(rat) if rat < 1.0 else half

    def fA(th, need_air=True):
        return np.stack(integrands_theta(th, S, need_air))

    I1 = _simpson(fA, 0.0, th_air, npts)
    I2 = (_simpson(lambda t: fA(t, False), th_air, th_c, npts)
          if th_c > th_air else np.zeros(3))
    I3 = (_simpson(lambda t: fA(t, False), th_c, half, npts)
          if half - th_c > 1e-9 else np.zeros(3))
    I4 = _simpson(lambda v: integrand_v(v, S), EDGE,
                  np.sqrt(u_max ** 2 - 1.0), npts)

    air = I1[2] + I2[2]
    sub_tot = I1[1] + I2[1]
    w_cone = I1[0] + I2[0]
    P = w_cone + I3[0] + I4
    return dict(air=air / P, sub=(sub_tot - air) / P, wg=I3[0] / P,
                spp=I4 / P, abs=(w_cone - sub_tot) / P, P_tot=P)


def critical_indices(S):
    """k_x/k0 values at which some layer's k_z vanishes, so that its Fresnel
    coefficients degenerate to +-1 and the cavity factor becomes 0/0.  The
    integrand has a finite limit there; it just cannot be evaluated on the
    point itself."""
    v = set()
    for n in S.bot_no + S.bot_ne + S.top_no + S.top_ne:
        c = complex(n)
        if abs(c.imag) < 1e-6 and c.real > 0:
            v.add(round(c.real, 9))
    return sorted(v)


def spectrum(S, n_eff, u_max=3.0, tol=2e-3):
    """Dissipated power density per unit in-plane wavevector k_x/k0,
    resolved into TM and TE.  Returns (dP/dn_eff)_TM, (dP/dn_eff)_TE, both
    normalised so that the total dissipated power integrates to 1."""
    no, ne = S.no_e, S.ne_e
    cv, ch = _weights(S)
    n_eff = np.asarray(n_eff, dtype=float).copy()
    for nc in critical_indices(S):          # nudge off the degenerate points
        bad = np.abs(n_eff - nc) < tol
        n_eff[bad] = nc + np.where(n_eff[bad] >= nc, tol, -tol)
    out = []
    for pol, nref in (('p', ne), ('s', no)):
        u = np.asarray(n_eff, dtype=float) / nref
        cE = np.sqrt(1.0 - u ** 2 + 0j)
        P = _prop(u, S)
        rb, _, rt = P[pol]
        if pol == 'p':
            Xp = (1 + rb) * (1 + rt) / (1 - rb * rt)
            Yp = (1 - rb) * (1 - rt) / (1 - rb * rt)
            Kpv = 0.75 * (ne / no) * np.real(u ** 2 / cE * Xp)
            Kph = 3.0 / (6 * (no / ne) ** 2 + 2) * np.real(cE * Yp)
            dens = u * (cv * Kpv + ch * Kph)
        else:
            Xs = (1 + rb) * (1 + rt) / (1 - rb * rt)
            Ksh = 3.0 / (2 * (ne / no) ** 2 + 6) * np.real(Xs / cE)
            dens = u * ch * Ksh
        out.append(np.nan_to_num(dens / nref, nan=0.0, posinf=0.0, neginf=0.0))
    tot = solve(S, u_max=u_max)['P_tot']
    return out[0] / tot, out[1] / tot

"""
Dipole CPS / power-dissipation solver for the planar bottom-emitting stack

        air | Ag 100 nm | organic d_org (emitter at centre) | ITO 50 nm | glass n_sub

Faithful re-implementation of Planar_sweep22_preprint.m (W.C. Lee), with the
one change the physics demands: the u-integrals are evaluated with a smooth
substitution + Simpson quadrature instead of a rectangle sum on a uniform
u-grid.  The uniform grid under-samples the guided-mode poles (their width in
u is 1e-5..1e-3, the grid spacing is 1e-3), which is what makes disp_matrix
"jitter" as the organic thickness is swept.

Substitutions that remove every integrable singularity of the integrand:
    u < 1 :  u = sin(th),      du = cos(th) dth   (kills the 1/sqrt(1-u^2) branch point)
    u > 1 :  u = sqrt(1+v^2),  du = v/sqrt(1+v^2) dv

Conventions (verified against the source script term by term):
    n~ = n + ik,  field ~ exp(+i kz z),  Im(kz) >= 0
    r_p = (n2^2 kz1 - n1^2 kz2)/(n2^2 kz1 + n1^2 kz2)   ->  r_p=+1, r_s=-1 at a perfect mirror
    t_p = 2 n1 n2 kz1 /(n2^2 kz1 + n1^2 kz2)            ->  E-field amplitude
"""
import numpy as np

# ---------------------------------------------------------------- TMM ------

def _kz(n, kx, k0):
    kz = np.sqrt((n * k0) ** 2 - kx ** 2 + 0j)
    return np.where(kz.imag < 0, -kz, kz)


def stack_rt(ns, ds, kx, k0):
    """r, t of a stack referenced to a source plane inside layer 0.

    ns[0]  : medium containing the source          ds[0]  : source -> first interface
    ns[-1] : exit medium (semi-infinite)           ds[-1] : ignored
    """
    L = len(ns)
    kz = [_kz(np.asarray(n, dtype=complex), kx, k0) for n in ns]
    n2 = [np.asarray(n, dtype=complex) ** 2 for n in ns]

    def fres(j):
        a, b = j, j + 1
        ds_ = kz[a] + kz[b]
        rs = (kz[a] - kz[b]) / ds_
        ts = 2 * kz[a] / ds_
        dp = n2[b] * kz[a] + n2[a] * kz[b]
        rp = (n2[b] * kz[a] - n2[a] * kz[b]) / dp
        tp = 2 * ns[a] * ns[b] * kz[a] / dp
        return rs, ts, rp, tp

    rs, ts, rp, tp = fres(L - 2)
    for j in range(L - 3, -1, -1):
        ph = np.exp(1j * kz[j + 1] * ds[j + 1])
        ph2 = ph * ph
        rsj, tsj, rpj, tpj = fres(j)
        den = 1.0 + rsj * rs * ph2
        ts = tsj * ts * ph / den
        rs = (rsj + rs * ph2) / den
        den = 1.0 + rpj * rp * ph2
        tp = tpj * tp * ph / den
        rp = (rpj + rp * ph2) / den

    ph0 = np.exp(1j * kz[0] * ds[0])
    return rs * ph0 * ph0, ts * ph0, rp * ph0 * ph0, tp * ph0


# ------------------------------------------------------- integrand K(u) ----

def outcoupled(u, d_org, n_ag, n_org, n_ito, d_ito, n_sub, d_ag, lam):
    """Incoherent substrate/air factor  T_sub-air / (1 - R_oled * R_sub-air)."""
    k0 = 2 * np.pi / lam
    kx = n_org * u * k0
    # whole OLED seen from inside the substrate: sub | ITO | organic | Ag | air
    rs_o, _, rp_o, _ = stack_rt([n_sub, n_ito, n_org, n_ag, 1.0],
                                [0.0, d_ito, d_org, d_ag, 0.0], kx, k0)
    R_p_o, R_s_o = np.abs(rp_o) ** 2, np.abs(rs_o) ** 2

    cS = np.sqrt(1.0 - (n_org * u / n_sub) ** 2 + 0j)     # cos(theta) in the substrate
    cA = np.sqrt(1.0 - (n_org * u) ** 2 + 0j)             # cos(theta) in air
    r_p = (1.0 * cS - n_sub * cA) / (1.0 * cS + n_sub * cA)
    r_s = (n_sub * cS - 1.0 * cA) / (n_sub * cS + 1.0 * cA)
    R_p_sa, R_s_sa = np.abs(r_p) ** 2, np.abs(r_s) ** 2
    T_p, T_s = 1.0 - R_p_sa, 1.0 - R_s_sa
    return (T_p / (1.0 - R_p_o * R_p_sa),
            T_s / (1.0 - R_s_o * R_s_sa))


# ------------------------------------------------- Jacobian-weighted forms --
# All 1/cos(theta) and 1/sqrt(u^2-1) factors are cancelled analytically against
# the substitution Jacobian, so nothing is ever divided by zero.

def _prop(u, d_org, st):
    """Common reflection/transmission products at the emitter plane."""
    lam, n_ag, n_org, n_ito, d_ito, n_sub, d_ag = st
    k0 = 2 * np.pi / lam
    kx = n_org * u * k0
    z0 = d_org / 2.0
    rs_b, ts_b, rp_b, tp_b = stack_rt([n_org, n_ito, n_sub],
                                      [d_org - z0, d_ito, 0.0], kx, k0)
    rs_t, ts_t, rp_t, tp_t = stack_rt([n_org, n_ag, 1.0],
                                      [z0, d_ag, 0.0], kx, k0)
    Xp = (1 + rp_b) * (1 + rp_t) / (1 - rp_b * rp_t)      # vertical, p
    Yp = (1 - rp_b) * (1 - rp_t) / (1 - rp_b * rp_t)      # horizontal, p
    Xs = (1 + rs_b) * (1 + rs_t) / (1 - rs_b * rs_t)      # horizontal, s
    Apv = (1 + rp_t) * tp_b / (1 - rp_b * rp_t)   # vertical dipole, p
    Aph = (1 - rp_t) * tp_b / (1 - rp_b * rp_t)   # horizontal dipole, p
    Ash = (1 + rs_t) * ts_b / (1 - rs_b * rs_t)   # horizontal dipole, s
    return Xp, Yp, Xs, Apv, Aph, Ash


def integrands_theta(th, d_org, st, need_air=True):
    """(W, S, A) d u / d th  on the propagating branch  u = sin(th)."""
    lam, n_ag, n_org, n_ito, d_ito, n_sub, d_ag = st
    u, c = np.sin(th), np.cos(th)
    Xp, Yp, Xs, Apv, Aph, Ash = _prop(u, d_org, st)

    # total dissipated power density (already x cos th)
    Kpv = 0.75 * u ** 2 * np.real(Xp)
    Kph = 0.375 * c ** 2 * np.real(Yp)
    Ksh = 0.375 * np.real(Xs)
    W = u * (Kpv / 3.0 + 2.0 / 3.0 * (Kph + Ksh))

    # transmitted into a semi-infinite substrate (zero beyond the substrate cone)
    rat = n_sub / n_org
    inside = u < rat
    cS = np.sqrt(np.clip(rat ** 2 - u ** 2, 0.0, None))
    safe_c = np.where(c > 0, c, 1.0)
    Kpv2 = 0.375 * cS * u ** 2 * np.abs(Apv) ** 2 / safe_c
    Kph2 = 0.1875 * cS * np.abs(Aph) ** 2 * c
    Ksh2 = 0.1875 * cS * np.abs(Ash) ** 2 / safe_c
    if abs(rat - 1.0) < 1e-12:            # cS == c exactly: cancel analytically
        Kpv2 = 0.375 * u ** 2 * np.abs(Apv) ** 2
        Kph2 = 0.1875 * np.abs(Aph) ** 2 * c ** 2
        Ksh2 = 0.1875 * np.abs(Ash) ** 2
    S = u * (Kpv2 / 3.0 + 2.0 / 3.0 * (Kph2 + Ksh2)) * inside

    # escaped into air through an incoherent substrate
    if need_air:
        with np.errstate(divide='ignore', invalid='ignore'):
            Fp, Fs = outcoupled(u, d_org, n_ag, n_org, n_ito, d_ito,
                                n_sub, d_ag, lam)
        esc = u <= 1.0 / n_org
        Fp = np.where(esc & np.isfinite(Fp), Fp, 0.0)
        Fs = np.where(esc & np.isfinite(Fs), Fs, 0.0)
        A = u * ((Kpv2 / 3.0 + 2.0 / 3.0 * Kph2) * Fp
                 + 2.0 / 3.0 * Ksh2 * Fs) * inside
    else:
        A = np.zeros_like(W)
    W, S, A = (np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
               for x in (W, S, A))
    return W, S, A


def integrand_v(v, d_org, st):
    """W du/dv on the evanescent branch  u = sqrt(1+v^2)."""
    u = np.sqrt(1.0 + v ** 2)
    Xp, Yp, Xs, Apv, Aph, Ash = _prop(u, d_org, st)
    # each K is already multiplied by v, which is exactly the Jacobian
    # u * du/dv = v that W would otherwise carry
    Kpv = 0.75 * u ** 2 * np.imag(Xp)
    Kph = -0.375 * v ** 2 * np.imag(Yp)
    Ksh = 0.375 * np.imag(Xs)
    return np.nan_to_num(Kpv / 3.0 + 2.0 / 3.0 * (Kph + Ksh),
                         nan=0.0, posinf=0.0, neginf=0.0)


# ----------------------------------------------------------- quadrature ----

def _simpson(f, a, b, n):
    """Composite Simpson on [a,b] with n panels (n even); f takes an array."""
    if n % 2:
        n += 1
    x = np.linspace(a, b, n + 1)
    y = f(x)
    w = np.ones(n + 1)
    w[1:-1:2] = 4.0
    w[2:-1:2] = 2.0
    return np.sum(w * y, axis=-1) * (b - a) / (3.0 * n)


def solve(d_org, n_sub, lam=550.0, n_ag=0.044 + 3.819j, n_org=1.8,
          n_ito=1.86 + 0.003j, d_ito=50.0, d_ag=100.0, u_max=3.0, npts=24000):
    """Five-channel power budget. Returns a dict of fractions summing to 1."""
    st = (lam, n_ag, n_org, n_ito, d_ito, n_sub, d_ag)
    rat = n_sub / n_org
    th_air = np.arcsin(min(1.0 / n_org, rat, 1.0))
    th_c = np.arcsin(min(rat, 1.0)) if rat < 1.0 else np.pi / 2 - 1e-5
    # u = 1 exactly makes k_z in the organic vanish and every Fresnel
    # coefficient degenerate to -1 (0/0).  The integrand tends to zero there,
    # so the branch points are approached to within EDGE instead of touched;
    # the excluded sliver is O(EDGE^2) ~ 1e-10.
    EDGE = 1e-5
    half = np.pi / 2 - EDGE

    def fA(th, need_air=True):
        W, S, A = integrands_theta(th, d_org, st, need_air)
        return np.stack([W, S, A])

    I1 = _simpson(fA, 0.0, th_air, npts)                       # 0 .. escape cone
    I2 = (_simpson(lambda t: fA(t, False), th_air, th_c, npts)
          if th_c > th_air else np.zeros(3))
    I3 = (_simpson(lambda t: fA(t, False), th_c, half, npts)
          if half - th_c > 1e-9 else np.zeros(3))
    v_max = np.sqrt(u_max ** 2 - 1.0)
    I4 = _simpson(lambda v: integrand_v(v, d_org, st), EDGE, v_max, npts)

    air = I1[2] + I2[2]
    sub_tot = I1[1] + I2[1]
    w_cone = I1[0] + I2[0]
    P_tot = w_cone + I3[0] + I4
    return dict(d_org=d_org, n_sub=n_sub,
                air=air / P_tot,
                sub=(sub_tot - air) / P_tot,
                wg=I3[0] / P_tot,
                spp=I4 / P_tot,
                abs=(w_cone - sub_tot) / P_tot,
                P_tot=P_tot)


# ------------------------------------------------------------- legacy ------

def legacy(d_org, n_sub, lam=550.0, n_ag=0.044 + 3.819j, n_org=1.8,
           n_ito=1.86 + 0.003j, d_ito=50.0, d_ag=100.0, u_max=3, N=1000):
    """Bit-for-bit reproduction of Planar_sweep22_preprint.m: uniform u-grid,
    rectangle sum, index-snapped channel boundaries.  Provided only so that the
    new quadrature can be checked against it (and its jitter demonstrated)."""
    st = (lam, n_ag, n_org, n_ito, d_ito, n_sub, d_ag)
    u = np.concatenate([np.arange(0, N) / N, np.arange(N + 1, u_max * N + 1) / N])

    Xp, Yp, Xs, Apv, Aph, Ash = _prop(u, d_org, st)
    cE = np.sqrt(1.0 - u ** 2 + 0j)
    Kpv = 0.75 * np.real(u ** 2 / cE * Xp)
    Kph = 0.375 * np.real(cE * Yp)
    Ksh = 0.375 * np.real(Xs / cE)
    W = u * (Kpv / 3.0 + 2.0 / 3.0 * (Kph + Ksh))

    rat = n_sub / n_org
    cS = np.sqrt(np.clip(rat ** 2 - u ** 2, 0.0, None))
    one_u2 = np.abs(1.0 - u ** 2)
    with np.errstate(divide='ignore', invalid='ignore'):
        Kpv2 = 0.375 * cS * u ** 2 * np.abs(Apv) ** 2 / one_u2
        Kph2 = 0.1875 * cS * np.abs(Aph) ** 2
        Ksh2 = 0.1875 * cS * np.abs(Ash) ** 2 / one_u2
    Sp = u * (Kpv2 / 3.0 + 2.0 / 3.0 * Kph2)
    Ss = u * (2.0 / 3.0 * Ksh2)
    with np.errstate(divide='ignore', invalid='ignore'):
        Fp, Fs = outcoupled(u, d_org, n_ag, n_org, n_ito, d_ito, n_sub, d_ag, lam)
    Fp = np.where(np.isfinite(Fp), Fp, 0.0)
    Fs = np.where(np.isfinite(Fs), Fs, 0.0)
    for a in (W, Sp, Ss):
        a[~np.isfinite(a)] = 0.0

    # MATLAB index-snapped boundaries (1-based counts -> python slice ends)
    m = -1 if n_sub > n_org else 0
    i_sub = int(np.ceil(N * n_sub / n_org)) + m
    i_air = min(i_sub, int(np.ceil(N * 1.0 / n_org)))

    air = np.sum(Sp[:i_air] * Fp[:i_air]) + np.sum(Ss[:i_air] * Fs[:i_air])
    sub_tot = np.sum(Sp[:i_sub]) + np.sum(Ss[:i_sub])
    sub = sub_tot - air
    wg = np.sum(W[i_sub:N])
    spp = np.sum(W[N:])
    ab = np.sum(W[:i_sub]) - sub_tot

    P = np.sum(W)
    air, sub = air / P, sub / P
    rem = 1.0 - (air + sub)
    tl = (wg + spp + ab) / P
    if tl > 0:
        wg, spp, ab = (rem * wg / P / tl, rem * spp / P / tl, rem * ab / P / tl)
    else:
        wg = spp = ab = 0.0
    return dict(d_org=d_org, n_sub=n_sub, air=air, sub=sub, wg=wg, spp=spp,
                abs=ab, closure=tl / (rem if rem else 1.0))

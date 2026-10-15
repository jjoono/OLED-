"""Forward models shared by every stage of the re-fit.

Both substrates are described with the same corrected transfer matrix
(lib/tmm_fix.tmm - propagation exponent e^{-i k0 q d} for the forward wave):

  Si  piece   air / roughness (Bruggeman 50 % Ag + void) / Ag / seed /
              native oxide / Si          -> Psi, Delta -> N, C, S
  glass piece air / roughness / Ag / seed / glass(1 mm) / air
              -> coherent front stack on a thick incoherent slab -> T, R

The glass slab is 1 mm, thousands of wavelengths thick, so its two surfaces add
in intensity, not amplitude: the front stack is coherent and the back surface
enters through the usual multiple-reflection sum.
"""
import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'lib'))

import numpy as np
import ce_osc as osc
from tmm_fix import tmm

HC = 1239.841984      # eV.nm


# ------------------------------------------------------------------ optics
def bruggeman(N_a, N_b, f_a=0.5):
    """Bruggeman EMA, volume fraction f_a of N_a; physical root has Im(N) >= 0."""
    ea, eb = N_a**2, N_b**2
    fb = 1.0 - f_a
    b = (2 * f_a - fb) * ea + (2 * fb - f_a) * eb
    s = np.sqrt((b**2 + 8 * ea * eb).astype(complex))
    N1 = np.sqrt((b + s) / 4.0); N1 = np.where(N1.imag < 0, -N1, N1)
    N2 = np.sqrt((b - s) / 4.0); N2 = np.where(N2.imag < 0, -N2, N2)
    return np.where(N1.real >= N2.real, N1, N2)


def osc_N(einf, oscs, wl):
    """oscs = [('drude', rho, tau) | ('gauss', A, Br, En) | ('tl', A, Br, Eo, Eg)]"""
    E = HC / np.asarray(wl, float)
    e1 = np.zeros_like(E); e2 = np.zeros_like(E)
    for o in oscs:
        if o[0] == 'drude':
            a, b = osc.drude_rt(E, o[1], o[2])
        elif o[0] == 'gauss':
            a, b = osc.gaussian(E, o[1], o[2], o[3])
        elif o[0] == 'tl':
            a, b = osc.tauc_lorentz(E, o[1], o[2], o[3], o[4])
        else:
            raise ValueError(o[0])
        e1 += a; e2 += b
    N = np.sqrt((einf + e1 + 1j * e2).astype(complex))
    return np.where(N.imag < 0, -N, N)


def t_amp(wl, layers, ds, theta_deg):
    """Transmission amplitudes in the same convention as tmm_fix.tmm."""
    th = np.deg2rad(theta_deg); k0 = 2 * np.pi / wl
    s0 = layers[0] * np.sin(th); n = len(wl)
    q = []
    for ni in layers:
        c = np.sqrt((ni**2 - s0**2).astype(complex))
        q.append(np.where(c.imag < 0, -c, c))
    def eye():
        M = np.zeros((n, 2, 2), complex); M[:, 0, 0] = 1; M[:, 1, 1] = 1; return M
    mm = lambda A, B: np.einsum('nij,njk->nik', A, B)
    def I_(r, t):
        M = np.zeros((n, 2, 2), complex)
        M[:, 0, 0] = 1 / t; M[:, 0, 1] = r / t; M[:, 1, 0] = r / t; M[:, 1, 1] = 1 / t
        return M
    def Ip(ni, nj, qi, qj):
        d = nj**2 * qi + ni**2 * qj
        return I_((nj**2 * qi - ni**2 * qj) / d, 2 * ni * nj * qi / d)
    def Is(qi, qj):
        d = qi + qj
        return I_((qi - qj) / d, 2 * qi / d)
    def Pr(qj, d):
        b = k0 * qj * d
        M = np.zeros((n, 2, 2), complex)
        M[:, 0, 0] = np.exp(-1j * b); M[:, 1, 1] = np.exp(1j * b)
        return M
    Mp, Ms = eye(), eye()
    for i, d in enumerate(ds):
        ni, nj, qi, qj = layers[i], layers[i + 1], q[i], q[i + 1]
        Mp = mm(mm(Mp, Ip(ni, nj, qi, qj)), Pr(qj, d))
        Ms = mm(mm(Ms, Is(qi, qj)), Pr(qj, d))
    ni, nj, qi, qj = layers[-2], layers[-1], q[-2], q[-1]
    Mp = mm(Mp, Ip(ni, nj, qi, qj)); Ms = mm(Ms, Is(qi, qj))
    return 1 / Mp[:, 0, 0], 1 / Ms[:, 0, 0]


# ------------------------------------------------------------------ SE side
def se_ncs(wl, layers, ds, angles):
    """-> N, C, S each (nwl, nang).  layers/ds are ambient-first."""
    N = np.empty((len(wl), len(angles)))
    C = np.empty_like(N); S = np.empty_like(N)
    for a, ang in enumerate(angles):
        rp, rs = tmm(wl, layers, ds, ang)
        rho = np.conj(rp / rs)          # instrument's Delta sign convention
        psi = np.arctan(np.abs(rho)); dl = np.angle(rho)
        N[:, a] = np.cos(2 * psi)
        C[:, a] = np.sin(2 * psi) * np.cos(dl)
        S[:, a] = np.sin(2 * psi) * np.sin(dl)
    return N, C, S


def mse_ncs(res, npar):
    """CompleteEASE-like MSE over an N/C/S residual vector."""
    n = res.size // 3
    return 1000.0 * np.sqrt(np.sum(res**2) / max(3 * n - npar, 1))


# --------------------------------------------------------------- T/R side
def stack_TR(wl, layers, ds, n_glass, theta_deg=0.0, tau=1.0):
    """Coherent front stack (ambient .. glass) on a thick incoherent slab.

    layers/ds describe the film stack only; ambient and glass are added here.
    tau is the slab's single-pass internal transmittance exp(-4 pi k d / lambda):
    1 mm of soda-lime is not quite transparent and the bare-substrate scan
    measures how far from 1 it is.
    """
    amb = np.ones(len(wl), complex)
    ng = np.asarray(n_glass, complex)
    front = [amb] + list(layers) + [ng]
    back = [ng] + list(layers)[::-1] + [amb]
    dsb = list(ds)[::-1]

    rp, rs = tmm(wl, front, ds, theta_deg)
    Rf = (np.abs(rp)**2 + np.abs(rs)**2) / 2
    tp, ts = t_amp(wl, front, ds, theta_deg)
    Tf = (np.abs(tp)**2 + np.abs(ts)**2) / 2 * ng.real

    rp2, rs2 = tmm(wl, back, dsb, theta_deg)
    Rb = (np.abs(rp2)**2 + np.abs(rs2)**2) / 2
    tp2, ts2 = t_amp(wl, back, dsb, theta_deg)
    Tb = (np.abs(tp2)**2 + np.abs(ts2)**2) / 2 / ng.real

    Rg = ((ng.real - 1) / (ng.real + 1))**2       # glass/air back surface
    tau = np.asarray(tau, float)
    den = 1.0 - Rb * Rg * tau**2
    return Tf * tau * (1 - Rg) / den, Rf + Tf * Tb * Rg * tau**2 / den


def n_from_glass_T(T):
    """Index of a transparent slab from its measured total transmittance.

    T = (1-R0)/(1+R0) with R0 = ((n-1)/(n+1))^2 - both surfaces, incoherent.
    """
    R0 = (1.0 - T) / (1.0 + T)
    s = np.sqrt(np.clip(R0, 0, 1))
    return (1 + s) / (1 - s)



if __name__ == '__main__':
    wl = np.array([550.0])
    ng = np.array([1.5230 + 0j])
    T, R = stack_TR(wl, [], [], ng)
    n_an = n_from_glass_T(T)
    print('bare slab n = 1.5230 ->  T = %.4f  R = %.4f  T+R = %.4f' % (T[0], R[0], T[0] + R[0]))
    print('inverted back from T  ->  n = %.4f' % n_an[0])
    Ag = np.array([0.13 + 3.90j]); air = np.array([1 + 0j])
    rp, rs = tmm(np.array([600.0]), [air, Ag, air], [200.0], 0.0)
    tp, ts = t_amp(np.array([600.0]), [air, Ag, air], [200.0], 0.0)
    print('200 nm Ag in air, 600 nm:  R = %.4f  T = %.2e  R+T = %.4f (<=1)'
          % (abs(rp[0])**2, abs(tp[0])**2, abs(rp[0])**2 + abs(tp[0])**2))

"""Thin-film collimator for top-emitting micro-OLEDs:  device | randomizing scatterer | DBR | air.

Light leaving the device is randomized by a thin scatterer, meets an angle-selective DBR
that passes only a near-normal cone, and the rejected light is returned through the
scatterer to the low-loss device reflector for another attempt (multi-pass extraction,
Kim, Lee, Yoo, Gather, preprint).  All media above the device have n = N_MED.

Device optics: CPS dipole model + TMM of the preprint (cps2), solved per wavelength.
DBR: isotropic TMM, TE/TM averaged, incidence medium N_MED, exit air.
Multi-pass: angle-binned state vector, scatterer = partial randomizer (m) with loss a_sc
per traversal, device reflectance R(theta) per wavelength.
"""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'cps'))
import cps2, materials as M

N_MED = 1.5
NB = 360                                    # 0.25 deg bins in the medium
TH = (np.arange(NB) + 0.5) * 90.0 / NB
THR = np.radians(TH)
W = np.cos(THR) * np.sin(THR); WN = W / W.sum()
AIR_EDGES = np.linspace(0, 90, 361)         # 0.25 deg air bins
TH_AIR = np.degrees(np.arcsin(np.clip(N_MED * np.sin(THR), 0, 1)))
ESC = N_MED * np.sin(THR) < 1.0


# ---------------- device ----------------
def device(lam, d_exit, d_refl, h=0.8, cathode='ito', d_ag=15.0, d_cpl=60.0):
    """Top emitter. EML 20 nm, dipole at centre. Reflective Ag anode on the far side.
    cathode='ito'   : 50 nm ITO-type transparent cathode (weak cavity)
    cathode='ag'    : thin Ag (d_ag) + capping layer n=1.9 (d_cpl)  (conventional strong cavity)"""
    above = [(1.8, 1.8, d_refl), (M.AG, M.AG, 100.0)]
    if cathode == 'ito':
        below = [(1.8, 1.8, d_exit), (M.ITO, M.ITO, 50.0)]
    else:
        below = [(1.8, 1.8, d_exit), (M.AG, M.AG, d_ag), (1.9, 1.9, d_cpl)]
    S = cps2.Stack(float(lam), (1.8, 1.8), 20.0, 10.0, above=above, below=below, n_sub=N_MED, h=h)
    r = cps2.solve_pol(S, npts=3000)
    R = np.clip(np.real(cps2.stack_reflectance(S, TH)), 0, 1)
    P = np.clip(np.real(cps2.sub_angular(S, TH)), 0, None)
    return float(r['air'] + r['sub']), P / P.sum(), R


# ---------------- DBR ----------------
def tmm_T(ns, ds, lam, th_in, n_in=N_MED, n_out=1.0):
    """Unpolarized power transmittance of a lossless multilayer (ns, ds in nm), from n_in at
    angle th_in (rad array) into n_out."""
    kx = n_in * np.sin(th_in)
    out = np.zeros_like(th_in)
    for pol in ('s', 'p'):
        def eta(n):
            c = np.sqrt((n ** 2 - kx ** 2).astype(complex))
            return c if pol == 's' else n ** 2 / c
        Mtot = np.tile(np.eye(2, dtype=complex), (th_in.size, 1, 1))
        for n, d in zip(ns, ds):
            c = np.sqrt((n ** 2 - kx ** 2).astype(complex)); delta = 2 * np.pi / lam * c * d; e = eta(n)
            m = np.empty((th_in.size, 2, 2), complex)
            m[:, 0, 0] = np.cos(delta); m[:, 0, 1] = 1j * np.sin(delta) / e
            m[:, 1, 0] = 1j * e * np.sin(delta); m[:, 1, 1] = np.cos(delta)
            Mtot = Mtot @ m
        e0, es = eta(n_in), eta(n_out)
        B = Mtot[:, 0, 0] + Mtot[:, 0, 1] * es
        Cc = Mtot[:, 1, 0] + Mtot[:, 1, 1] * es
        t = 2 * e0 / (e0 * B + Cc)
        T = np.real(es) / np.real(e0) * np.abs(t) ** 2
        T[np.imag(np.sqrt((n_out ** 2 - kx ** 2).astype(complex))) > 1e-9] = 0.0   # TIR
        out += 0.5 * np.nan_to_num(T)
    return np.clip(out, 0, 1)


def dbr_stack(lam_b, N, nH=2.40, nL=1.46):
    """(L/2 H L/2)^N symmetric period, quarter-wave at lam_b: an edge filter whose stop band
    sits on the long-wavelength side of the emitter at normal incidence and blue-shifts onto
    it at oblique incidence."""
    dH, dL = lam_b / (4 * nH), lam_b / (4 * nL)
    ns, ds = [], []
    for _ in range(N):
        ns += [nL, nH, nL]; ds += [dL / 2, dH, dL / 2]
    return ns, ds


# ---------------- multi-pass ----------------
def multipass(P, R, T, m=1.0, a_sc=0.0, K=400):
    """Returns air-angle power histogram (0.25 deg bins) per unit power delivered to the medium."""
    def scat(v):
        return (1 - a_sc) * ((1 - m) * v + m * v.sum() * WN)
    v = scat(P.copy())
    air = np.zeros(NB)
    for _ in range(K):
        if v.sum() < 1e-7:
            break
        air += T * v
        v = scat(R * scat((1 - T) * v))
    hist = np.bincount(np.clip(np.searchsorted(AIR_EDGES, TH_AIR[ESC]) - 1, 0, 359),
                       weights=air[ESC], minlength=360)
    return hist


def planar(P, R, T, K=200):
    """No scatterer: angle-preserving multi-pass (planar cover / micro-OLED reference)."""
    v = P.copy(); air = np.zeros(NB)
    for _ in range(K):
        if v.sum() < 1e-7:
            break
        air += T * v
        v = R * (1 - T) * v
    return np.bincount(np.clip(np.searchsorted(AIR_EDGES, TH_AIR[ESC]) - 1, 0, 359),
                       weights=air[ESC], minlength=360)


def fresnel_T():
    return tmm_T([], [], 530.0, THR)        # bare n=1.5 -> air interface


def spectrum(lam, lam0=530.0, fwhm=30.0):
    s = np.exp(-4 * np.log(2) * ((lam - lam0) / fwhm) ** 2)
    return s / s.sum()


# ---------------- metrics ----------------
AIR_C = 0.5 * (AIR_EDGES[:-1] + AIR_EDGES[1:])


def cone(hist, beta):
    return hist[AIR_C < beta].sum()


def axial_radiance(hist, half=2.5):
    """Power per steradian within `half` deg of the normal (normalized so a Lambertian source of
    unit total power gives 1/pi)."""
    sel = AIR_C < half
    sr = 2 * np.pi * (1 - np.cos(np.radians(half)))
    return hist[sel].sum() / sr

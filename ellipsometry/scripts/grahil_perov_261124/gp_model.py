"""Model library for the GraHIL / perovskite-NC ellipsometry (gp_*.py).

Optics: corrected TMM (lib/tmm_fix), instrument Delta convention conj(rp/rs)
as in the earlier sessions.  Si and its native oxide are the CompleteEASE
SI_JAW / NTVE_JAW tables (decoded from the lab .mod, 188 nm upward).

What is here
  KKBasis      a Kramers-Kronig-consistent, model-free dispersion: eps2 is a
               non-negative sum of Gaussians on a fixed photon-energy grid,
               eps1 their exact KK partner + a UV pole (+ eps_inf).  Linear in
               the amplitudes, so cheap; the grid can be dense where the
               material has sharp structure (a band edge).
  ncs()        N, C, S of a stack.  Optional lateral thickness non-uniformity
               (the stack is averaged incoherently over a Gaussian spread of
               the film thicknesses), which is what depolarizes: the result
               then has N^2 + C^2 + S^2 < 1 and is compared with the measured
               N, C, S scaled by sqrt(1 - depolarization).
  ncs_coh()    the same spread but INSIDE the coherence area: amplitudes are
               averaged, which damps interference fringes without depolarizing
               (roughness with a lateral scale that is not small against the
               wavelength - grains, film undulation).
  graded()     a layer cut into sublayers whose composition runs linearly from
               bottom to top (Bruggeman mix of two media)."""
import os as _os
import sys, numpy as np
_HERE = _os.path.dirname(_os.path.abspath(__file__))
for p in ('../tooling_261117', '../joint_nk_260910', '../lib'):
    sys.path.insert(0, _os.path.join(_HERE, p))
import tl_ce_genosc as C
import jnk_model as M, ce_osc as osc
from tmm_fix import tmm

HC = 1239.841984
bruggeman = M.bruggeman
mse_ncs = M.mse_ncs


def n_si(wl):   return C.N_ce('SI_JAW', wl)
def n_ox(wl):   return C.N_ce('NTVE_JAW', wl)
def n_air(wl):  return np.ones(len(wl), complex)


class KKBasis:
    """eps(E) = einf + A_uv / (E_uv^2 - E^2) + sum_j a_j G_j(E),  a_j >= 0.

    G_j: CompleteEASE Gaussian of unit amplitude at centre c_j, FWHM w_j
    (eps2 and its KK eps1).  centres / widths are fixed; parameters are
    [einf?, A_uv, a_1..a_n]."""
    def __init__(self, wl, centres, widths, e_uv=11.0):
        self.wl = np.asarray(wl, float); self.E = HC / self.wl
        self.c = np.asarray(centres, float); self.w = np.broadcast_to(np.asarray(widths, float), self.c.shape).copy()
        self.e_uv = e_uv
        G1 = np.empty((len(self.E), len(self.c))); G2 = np.empty_like(G1)
        for j, (cj, wj) in enumerate(zip(self.c, self.w)):
            G1[:, j], G2[:, j] = osc.gaussian(self.E, 1.0, wj, cj)
        self.G = G1 + 1j * G2
        self.pole = 1.0 / (e_uv**2 - self.E**2)
        self.n = len(self.c)

    def eps(self, a, a_uv, einf=1.0):
        return einf + a_uv * self.pole + self.G @ a

    def N(self, a, a_uv, einf=1.0):
        n = np.sqrt(self.eps(a, a_uv, einf)); return np.where(n.imag < 0, -n, n)

    def on(self, wl):
        """the same basis evaluated on another wavelength grid"""
        return KKBasis(wl, self.c, self.w, self.e_uv)


def grid(lo, hi, step, fine=None):
    """Gaussian centres from lo to hi eV; fine = (lo, hi, step) inserts a denser stretch.
    Widths = 1.6 x the local spacing, so neighbours overlap into a smooth eps2."""
    c = list(np.arange(lo, hi + 1e-9, step))
    if fine:
        c = [x for x in c if not (fine[0] - step / 2 < x < fine[1] + step / 2)] + list(np.arange(fine[0], fine[1] + 1e-9, fine[2]))
    c = np.array(sorted(c)); sp = np.gradient(c)
    return c, 1.6 * sp


def graded(Na, Nb, f_bot, f_top, d, m):
    """m sublayers, top first; volume fraction of Na goes f_top -> f_bot. Returns (layers, thicknesses)."""
    L, D = [], []
    for i in range(m):
        f = f_top + (f_bot - f_top) * (i + 0.5) / m
        L.append(bruggeman(Na, Nb, f)); D.append(d / m)
    return L, D


_GH = {1: (np.array([0.0]), np.array([1.0]))}
for _k in (3, 5, 7):
    _x, _w = np.polynomial.hermite_e.hermegauss(_k); _GH[_k] = (_x, _w / _w.sum())


def ncs(wl, ang, layers, ds, sig=0.0, vary=None, nodes=5):
    """N, C, S [n_wl, n_ang].  layers: complex index arrays, ambient .. substrate.
    sig: relative 1-sigma spread of the thicknesses listed in vary (indices into ds; default all)."""
    x, w = _GH[nodes if sig > 0 else 1]
    vary = range(len(ds)) if vary is None else vary
    out = []
    for th in ang:
        Ip = Is = 0.0; X = 0.0
        for xi, wi in zip(x, w):
            dd = list(ds)
            for j in vary: dd[j] = max(ds[j] * (1 + sig * xi), 0.0)
            rp, rs = tmm(wl, layers, dd, th)
            Ip = Ip + wi * np.abs(rp)**2; Is = Is + wi * np.abs(rs)**2; X = X + wi * np.conj(rp) * rs
        I = Ip + Is
        out.append(((Is - Ip) / I, 2 * X.real / I, 2 * X.imag / I))
    return tuple(np.stack([o[k] for o in out], 1) for k in range(3))


def ncs_coh(wl, ang, layers, ds, j, sig, nodes=11, extra=None):
    """N, C, S with the height of the top surface spread laterally INSIDE the coherence area.

    The thickness ds[j] (the top film; j indexes ds) takes the values d + h, h ~ N(0, sig), and the
    AMPLITUDES r_p, r_s are averaged, each referred to one fixed plane above the sample (an air gap
    c - h is put on top, so the phase of the front reflection moves with the surface as it does in
    reality).  This is the scalar (Kirchhoff) picture of roughness whose lateral scale is not small
    against the wavelength: it damps the interference between the front and the buried interfaces
    and sends the lost light into diffuse scatter, but leaves the specular beam fully polarized - no
    depolarization, unlike an incoherent average (ncs(..., sig=)).  On an opaque film it does nothing
    to rho = r_p / r_s, since both get the same factor.
    extra: optional (weights, [(layers_k, ds_k), ...]) is not used here; see ncs_mix for area mixing."""
    if sig <= 0:
        return ncs(wl, ang, layers, ds)
    x, w = np.polynomial.hermite_e.hermegauss(nodes); w = w / w.sum()
    c = sig * x.max()
    L = [layers[0], layers[0]] + list(layers[1:])
    out = []
    for th in ang:
        rp = 0.0; rs = 0.0
        for xi, wi in zip(x, w):
            h = sig * xi; dd = list(ds); dd[j] = max(ds[j] + h, 0.0)
            a, b_ = tmm(wl, L, [c - h] + dd, th)
            rp = rp + wi * a; rs = rs + wi * b_
        Ip, Is = np.abs(rp)**2, np.abs(rs)**2; X = np.conj(rp) * rs; I = Ip + Is
        out.append(((Is - Ip) / I, 2 * X.real / I, 2 * X.imag / I))
    return tuple(np.stack([o[k] for o in out], 1) for k in range(3))


def rr_coh2(wl, ang, layers, ds, jt, jb, sig_s, sig_i, nodes=7):
    """Coherently averaged r_p, r_s [n_wl, n_ang] with TWO independent Gaussian height spreads:
    sig_s of the top surface and sig_i of the internal boundary between ds[jt] (upper film) and
    ds[jb] (the film under it).  Surface up by h: upper film + h, air gap - h.  Boundary up by g:
    lower film + g, upper film - g.  Amplitudes are referred to one fixed plane, so this also
    carries the Debye-Waller loss of the specular beam (compare absolute reflectance with it)."""
    xs, ws = (np.array([0.0]), np.array([1.0])) if sig_s <= 0 else np.polynomial.hermite_e.hermegauss(nodes)
    xi_, wi_ = (np.array([0.0]), np.array([1.0])) if sig_i <= 0 else np.polynomial.hermite_e.hermegauss(nodes)
    ws = ws / ws.sum(); wi_ = wi_ / wi_.sum()
    c = sig_s * xs.max()
    L = [layers[0], layers[0]] + list(layers[1:])
    RP, RS = [], []
    for th in ang:
        rp = 0.0; rs = 0.0
        for x1, w1 in zip(xs, ws):
            h = sig_s * x1
            for x2, w2 in zip(xi_, wi_):
                g = sig_i * x2; dd = list(ds)
                dd[jt] = max(ds[jt] + h - g, 0.0); dd[jb] = max(ds[jb] + g, 0.0)
                a_, b_ = tmm(wl, L, [c - h] + dd, th)
                rp = rp + w1 * w2 * a_; rs = rs + w1 * w2 * b_
        RP.append(rp); RS.append(rs)
    return np.stack(RP, 1), np.stack(RS, 1)


def ncs_from_r(rp, rs):
    Ip, Is = np.abs(rp)**2, np.abs(rs)**2; X = np.conj(rp) * rs; I = Ip + Is
    return (Is - Ip) / I, 2 * X.real / I, 2 * X.imag / I


def tmm_uni(wl, eo, ee, ds, theta):
    """Stack of uniaxial layers, optic axis along z (isotropic: eo = ee); eps lists ambient..substrate.
    Same propagation and sign convention as tmm_fix.tmm (copied from tooling_261117/tl_bphen10)."""
    sq = lambda e: (lambda n: np.where(n.imag < 0, -n, n))(np.sqrt(np.asarray(e, complex)))
    th = np.deg2rad(theta); k0 = 2 * np.pi / wl
    s2 = (sq(eo[0]) * np.sin(th))**2
    qs = [sq(o - s2) for o in eo]
    qp = [sq(o * (1 - s2 / e)) for o, e in zip(eo, ee)]
    eta_s = qs; eta_p = [o / q for o, q in zip(eo, qp)]
    n = len(wl)
    def run(eta, q):
        Mx = np.zeros((n, 2, 2), complex); Mx[:, 0, 0] = 1; Mx[:, 1, 1] = 1
        mm = lambda A, B: np.einsum('nij,njk->nik', A, B)
        for i in range(len(eo) - 1):
            r = (eta[i] - eta[i + 1]) / (eta[i] + eta[i + 1])
            I = np.zeros((n, 2, 2), complex); I[:, 0, 0] = 1; I[:, 0, 1] = r; I[:, 1, 0] = r; I[:, 1, 1] = 1
            Mx = mm(Mx, I)
            if i < len(ds):
                ph = k0 * q[i + 1] * ds[i]
                Pm = np.zeros((n, 2, 2), complex); Pm[:, 0, 0] = np.exp(-1j * ph); Pm[:, 1, 1] = np.exp(1j * ph)
                Mx = mm(Mx, Pm)
        return Mx[:, 1, 0] / Mx[:, 0, 0]
    return -run(eta_p, qp), run(eta_s, qs)


def data_ncs(d, with_dep=False):
    """measured N, C, S; with_dep scales them by sqrt(1 - depolarization) (clipped to [0, 1])."""
    if not with_dep:
        return d['N'], d['C'], d['S']
    P = np.sqrt(np.clip(1 - d['dep'] / 100.0, 0.0, 1.0))
    return d['N'] * P, d['C'] * P, d['S'] * P


def resid(model, meas, w=None):
    r = np.stack([model[k] - meas[k] for k in range(3)])          # [3, n_wl, n_ang]
    if w is not None: r = r * w[None, :, None]
    return r.ravel()


def psidelta(N, C, S):
    P = np.sqrt(N**2 + C**2 + S**2)
    return 0.5 * np.degrees(np.arccos(np.clip(N / P, -1, 1))), np.degrees(np.arctan2(S, C)) % 360

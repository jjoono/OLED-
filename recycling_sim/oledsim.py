"""
oledsim.py -- wave optics (transfer matrix + classical dipole model) and a
multi-pass recycling model for top-emitting OLEDs with angular photon recycling.

Conventions
-----------
* time dependence exp(-i w t); a "forward" wave in a layer is exp(+i k_z z) with
  Im(k_z) >= 0 (decaying / outgoing).
* beta = k_x / k0 : in-plane wavevector normalised to the vacuum wavenumber.
* All amplitudes are TANGENTIAL electric-field amplitudes, for both s and p.
  Normalised admittances (H_t = eta * E_t for a forward wave):
      s (TE):  eta = q,            q = k_z/k0 = sqrt(no^2 - beta^2)
      p (TM):  eta = no^2 / q,     q = k_z/k0 = sqrt(no^2 - (no/ne)^2 beta^2)
  (uniaxial layers, optic axis along the surface normal; isotropic: no = ne).
  z-flux of a forward wave with tangential amplitude a:  S = 0.5 Re(eta) |a|^2.
* Interface (j -> k):  r = (eta_j - eta_k)/(eta_j + eta_k),  t = 2 eta_j/(eta_j + eta_k).
  NOTE: for p polarisation this tangential-E reflection coefficient is the
  NEGATIVE of the Born&Wolf r_p used in the MATLAB code of this repository
  (TMF_birefringence_whole.m).  The dipole formulas below are written for the
  tangential-E convention, so that  (1 - r_p^BW) = (1 + r_p^tan)  etc.

Dipole model (Neyts 1998; Furno et al. PRB 85, 115205 (2012)), emitter layer
treated as isotropic.  Power densities are per unit of (u^2), u = beta/n_EML,
normalised such that a dipole in an infinite EML medium gives 1 for each
orientation:   F = int 2u du K(u).
"""
from __future__ import annotations

import numpy as np
import scipy.io as sio

# ----------------------------------------------------------------------------
# Material database
# ----------------------------------------------------------------------------
_GRIDS = {301: np.arange(400.0, 701.0), 401: np.arange(400.0, 801.0), 351: np.arange(400.0, 751.0)}


class MaterialDB:
    """Reads nk_JH_total.mat (struct 'material', 'spectrum', 'color').

    Arrays of length 301 are on 400:1:700 nm, 401 on 400:1:800 nm (verified with
    the CIE y_bar(550 nm) = 0.99495 entry and with duplicated spectra).
    """

    def __init__(self, path: str):
        m = sio.loadmat(path, squeeze_me=True, struct_as_record=False)
        self._mat = m["material"]
        self._spec = m["spectrum"]
        self._color = m["color"]

    def names(self):
        return list(self._mat._fieldnames)

    @staticmethod
    def _interp(a, wl):
        a = np.asarray(a).astype(complex).ravel()
        g = _GRIDS[a.size]
        return np.interp(wl, g, a.real) + 1j * np.interp(wl, g, a.imag)

    def nk(self, name, wl):
        """complex refractive index n + ik on wavelength grid wl (nm)."""
        if isinstance(name, (int, float, complex, np.number)):
            return np.full(np.shape(wl), complex(name))
        return self._interp(getattr(self._mat, name), wl)

    def spectrum(self, name, wl):
        s = self._interp(getattr(self._spec, name), wl).real
        return np.clip(s, 0, None)

    def V(self, wl):
        """CIE 1931 photopic luminosity function (= y_bar), 400-700 nm."""
        y = np.asarray(self._color.y_bar).astype(float)
        return np.interp(wl, _GRIDS[301], y, left=0.0, right=0.0)

    def xyz_bar(self, wl):
        g = _GRIDS[301]
        c = self._color
        return (np.interp(wl, g, np.asarray(c.x_bar, float), left=0, right=0),
                np.interp(wl, g, np.asarray(c.y_bar, float), left=0, right=0),
                np.interp(wl, g, np.asarray(c.z_bar, float), left=0, right=0))


def gaussian_spectrum(wl, peak, fwhm, skew=0.0, shoulder=0.0, shoulder_shift=1350.0):
    """Model narrowband emitter (energy spectrum, arbitrary units).

    A Gaussian in photon energy with FWHM given in nm at the peak, plus an optional
    vibronic shoulder red-shifted by `shoulder_shift` cm^-1 with relative height
    `shoulder` (typical MR-TADF: ~0.15-0.3, ~1300-1400 cm^-1)."""
    E = 1e7 / wl  # cm^-1
    E0 = 1e7 / peak
    dE = 1e7 / (peak - fwhm / 2) - 1e7 / (peak + fwhm / 2)
    sig = dE / (2 * np.sqrt(2 * np.log(2)))
    s = np.exp(-0.5 * ((E - E0) / sig) ** 2)
    if shoulder:
        s = s + shoulder * np.exp(-0.5 * ((E - (E0 - shoulder_shift)) / (1.15 * sig)) ** 2)
    return s / s.max()


def energy_to_photon_spectrum(wl, sE):
    """convert an energy spectral density to a photon-number distribution (sum = 1)."""
    sp = np.asarray(sE, float) * wl
    return sp / sp.sum()


# ----------------------------------------------------------------------------
# Transfer-matrix core
# ----------------------------------------------------------------------------
def _col(x):
    """reshape a per-wavelength quantity to (nwl, 1); scalars pass through."""
    x = np.asarray(x)
    return x[:, None] if x.ndim == 1 else x


def _branch(q):
    q = np.where(q.imag < 0, -q, q)
    q = np.where((q.imag == 0) & (q.real < 0), -q, q)
    return q


def kz_norm(no, ne, beta, pol):
    """k_z/k0 for a (uniaxial) medium. no, ne: (nwl,1) or scalar; beta broadcastable."""
    no = _col(no)
    ne = _col(ne)
    if pol == "s":
        return _branch(np.sqrt(no ** 2 - beta ** 2 + 0j))
    return _branch(np.sqrt(no ** 2 - (no / ne) ** 2 * beta ** 2 + 0j))


def admittance(no, q, pol):
    no = _col(no)
    return q if pol == "s" else no ** 2 / q


class Layer:
    """thin-film layer: ordinary/extraordinary complex index arrays (nwl,) and thickness (nm)."""

    def __init__(self, no, d, ne=None, name=""):
        self.no = np.asarray(no, complex)
        self.ne = self.no if ne is None else np.asarray(ne, complex)
        self.d = float(d)
        self.name = name

    def __repr__(self):
        return f"Layer({self.name}, {self.d:g} nm)"


class Medium:
    def __init__(self, no, ne=None, name=""):
        self.no = np.asarray(no, complex)
        self.ne = self.no if ne is None else np.asarray(ne, complex)
        self.name = name


def stack_r_tau(layers, m_in: Medium, m_out: Medium, beta, wl, pol):
    """Stable (Airy/Rouard) recursion.

    Returns (r, tau, eta_in, eta_out, q_in) with shapes (nwl, nb):
      r   : tangential-E reflection coefficient at the first interface (seen from m_in)
      tau : tangential-E amplitude transmitted into m_out per unit forward amplitude
            incident on the first interface.
    Works for evanescent beta (no exponential growth)."""
    wl = np.asarray(wl, float)
    k0 = (2 * np.pi / wl)[:, None]
    B = np.asarray(beta, float)
    if B.ndim == 1:
        B = B[None, :]
    media = [m_in] + list(layers) + [m_out]
    q = [kz_norm(m.no, m.ne, B, pol) for m in media]
    eta = [admittance(m.no, qq, pol) for m, qq in zip(media, q)]
    n = len(media)
    # last interface
    G = (eta[n - 2] - eta[n - 1]) / (eta[n - 2] + eta[n - 1])
    T = 2 * eta[n - 2] / (eta[n - 2] + eta[n - 1])
    for j in range(n - 3, -1, -1):
        lay = media[j + 1]
        ph1 = np.exp(1j * k0 * q[j + 1] * lay.d)
        ph2 = ph1 * ph1
        rj = (eta[j] - eta[j + 1]) / (eta[j] + eta[j + 1])
        tj = 2 * eta[j] / (eta[j] + eta[j + 1])
        den = 1 + rj * G * ph2
        T = tj * ph1 * T / den
        G = (rj + G * ph2) / den
    return G, T, eta[0], eta[-1], q[0]


def stack_RTA(layers, m_in: Medium, m_out: Medium, beta, wl, pol):
    """Reflectance, transmittance and absorptance of each layer for a plane wave
    incident from m_in (must be non-absorbing, beta < n_in).

    Characteristic-matrix back-propagation of tangential (E, H) from the exit side.
    Returns R (nwl,nb), T (nwl,nb), A (nlayers, nwl, nb)."""
    wl = np.asarray(wl, float)
    k0 = (2 * np.pi / wl)[:, None]
    B = np.asarray(beta, float)
    if B.ndim == 1:
        B = B[None, :]
    q_out = kz_norm(m_out.no, m_out.ne, B, pol)
    eta_out = admittance(m_out.no, q_out, pol)
    shape = np.broadcast(k0, B).shape
    E = np.ones(shape, complex)
    H = eta_out * E
    S = [0.5 * np.real(E * np.conj(H))]  # flux at the exit side of each interface, from the back
    for lay in reversed(list(layers)):
        qj = kz_norm(lay.no, lay.ne, B, pol)
        ej = admittance(lay.no, qj, pol)
        phi = k0 * qj * lay.d
        c, s_ = np.cos(phi), np.sin(phi)
        E, H = c * E - 1j * s_ * H / ej, -1j * ej * s_ * E + c * H
        S.append(0.5 * np.real(E * np.conj(H)))
    q_in = kz_norm(m_in.no, m_in.ne, B, pol)
    eta_in = admittance(m_in.no, q_in, pol)
    a_inc = 0.5 * (E + H / eta_in)
    a_ref = 0.5 * (E - H / eta_in)
    S_inc = 0.5 * np.real(eta_in) * np.abs(a_inc) ** 2
    R = np.abs(a_ref / a_inc) ** 2
    S = np.array(S[::-1])  # S[0] = entrance, S[-1] = exit
    T = S[-1] / S_inc
    A = (S[:-1] - S[1:]) / S_inc
    return R, T, A


# ----------------------------------------------------------------------------
# Dipole model
# ----------------------------------------------------------------------------
class DipoleResult(dict):
    __getattr__ = dict.__getitem__


def dipole_emission(wl, n_eml, z_up, z_dn, up_layers, top: Medium, dn_layers, bot: Medium,
                    theta_h=2 / 3, u_max=3.0, du=5e-4, keep_u=True):
    """Classical dipole emission in a planar stack.

    wl        : (nwl,) wavelengths in nm
    n_eml     : (nwl,) complex index of the (isotropic) emitting layer
    z_up/z_dn : distance (nm) from the emitter plane to the upper / lower boundary of the EML
    up_layers : layers from the EML upward (excluding the EML), ending at medium `top`
    dn_layers : layers from the EML downward, ending at medium `bot`
    theta_h   : horizontal dipole ratio (fraction of horizontal dipoles, 2/3 = isotropic)

    Returns a DipoleResult with (per wavelength)
      F        : Purcell factor (total decay-rate enhancement for the orientation mix)
      top, bot : fraction of the radiated power transmitted into the top / bottom media
      and, if keep_u, the u grid and per-polarisation transmitted densities for angular
      distributions (Pt_s, Pt_p, Pb_s, Pb_p) and the total density Ktot.
    """
    wl = np.asarray(wl, float)
    k0 = (2 * np.pi / wl)[:, None]
    ne_ = np.real(np.asarray(n_eml, complex))[:, None]  # source normalisation index
    eml = Medium(np.real(np.asarray(n_eml, complex)) + 1j * np.imag(np.asarray(n_eml, complex)))
    u = (np.arange(int(round(u_max / du))) + 0.5) * du
    U = u[None, :]
    beta = ne_ * U  # (nwl, nu)
    c = np.sqrt(1 - U ** 2 + 0j)
    c = _branch(c)
    qe = ne_ * c  # k_z/k0 in the EML (lossless approximation for the source plane)
    ph_up = np.exp(1j * k0 * qe * z_up)
    ph_dn = np.exp(1j * k0 * qe * z_dn)
    out = DipoleResult()
    dens = {}
    for pol in ("s", "p"):
        ru, tu, _, eta_t, _ = stack_r_tau(up_layers, eml, top, beta, wl, pol)
        rd, td, _, eta_b, _ = stack_r_tau(dn_layers, eml, bot, beta, wl, pol)
        ru = ru * ph_up ** 2
        tu = tu * ph_up
        rd = rd * ph_dn ** 2
        td = td * ph_dn
        den = 1 - ru * rd
        if pol == "s":
            K_h = 3 / 8 * np.real((1 / c) * (1 + ru) * (1 + rd) / den)
            A2 = 3 / 8 / (ne_ * np.abs(c) ** 2)
            Pt = 0.5 * np.real(eta_t) * np.abs(tu) ** 2 * A2 * np.abs(1 + rd) ** 2 / np.abs(den) ** 2
            Pb = 0.5 * np.real(eta_b) * np.abs(td) ** 2 * A2 * np.abs(1 + ru) ** 2 / np.abs(den) ** 2
            dens["s"] = dict(Kh=K_h, Pt_h=Pt, Pb_h=Pb)
        else:
            K_h = 3 / 8 * np.real(c * (1 + ru) * (1 + rd) / den)
            K_v = 3 / 4 * np.real((U ** 2 / c) * (1 - ru) * (1 - rd) / den)
            A2h = 3 / 8 * np.abs(c) ** 2 / ne_
            A2v = 3 / 4 * U ** 2 / ne_
            g_t = 0.5 * np.real(eta_t) * np.abs(tu) ** 2 / np.abs(den) ** 2
            g_b = 0.5 * np.real(eta_b) * np.abs(td) ** 2 / np.abs(den) ** 2
            dens["p"] = dict(Kh=K_h, Kv=K_v,
                             Pt_h=g_t * A2h * np.abs(1 + rd) ** 2, Pb_h=g_b * A2h * np.abs(1 + ru) ** 2,
                             Pt_v=g_t * A2v * np.abs(1 - rd) ** 2, Pb_v=g_b * A2v * np.abs(1 - ru) ** 2)
    th = theta_h
    w = 2 * U * du
    Ktot = (1 - th) * dens["p"]["Kv"] + th * (dens["p"]["Kh"] + dens["s"]["Kh"])
    Pt_s = th * dens["s"]["Pt_h"]
    Pt_p = (1 - th) * dens["p"]["Pt_v"] + th * dens["p"]["Pt_h"]
    Pb_s = th * dens["s"]["Pb_h"]
    Pb_p = (1 - th) * dens["p"]["Pb_v"] + th * dens["p"]["Pb_h"]
    F = np.sum(Ktot * w, axis=1)
    out["F"] = F
    out["top"] = np.sum((Pt_s + Pt_p) * w, axis=1) / F
    out["bot"] = np.sum((Pb_s + Pb_p) * w, axis=1) / F
    out["Ktot_u"] = Ktot / F[:, None]
    if keep_u:
        out["u"] = u
        out["du"] = du
        out["n_eml"] = ne_[:, 0]
        out["Pt_s"], out["Pt_p"] = Pt_s / F[:, None], Pt_p / F[:, None]
        out["Pb_s"], out["Pb_p"] = Pb_s / F[:, None], Pb_p / F[:, None]
    return out


def radiant_intensity(res: DipoleResult, side: str, n_out, theta_deg):
    """Photon radiant intensity per unit solid angle in the outer medium (per unit
    radiated power, per polarisation) at angles theta (deg) measured in that medium.

    I(theta) = P(u(theta)) * (n_out/n_eml)^2 * cos(theta) / pi,  u = (n_out/n_eml) sin(theta).
    Returns (I_s, I_p) arrays (nwl, ntheta)."""
    th = np.deg2rad(np.asarray(theta_deg, float))
    n_out = np.real(_col(n_out))
    ne_ = res["n_eml"][:, None]
    uq = (n_out / ne_) * np.sin(th)[None, :]
    key = "Pt" if side == "top" else "Pb"
    I = []
    for pol in ("s", "p"):
        P = res[f"{key}_{pol}"]
        Pi = np.array([np.interp(uq[i], res["u"], P[i], right=0.0) for i in range(P.shape[0])])
        I.append(Pi * (n_out / ne_) ** 2 * np.cos(th)[None, :] / np.pi)
    return I[0], I[1]


# ----------------------------------------------------------------------------
# Photometry
# ----------------------------------------------------------------------------
HC_E = 1239.841984  # eV*nm, i.e. hc/(e*lambda) in volts when lambda in nm
KM = 683.002  # lm/W


def current_efficiency(wl, V, N_onaxis):
    """cd/A from on-axis photon radiance per electron per sr, N(lambda) summed on wl grid
    (N_onaxis is photons per electron per sr per wavelength bin)."""
    return KM * np.sum(V * (HC_E / wl) * N_onaxis)


def lambertian_ce(wl, V, s_ph, eqe=1.0):
    return KM / np.pi * eqe * np.sum(V * (HC_E / wl) * s_ph)

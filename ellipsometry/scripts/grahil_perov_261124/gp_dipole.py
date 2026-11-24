"""Dipole emission in a planar stack: where the power goes.

A dipole in layer j at depth z radiates into plane waves of in-plane index
s = k_par / k_0 (s = n sin(theta) in any layer).  For each s, polarization and
dipole orientation the stack is solved exactly (source terms + multiple
reflection at the dipole plane, then tangential fields carried through the
layers), and the net power flow through every interface is taken:

    K(s)        power leaving the dipole plane (upward + downward)
    F_first     power entering the first medium (semi-infinite, e.g. glass)
    F_last      power entering the last medium (e.g. the metal)
    A_i         power absorbed in layer i  (difference of the flows at its faces)

normalised so that a dipole in an infinite lossless medium radiates 1.

Conventions
  media 0..M, 0 and M semi-infinite; z increases from medium 0 to medium M.
  U = E_y (TE) or H_y (TM);  V = eta (a - b) with eta = q (TE), q / eps (TM),
  q = sqrt(eps - s^2), Im q >= 0.  U and V are continuous; the power flow
  along +z is Re(U* V).
  Source amplitudes (U waves, upward / downward): the plane-wave expansion of
  the dipole field in the emitter layer (q_j, complex if the layer absorbs),
  scaled to unit power in a lossless medium of index n_e = Re n_j:
      TE, horizontal   A = sqrt(3 s / (8 n_e)) / q_j        (+, +)
      TM, horizontal   A = sqrt(3 s / (8 n_e))              (+, -)
      TM, vertical     A = sqrt(3 s^3 / (4 n_e)) / q_j      (+, +)

s is integrated from 0 to s_max.  Above the largest real index of the stack
every wave is evanescent everywhere and what is left is near-field transfer to
absorbers, which for an emitter inside its own absorbing medium diverges and
is not radiation; s_max is the caller's choice (the largest real index)."""
import numpy as np

ORI = (('TE', 'h'), ('TM', 'h'), ('TM', 'v'))


def _q(eps, s):
    q = np.sqrt(eps - s**2 + 0j); return np.where(q.imag < 0, -q, q)


def _refl(eta, phi, order):
    """U-amplitude reflection of the stack of media `order` (list of indices, first = the medium the wave is in),
    referred to the first interface.  eta[i], phi[i] (phase thickness; unused for the first and last) arrays over s."""
    r = (eta[order[-2]] - eta[order[-1]]) / (eta[order[-2]] + eta[order[-1]])
    for k in range(len(order) - 2, 0, -1):
        i, j = order[k - 1], order[k]
        r12 = (eta[i] - eta[j]) / (eta[i] + eta[j]); e = np.exp(2j * phi[j])
        r = (r12 + r * e) / (1 + r12 * r * e)
    return r


def _carry(U, V, eta, phi):
    """(U, V) carried across a layer of phase thickness phi in the direction of the phase"""
    c, sn = np.cos(phi), np.sin(phi)
    return U * c + 1j * V * sn / eta, V * c + 1j * eta * U * sn


def dipole(wl, n, d, j, z, s, n_src=None):
    """wl: one wavelength (nm).  n: complex indices of media 0..M.  d: thicknesses of 1..M-1 (nm).
    j: emitter layer (1..M-1), z: depth in it from the side of medium 0 (nm).  s: array of in-plane indices.
    -> dict[(pol, ori)] = dict(K, F_first, F_last, A=[per layer 1..M-1], rr=(r_minus, r_plus))  (all arrays over s)."""
    n = np.asarray(n, complex); M = len(n) - 1; k0 = 2 * np.pi / wl; s = np.asarray(s, float)
    eps = n**2; q = [_q(e, s) for e in eps]
    dd = [0.0] + list(d) + [0.0]
    ne = float(np.real(n[j])) if n_src is None else float(n_src)
    out = {}
    for pol in ('TE', 'TM'):
        eta = [q[i] if pol == 'TE' else q[i] / eps[i] for i in range(M + 1)]
        phi = [k0 * q[i] * dd[i] for i in range(M + 1)]
        # reflection seen from the emitter layer at its two faces, then referred to the dipole plane
        r_dn = _refl(eta, phi, list(range(j, -1, -1))) * np.exp(2j * k0 * q[j] * z)
        r_up = _refl(eta, phi, list(range(j, M + 1))) * np.exp(2j * k0 * q[j] * (dd[j] - z))
        for ori in ('h', 'v'):
            if pol == 'TE' and ori == 'v': continue
            if pol == 'TE':   Au = np.sqrt(3 * s / (8 * ne)) / q[j]; Ad = Au
            elif ori == 'h':  Au = np.sqrt(3 * s / (8 * ne)) + 0j; Ad = -Au
            else:             Au = np.sqrt(3 * s**3 / (4 * ne)) / q[j]; Ad = Au
            den = 1 - r_up * r_dn
            a_p = (Au + r_dn * Ad) / den; b_p = r_up * a_p            # just above the dipole plane: up, down
            b_m = (Ad + r_up * Au) / den; a_m = r_dn * b_m            # just below
            # upward side: carry (U, V) to each interface up to medium M
            U, V = a_p + b_p, eta[j] * (a_p - b_p); flow = [np.real(np.conj(U) * V)]
            U, V = _carry(U, V, eta[j], k0 * q[j] * (dd[j] - z)); flow.append(np.real(np.conj(U) * V))
            for i in range(j + 1, M):
                U, V = _carry(U, V, eta[i], phi[i]); flow.append(np.real(np.conj(U) * V))
            up = flow                                                 # [at plane, top of j, top of j+1, ..., into M]
            # downward side: power flow toward medium 0 is -Re(U* V)
            U, V = a_m + b_m, eta[j] * (a_m - b_m); flow = [-np.real(np.conj(U) * V)]
            U, V = _carry(U, V, eta[j], -k0 * q[j] * z); flow.append(-np.real(np.conj(U) * V))
            for i in range(j - 1, 0, -1):
                U, V = _carry(U, V, eta[i], -phi[i]); flow.append(-np.real(np.conj(U) * V))
            dn = flow                                                 # [at plane, bottom of j, bottom of j-1, ..., into 0]
            A = [None] * (M + 1)
            A[j] = (up[0] - up[1]) + (dn[0] - dn[1])
            for k, i in enumerate(range(j + 1, M)): A[i] = up[k + 1] - up[k + 2]
            for k, i in enumerate(range(j - 1, 0, -1)): A[i] = dn[k + 1] - dn[k + 2]
            out[(pol, ori)] = dict(K=up[0] + dn[0], F_first=dn[-1], F_last=up[-1], A=A[1:M], rr=(r_dn, r_up))
    return out


def plane_wave(wl, n, d, s, side='first'):
    """A plane wave arriving from the first (or last) medium at in-plane index s: reflectance and the
    absorptance of every layer, per polarization.  -> dict[pol] = dict(R, T, A=[layers 1..M-1])."""
    n = np.asarray(n, complex)
    if side == 'last': n = n[::-1]; d = list(d)[::-1]
    M = len(n) - 1; k0 = 2 * np.pi / wl; s = np.asarray(s, float); eps = n**2; q = [_q(e, s) for e in eps]
    dd = [0.0] + list(d) + [0.0]; out = {}
    for pol in ('TE', 'TM'):
        eta = [q[i] if pol == 'TE' else q[i] / eps[i] for i in range(M + 1)]
        phi = [k0 * q[i] * dd[i] for i in range(M + 1)]
        r = _refl(eta, phi, list(range(0, M + 1)))
        U, V = 1 + r, eta[0] * (1 - r); f0 = np.real(np.conj(U) * V); inc = np.real(eta[0]); flow = [f0]
        for i in range(1, M):
            U, V = _carry(U, V, eta[i], phi[i]); flow.append(np.real(np.conj(U) * V))
        with np.errstate(divide='ignore', invalid='ignore'):
            A = [(flow[i - 1] - flow[i]) / inc for i in range(1, M)]; T = flow[-1] / inc
        A = A[::-1] if side == 'last' else A
        out[pol] = dict(R=np.abs(r)**2, T=T, A=A)
    return out


def isotropic(res, horiz=2 / 3):
    """orientation average: horiz = share of the horizontal dipoles (2/3 = isotropic).  -> same keys, per polarization"""
    v = 1 - horiz; out = {}
    for pol in ('TE', 'TM'):
        parts = [(horiz, res[(pol, 'h')])] + ([(v, res[(pol, 'v')])] if pol == 'TM' else [])
        out[pol] = dict(K=sum(w * p['K'] for w, p in parts), F_first=sum(w * p['F_first'] for w, p in parts), F_last=sum(w * p['F_last'] for w, p in parts),
                        A=[sum(w * p['A'][i] for w, p in parts) for i in range(len(parts[0][1]['A']))])
    return out


def s_grid(s_max, n=1500, crit=()):
    """in-plane indices on [0, s_max]: uniform in angle-like variable, refined toward the listed critical indices,
    never exactly on one (the source terms are singular at s = n_e)."""
    t = (np.arange(n) + 0.5) / n; s = s_max * np.sin(0.5 * np.pi * t)                 # denser toward s_max
    extra = [c + sg * c * 10.0**(-e) for c in crit if 0 < c <= s_max * (1 + 1e-12) for e in np.arange(0.8, 8.01, 0.1) for sg in (-1, 1)]
    s = np.unique(np.r_[s, [x for x in extra if 0 < x < s_max]])
    return s


def integrate(y, s):
    return np.trapezoid(y, s, axis=-1)

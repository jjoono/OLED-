"""Angle-resolved electroluminescence of the green PeLED of Nature 611, 688:

  glass (incoherent) / FTO 350 / BufHIL 75 / perovskite 270 / ZADN 45 / LiF 1.2 / Al 100

Far field by reciprocity: the power a dipole at depth z radiates into the
glass at angle theta_g with polarization s or p is proportional to |E(z)|^2 of
a unit plane wave coming in from the glass at that angle and polarization
(E_y for s; E_x and E_z for p; isotropic emitters take the mean of the three
components).  Emitters are spread over the perovskite thickness with a chosen
profile and over the internal emission spectrum; the perovskite's own k
(re-absorption) is in the fields.  Then glass -> air: Fresnel transmission,
the solid-angle change cos(theta_a) / (n_g^2 cos(theta_g)), and optionally the
incoherent multiple reflections between the glass/air face and the stack.

Output is radiant intensity I(theta_a) per unit emitting area, to compare with
Extended Data Fig. 7g, where a Lambertian source is cos(theta)."""
import os as _os
import sys, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
HC = 1239.841984


def n_al(wl):
    """Aluminium, Lorentz-Drude fit of Rakic et al., Appl. Opt. 37, 5271 (1998)."""
    E = HC / np.asarray(wl, float); wp = 14.98
    f = [0.523, 0.227, 0.050, 0.166, 0.030]; g = [0.047, 0.333, 0.312, 1.351, 3.382]; w0 = [0.0, 0.162, 1.544, 1.808, 3.473]
    eps = 1.0 - f[0] * wp**2 / (E * (E + 1j * g[0]))
    for fj, gj, wj in zip(f[1:], g[1:], w0[1:]):
        eps = eps + fj * wp**2 / (wj**2 - E**2 - 1j * E * gj)
    n = np.sqrt(eps); return np.where(n.imag < 0, -n, n)


def fields(wl, N, ds, th0, zj, j):
    """Plane wave, unit |E|, from medium 0 at angle th0 (rad, real medium 0).
    N: complex index arrays for media 0..M (0 and M semi-infinite); ds: thicknesses of 1..M-1.
    Returns r_s, r_p and the field intensities in layer j at depths zj (nm from its side facing medium 0):
    Ey2 [n_wl, n_z] (s) and Ex2, Ez2 (p)."""
    wl = np.asarray(wl, float); k0 = 2 * np.pi / wl
    eps = [np.asarray(n, complex)**2 * np.ones(len(wl)) for n in N]
    s0 = np.sqrt(eps[0]) * np.sin(th0)
    q = []
    for e in eps:
        c = np.sqrt(e - s0**2); q.append(np.where(c.imag < 0, -c, c))
    M = len(N) - 1; d = [0.0] + list(ds) + [0.0]
    out = {}
    for pol in ('s', 'p'):
        eta = [q[i] if pol == 's' else q[i] / eps[i] for i in range(M + 1)]
        A = [None] * (M + 1); Bc = [None] * (M + 1)
        A[M] = np.ones(len(wl), complex); Bc[M] = np.zeros(len(wl), complex)
        for i in range(M - 1, -1, -1):
            u = A[i + 1] + Bc[i + 1]; v = eta[i + 1] * (A[i + 1] - Bc[i + 1])
            ph = np.exp(1j * k0 * q[i] * d[i])
            A[i] = 0.5 * (u + v / eta[i]) / ph; Bc[i] = 0.5 * (u - v / eta[i]) * ph
        norm = A[0]; r = Bc[0] / norm
        a = (A[j] / norm)[:, None]; b = (Bc[j] / norm)[:, None]
        ef = np.exp(1j * (k0 * q[j])[:, None] * zj[None, :]); eb = 1.0 / ef
        if pol == 's':
            out['rs'] = r; out['Ey2'] = np.abs(a * ef + b * eb)**2
        else:                                   # amplitudes are H_y; incident H = n_0 for unit |E|
            n0 = np.sqrt(eps[0])[:, None]
            out['rp'] = r
            out['Ex2'] = np.abs((q[j] / eps[j])[:, None] * n0 * (a * ef - b * eb))**2
            out['Ez2'] = np.abs((s0 / eps[j])[:, None] * n0 * (a * ef + b * eb))**2
    return out


def angular(wl, spec, stack, ds, j, theta_air, profile=None, nz=40, n_glass=1.5, recycle_sub=True, horiz_frac=2 / 3):
    """stack: complex index arrays glass .. Al (glass first, last = Al semi-infinite or Al + air);
    j: index of the emitting layer; spec: internal emission spectrum on wl (photon weights);
    profile: function z/d in [0, 1] (0 = HIL side) -> weight; horiz_frac: 2/3 = isotropic dipoles.
    Returns I(theta_air) [n_theta] spectrally integrated, and the external spectrum per angle [n_theta, n_wl]."""
    dj = ds[j - 1]; z = (np.arange(nz) + 0.5) / nz * dj
    wz = np.ones(nz) if profile is None else np.array([profile(t) for t in z / dj]); wz = wz / wz.sum()
    ext = []
    for ta in np.deg2rad(theta_air):
        tg = np.arcsin(np.sin(ta) / n_glass)
        F = fields(wl, stack, ds, tg, z, j)
        # isotropic: 1/3 each of x, y, z; general: horizontal share split between x and y
        h, v = horiz_frac / 2, 1 - horiz_frac
        Is = h * (F['Ey2'] @ wz); Ip = h * (F['Ex2'] @ wz) + v * (F['Ez2'] @ wz)
        cg, ca = np.cos(tg), np.cos(ta)
        Rs = ((n_glass * cg - ca) / (n_glass * cg + ca))**2; Rp = ((cg - n_glass * ca) / (cg + n_glass * ca))**2
        geo = ca / (n_glass**2 * cg)
        ms = 1 / (1 - Rs * np.abs(F['rs'])**2) if recycle_sub else 1.0
        mp = 1 / (1 - Rp * np.abs(F['rp'])**2) if recycle_sub else 1.0
        ext.append(n_glass * geo * (Is * (1 - Rs) * ms + Ip * (1 - Rp) * mp) * spec)
    ext = np.array(ext)
    I = np.trapezoid(ext, wl, axis=1)
    return I, ext

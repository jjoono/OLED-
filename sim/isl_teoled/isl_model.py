"""Top-emitting OLED with an internal scattering layer (ISL) in front of the bottom Ag reflector,
compared with the planar TEOLED (coherent Ag reflector), at 550 nm.

  ISL device (bottom -> top):
     Ag (opaque) | ISL: matrix n_m, scattering optical thickness S, HG g, albedo  (incoherent, thick)
     | IZO 50 | ETL Bebq2 d_E | EML TCTA:B3PyMPM 25 (dipole at centre) | TAPC d_H | HATCN 10
     | Ag 15 (measured on HATCN) | TFE (Al2O3 30 / ParyC 500 / Al2O3 30 / ParyC 500 / Al2O3 30, optional)
     | lens body n_L (incoherent) | hexagonal close-packed hemisphere MLA (n_L) | air
  planar reference: Ag 100 | Bebq2 d_E | EML | TAPC d_H | HATCN 10 | Ag 15 | TFE | lens n_L | MLA | air

Two incoherent compartments, 1-degree angular bins:
  U = lens body (top boundary: MLA, B_T / B_R from the ray tracer; bottom: the coherent stack)
  D = ISL (bottom: slab + Ag mirror combined by the adding method; top: the coherent stack)
The coherent stack (IZO ... TFE) couples D and U at fixed n_eff = n sin(theta).
Dipole emission into D and into U from the CPS model (cps2.sub_angular on the stack and on its mirror).

   a = flux arriving at the MLA (up, in U)     b = flux leaving the stack downwards into D
   c = Reff b                                  (up, arriving at the stack from D)
   b = sD + T_UD BR a + Rd c,   a = sU + Ru BR a + T_DU c,   escape = BT . a
"""
import sys, numpy as np, scipy.io as sio
sys.path[:0] = ['/home/user/OLED-/sim/fig3c', '/home/user/OLED-/sim/mla', '/home/user/OLED-/sim/scatter']
import cps2

LAM = 550.0
m = sio.loadmat('/home/user/OLED-/sim/table1_author/nk_JH_total.mat', squeeze_me=True,
                struct_as_record=False)['material']
nk = lambda f: complex(np.asarray(getattr(m, f))[150])
AG = nk('l_Ag_McPeak'); AG15 = nk('l_Ag_15_onHATCN'); IZO = nk('l_IZO'); HAT = nk('l_HATCN')
BEB = nk('l_bebq2_JH'); AL2O3 = nk('l_Al2O3'); PARY = nk('l_ParyC')
TAPC = (nk('l_TAPC_o_JO').real, nk('l_TAPC_e_JO').real); EML = (nk('l_TCTA_B3_o_JO').real, nk('l_TCTA_B3_e_JO').real)
TFE = [(AL2O3, 30.0), (PARY, 500.0), (AL2O3, 30.0), (PARY, 500.0), (AL2O3, 30.0)]
TH = np.arange(90) + 0.5
EDGES = np.arange(91.0)


def iso(n, d):
    return (n, n, d)


# ---------- coherent stack between two semi-infinite incoherent media ----------
def _rt(layers, n_in, n_out, neff, pol):
    """power R, T for a plane wave from n_in; layers [(n (n+ik), d)], characteristic matrix in n-ik."""
    cj = lambda n: np.conj(complex(n))
    def c(n): return np.sqrt(1 - (neff / n) ** 2 + 0j)
    def eta(n):
        cc = c(n); cc = np.where(np.imag(cc) > 0, np.conj(cc), cc)
        return n * cc if pol == 's' else n / cc
    M = np.broadcast_to(np.eye(2, dtype=complex), neff.shape + (2, 2)).copy()
    for n, d in layers:
        n = cj(n); cc = c(n); cc = np.where(np.imag(cc) > 0, np.conj(cc), cc)
        dl = 2 * np.pi * n * cc * d / LAM; e = eta(n)
        L = np.empty_like(M)
        L[..., 0, 0] = np.cos(dl); L[..., 0, 1] = 1j * np.sin(dl) / e
        L[..., 1, 0] = 1j * e * np.sin(dl); L[..., 1, 1] = np.cos(dl)
        M = M @ L
    eo, ei = eta(cj(n_out)), eta(cj(n_in))
    B = M[..., 0, 0] + M[..., 0, 1] * eo; C = M[..., 1, 0] + M[..., 1, 1] * eo
    r = (ei * B - C) / (ei * B + C); t = 2 * ei / (ei * B + C)
    R = np.abs(r) ** 2
    T = np.real(eo) / np.real(ei) * np.abs(t) ** 2
    T = np.where(neff < np.real(n_out), T, 0.0)
    return R, T


def rt_unpol(layers, n_in, n_out, th_deg):
    neff = np.real(n_in) * np.sin(np.radians(th_deg))
    Rs, Ts = _rt(layers, n_in, n_out, neff, 's'); Rp, Tp = _rt(layers, n_in, n_out, neff, 'p')
    return 0.5 * (Rs + Rp), 0.5 * (Ts + Tp)


def transfer(layers, n_in, n_out):
    """specular R (diag, in n_in bins) and transmission matrix T[out_bin, in_bin] into n_out."""
    R, T = rt_unpol(layers, n_in, n_out, TH)
    Tm = np.zeros((90, 90))
    s = np.real(n_in) * np.sin(np.radians(TH)) / np.real(n_out)
    ok = s < 1
    j = np.floor(np.degrees(np.arcsin(s[ok]))).astype(int).clip(0, 89)
    Tm[j, np.where(ok)[0]] = T[ok]
    return R, Tm


# ---------- dipole emission into the two outer media ----------
def emission(below, above, n_below, n_above, h=0.76):
    """fractions of dissipated power per 1-deg bin entering the medium below and the medium above.
    below/above: lists of (n_o, n_e, d) outward from the EML."""
    Sd = cps2.Stack(LAM, EML, 25.0, 12.5, above=above, below=below, n_sub=n_below, n_top=n_above, h=h)
    Su = cps2.Stack(LAM, EML, 25.0, 12.5, above=below, below=above, n_sub=n_above, n_top=n_below, h=h)
    Pd = cps2.solve(Sd)['P_tot']; Pu = cps2.solve(Su)['P_tot']
    sD = np.real(cps2.sub_angular(Sd, TH)) * (180 / np.pi) * np.radians(1.0) / Pd if n_below is not None else 0
    sU = np.real(cps2.sub_angular(Su, TH)) * (180 / np.pi) * np.radians(1.0) / Pu
    return np.nan_to_num(sD), np.nan_to_num(sU), Pd, Pu


def mirror_R(n_m):
    R, _ = rt_unpol([], n_m, AG, TH)
    return R


def solve_isl(sD, sU, Reff, Rd, T_DU, Ru, T_UD, BT, BR):
    I = np.eye(90)
    A11 = I - Ru[:, None] * BR;           A12 = -T_DU @ Reff
    A21 = -T_UD @ BR;                     A22 = I - Rd[:, None] * Reff
    X = np.linalg.solve(np.block([[A11, A12], [A21, A22]]), np.concatenate([sU, sD]))
    a = X[:90]
    return float(BT @ a)


def solve_planar(sU, Ru, BT, BR):
    a = np.linalg.solve(np.eye(90) - Ru[:, None] * BR, sU)
    return float(BT @ a)

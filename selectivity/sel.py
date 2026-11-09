"""Directional selectivity of multi-pass OLED extraction: device | randomizer | angular filter | air.

Light delivered to the substrate passes a randomizing layer (an index-matched volume scatterer,
or an idealized partial randomizer), meets an angular filter on the flat exit face, and the
rejected light is returned to the OLED, reflected with R_LED(theta, lambda) and sent back for
another attempt.  This is the trans-scale multi-pass series of the round-trip-loss paper
(Kim, Lee, Yoo, Gather; Suppl. Note 1) with the outcoupling structure replaced by
"randomizer + filter", solved per wavelength in matrix form:

    u  = (I - Rs Df)^-1 Ts v           upward flux reaching the filter (all internal bounces)
    BT = diag(Tf) U,   BR = Rs + Ts Df U
    x  = (I - diag(R_LED) BR)^-1 P_sub  flux entering the randomizer, summed over passes
    air (by substrate angle at the exit face) = BT x

Ts, Rs: transmission / reflection matrices of the randomizer (column = input bin, symmetric slab),
Tf, Rf: filter transmittance into air / reflectance back into the substrate.
Device optics: CPS dipole model + TMM (cps2), dispersive optical constants of the group library.
"""
import os, sys, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'collimation', 'cps'))
import cps2

D = dict(np.load(os.path.join(HERE, 'data', 'inputs.npz')))
LAM1 = D['lam']                                   # 400..700, 1 nm
LAM = np.arange(400, 701, 2).astype(float)        # working grid
IL = (LAM - 400).astype(int)

NB = 180                                          # 0.5 deg bins in the substrate
EDG = np.linspace(0, 90, NB + 1)
TH = 0.5 * (EDG[:-1] + EDG[1:])
S2 = np.sin(np.radians(EDG)) ** 2
WL = np.diff(S2); WL = WL / WL.sum()              # Lambertian power per bin
BETAS = (5.0, 10.0, 15.0, 20.0, 30.0, 45.0, 60.0, 90.0)   # air-cone half angles; 5 deg -> on-axis


def nk(name, j=None):
    v = D['n_' + name]
    return v[IL] if j is None else v[j]


# ---------------------------------------------------------------- spectra
def gauss(lam0, fwhm):
    s = np.exp(-4 * np.log(2) * ((LAM - lam0) / fwhm) ** 2)
    return s / s.sum()


def real_spec(name):
    s = D['s_' + name][IL] * LAM                  # photon-number weighting (library spectra are energy)
    return s / s.sum()


# ---------------------------------------------------------------- devices
def stack(name, j):
    """Bottom-emitting stacks.  Returns (layers above EML [to rear], below [to substrate], n_sub)."""
    org = lambda m, d: (nk(m + '_o', j).real, nk(m + '_e', j).real, d)
    met = lambda m, d: (nk(m, j), nk(m, j), d)
    if name in ('Al', 'Ag'):                       # platform of the MLA study: ETL 61.3 / HTL 97.5, ITO 150, glass
        above = [org('B3', 61.3), met(name, 100.0)]
        below = [org('TCTA', 10.0), org('TAPC', 97.5), met('ITO', 150.0)]
        return above, below, 1.51
    if name in ('LL15', 'LL18'):                   # low-loss stack of Ext. Data Table 1: Ag/B3 200/EML/TAPC 180/ITO 50
        above = [org('B3', 200.0), met('Ag', 100.0)]
        below = [org('TAPC', 180.0), met('ITO', 50.0)]
        return above, below, (1.51 if name == 'LL15' else 1.80)
    if name == 'DBRrear':                          # ITO 50 cathode + chirped 10-pair ZnS/LiF rear DBR (Fig. 2e)
        above = [org('B3', 61.3), met('ITO', 50.0)]
        for k in range(10):
            f = k / 9
            above += [met('ZnS', 57 + 23 * f), met('LiF', 88 + 36 * f)]
        below = [org('TCTA', 10.0), org('TAPC', 97.5), met('ITO', 150.0)]
        return above, below, 1.51
    raise KeyError(name)


STACKS = {'Al': 'Al cathode, ITO 150 nm, glass', 'Ag': 'Ag cathode, ITO 150 nm, glass',
          'DBRrear': 'ITO cathode + ZnS/LiF rear DBR, glass', 'LL15': 'low-loss Ag/ITO 50 nm, n_sub 1.51',
          'LL18': 'low-loss Ag/ITO 50 nm, n_sub 1.80'}


def device(name, theta_h=0.865):
    """eta_sub(lam), P_sub(lam, bin) (normalized), R_LED(lam, bin) on the working grid."""
    es, P, R = [], [], []
    for j, lam in zip(IL, LAM):
        above, below, ns = stack(name, j)
        S = cps2.Stack(float(lam), (nk('EML_o', j).real, nk('EML_e', j).real), 25.0, 12.5,
                       above=above, below=below, n_sub=ns, h=theta_h)
        r = cps2.solve_pol(S, npts=4000)
        es.append(r['air'] + r['sub'])
        th_f = np.linspace(0, 90, 1801)[:-1] + 0.025
        pf = np.clip(np.real(cps2.sub_angular(S, th_f)), 0, None)
        p = pf.reshape(NB, 10).sum(1); P.append(p / p.sum())
        R.append(np.clip(np.real(cps2.stack_reflectance(S, TH)), 0, 1))
    return dict(eta_sub=np.array(es), P=np.array(P), R=np.array(R), n_sub=stack(name, 0)[2])


# ---------------------------------------------------------------- filters (TMM)
def tmm_RT(layers, lam, th_deg, n_in, n_out=1.0):
    """Unpolarized R, T of a multilayer [(n complex, d nm), ...] from n_in at angles th_deg into n_out."""
    kx = n_in * np.sin(np.radians(th_deg)).astype(complex)
    R = np.zeros(th_deg.size); T = np.zeros(th_deg.size)
    for pol in ('s', 'p'):
        def eta(n):
            c = np.sqrt(n ** 2 - kx ** 2 + 0j)
            c = np.where(c.imag < 0, -c, c)
            return c if pol == 's' else n ** 2 / c, c
        M = np.tile(np.eye(2, dtype=complex), (th_deg.size, 1, 1))
        for n, d in layers:
            e, c = eta(complex(n)); dl = 2 * np.pi / lam * c * d
            m = np.empty((th_deg.size, 2, 2), complex)
            m[:, 0, 0] = np.cos(dl); m[:, 0, 1] = 1j * np.sin(dl) / e
            m[:, 1, 0] = 1j * e * np.sin(dl); m[:, 1, 1] = np.cos(dl)
            M = M @ m
        e0, _ = eta(complex(n_in)); es, cs = eta(complex(n_out))
        B = M[:, 0, 0] + M[:, 0, 1] * es; C = M[:, 1, 0] + M[:, 1, 1] * es
        r = (e0 * B - C) / (e0 * B + C); t = 2 * e0 / (e0 * B + C)
        Tp = np.real(es) / np.real(e0) * np.abs(t) ** 2
        Tp[np.abs(cs.real) < 1e-9] = 0.0
        R += 0.5 * np.abs(r) ** 2; T += 0.5 * np.nan_to_num(Tp)
    return np.clip(R, 0, 1), np.clip(T, 0, 1)


DBR_MATS = {'TiO2/SiO2': ('TiO2', 'SiO2'), 'SiNx/SiO2': ('SiNx', 'SiO2'), 'ZnS/LiF': ('ZnS', 'LiF')}


def edge_dbr(mat, N, lam_b):
    """(L/2 H L/2)^N quarter-wave at lam_b: stop band on the red side of the emitter at normal
    incidence that blue-shifts onto it at oblique incidence (passes a forward cone)."""
    H, L = DBR_MATS[mat]
    jb = int(np.clip(round(lam_b) - 400, 0, 300))
    nH, nL = nk(H, jb).real, nk(L, jb).real
    dH, dL = lam_b / (4 * nH), lam_b / (4 * nL)
    return [(L, dL / 2), (H, dH), (L, dL / 2)] * N


def filter_RT(spec, n_sub):
    """spec: 'fresnel' | ('ideal', beta) | ('edge', mat, N, lam_b) | ('layers', [(mat, d), ...]).
    Returns Rf, Tf arrays (lam, bin)."""
    Rf = np.zeros((LAM.size, NB)); Tf = np.zeros((LAM.size, NB))
    if isinstance(spec, tuple) and spec[0] == 'ideal':
        tc = np.degrees(np.arcsin(np.sin(np.radians(spec[1])) / n_sub))
        f = np.clip((np.sin(np.radians(np.minimum(EDG[1:], tc))) ** 2 - S2[:-1]) / np.diff(S2), 0, 1)
        Tf[:] = f; Rf[:] = 1 - f
        return Rf, Tf
    if spec == 'fresnel':
        layers = []
    elif spec[0] == 'edge':
        layers = edge_dbr(*spec[1:])
    else:
        layers = spec[1]
    for i, (j, lam) in enumerate(zip(IL, LAM)):
        lay = [(nk(m, j), d) for m, d in layers]
        Rf[i], Tf[i] = tmm_RT(lay, lam, TH, n_sub)
    return Rf, Tf


# ---------------------------------------------------------------- randomizers
def mixer(m):
    """Idealized partial randomizer: a fraction m of each traversal is redistributed Lambertian."""
    Ts = (1 - m) * np.eye(NB) + m * np.outer(WL, np.ones(NB))
    return Ts, np.zeros((NB, NB))


# ---------------------------------------------------------------- multi-pass solve
def solve(dev, Rf, Tf, Ts, Rs):
    """Air output per wavelength, binned by substrate angle at the exit face: (lam, bin)."""
    out = np.zeros((LAM.size, NB)); I = np.eye(NB)
    for i in range(LAM.size):
        Df = Rf[i]
        U = np.linalg.solve(I - Rs * Df[None, :], Ts)
        BR = Rs + (Ts * Df[None, :]) @ U
        x = np.linalg.solve(I - dev['R'][i][:, None] * BR, dev['P'][i])
        out[i] = Tf[i] * (U @ x)
    return out


def cone_frac(n_sub):
    """F[k, bin]: fraction of each substrate bin that exits inside air cone BETAS[k] (flat exit)."""
    F = np.zeros((len(BETAS), NB))
    for k, b in enumerate(BETAS):
        tc = np.degrees(np.arcsin(min(1.0, np.sin(np.radians(b)) / n_sub)))
        F[k] = np.clip((np.sin(np.radians(np.minimum(EDG[1:], tc))) ** 2 - S2[:-1]) / np.diff(S2), 0, 1)
    return F


def cones(air, n_sub):
    """(lam, len(BETAS)) power per unit substrate-delivered power inside each air cone."""
    return air @ cone_frac(n_sub).T


def metrics(eta_sub, cone_l, spec):
    """EQE into each cone (photon basis, PLQY = 1) and the forward 'Lambertian-equivalent' EQE."""
    e = spec @ (eta_sub[:, None] * cone_l)
    out = {f'EQE_{b:g}': float(v) for b, v in zip(BETAS, e)}
    out['EQE_total'] = out['EQE_90']
    out['EQE_fwd'] = out['EQE_5'] / np.sin(np.radians(5.0)) ** 2   # on-axis luminance in Lambertian-EQE units
    return out

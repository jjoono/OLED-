"""Spectral TEOLED model, Ir(mphmq)2tmd in Bebq2, top-emitting through thin Ag on HATCN, Parylene-C TFE, n=1.65 MLA.
  Ag 100 | [Al a] | Liq 1 | Bebq2 d_E | Bebq2:Ir(mphmq)2tmd 25 (dipole at z0) | TAPC d_H | HATCN d_HAT | Ag d_Ag
  | Al2O3 30 | ParyC 500 | Al2O3 30 | ParyC 500 | Al2O3 30 | lens body n 1.65 | hemisphere MLA | air
eta_out(lambda) = full recycling series (CPS emission into the lens body, stack reflectance seen from the lens, MLA BSDF);
EQE = eta_rad * sum_lambda S(lambda) eta_out(lambda)   (Purcell change of eta_rad neglected)."""
import sys, numpy as np, scipy.io as sio
sys.path[:0] = ['/home/user/OLED-/sim/fig3c']
import cps2
M = sio.loadmat('/home/user/OLED-/sim/table1_author/nk_JH_total.mat', squeeze_me=True, struct_as_record=False)
mat, spec = M['material'], M['spectrum']
LAMS = np.linspace(400, 800, 401)
def nk(f): return np.asarray(getattr(mat, f), dtype=complex)
AG, AG15, AL, HAT = nk('l_Ag_McPeak'), nk('l_Ag_15_onHATCN'), nk('l_Al_JO'), nk('l_HATCN')
AG12T5 = nk('l_Ag_12nm_test5')
BEB = np.real(nk('l_bebq2_JH')) + 1j * np.imag(nk('l_bebq2_JH'))
AL2O3, PARY = nk('l_Al2O3'), nk('l_ParyC')
ZNS = nk('l_ZnS')
NPB = np.real(nk('l_NPB'))
IZOk = nk('l_IZO'); LIF = 1.408 + 0.0002j
_d = np.load('/home/user/OLED-/sim/dbr_design/dbr_maps.npz')
DBR = {'ch': list(_d['d_ch']), 'qw': list(_d['d_qw']), 'chr': list(_d['d_ch'])[::-1]}
TAPCo, TAPCe = np.real(nk('l_TAPC_o_JO')), np.real(nk('l_TAPC_e_JO'))
LIQ = 1.70; NL = 1.65; ETA_RAD = 0.96; THETA = 0.82
S = np.clip(np.real(np.asarray(spec.l_I_Irmphmq2tmd_measure_JH, dtype=float)), 0, None)
TH = np.arange(90) + 0.5
bs = np.load('/home/user/OLED-/sim/teoled_opt/mla_hemi_1.65.npz') if __import__('os').path.exists('/home/user/OLED-/sim/teoled_opt/mla_hemi_1.65.npz') else None


def idx(lam): return int(round(lam - 400))


def _rt(layers, n_in, n_out, neff, pol, lam):
    def c(n):
        cc = np.sqrt(1 - (neff / n) ** 2 + 0j); return np.where(np.imag(cc) > 0, np.conj(cc), cc)
    def eta(n): return n * c(n) if pol == 's' else n / c(n)
    Mx = np.broadcast_to(np.eye(2, dtype=complex), neff.shape + (2, 2)).copy()
    for n, d in layers:
        n = np.conj(complex(n)); cc = c(n); dl = 2 * np.pi * n * cc * d / lam; e = eta(n)
        L = np.empty_like(Mx)
        L[..., 0, 0] = np.cos(dl); L[..., 0, 1] = 1j * np.sin(dl) / e; L[..., 1, 0] = 1j * e * np.sin(dl); L[..., 1, 1] = np.cos(dl)
        Mx = Mx @ L
    eo, ei = eta(np.conj(complex(n_out))), eta(np.conj(complex(n_in)))
    B = Mx[..., 0, 0] + Mx[..., 0, 1] * eo; C = Mx[..., 1, 0] + Mx[..., 1, 1] * eo
    return np.abs((ei * B - C) / (ei * B + C)) ** 2


def R_from_lens(layers, lam, n_out=1.0):
    neff = NL * np.sin(np.radians(TH))
    return 0.5 * (_rt(layers, NL, n_out, neff, 's', lam) + _rt(layers, NL, n_out, neff, 'p', lam))


def stack(lam, dE, dH, dHAT, dAg, agtop='meas', al=0.0, z0=12.5, tfe=True, n_etl=None, n_eml=None, bottom='Ag', zns=0.0, htl='TAPC', mirror='Ag', izo=100.0, dbr='ch', n_back=1.0, topm='Ag', izo_nk=None):
    i = idx(lam)
    ne = BEB[i] if n_etl is None else n_etl
    mb = AG[i] if bottom == 'Ag' else AL[i]
    ag_t = (IZOk[i] if izo_nk is None else izo_nk[i]) if topm == 'IZO' else agtop[i] if isinstance(agtop, np.ndarray) else {'meas': AG15[i], 'bulk': AG[i], 'test5': AG12T5[i]}[agtop] if isinstance(agtop, str) else AG15[i] + agtop * (AG[i] - AG15[i])
    above = [(ne, ne, dE), (LIQ, LIQ, 1.0)] + ([(AL[i], AL[i], al)] if al > 0 else [])
    if mirror == 'Ag':
        above += [(mb, mb, 100.0)]
    else:                                   # IZO electrode on a ZnS/LiF DBR (ZnS first from the IZO side)
        above += [(IZOk[i], IZOk[i], izo)]
        for j, d in enumerate(DBR[dbr]):
            m_ = ZNS[i] if j % 2 == 0 else LIF
            above += [(m_, m_, d)]
        if mirror == 'DBR+Ag':
            above += [(AG[i], AG[i], 100.0)]
    hl = (TAPCo[i], TAPCe[i], dH) if htl == 'TAPC' else (NPB[i], NPB[i], dH)
    below = [hl, (HAT[i], HAT[i], dHAT), (ag_t, ag_t, dAg)]
    if zns > 0:
        below += [(ZNS[i], ZNS[i], zns)]
    if tfe == 'al2o3':
        below += [(AL2O3[i], AL2O3[i], 30.)]
    elif tfe:
        below += [(AL2O3[i], AL2O3[i], 30.), (PARY[i], PARY[i], 500.), (AL2O3[i], AL2O3[i], 30.), (PARY[i], PARY[i], 500.), (AL2O3[i], AL2O3[i], 30.)]
    return above, below, i


def eta_out(lam, dE, dH, dHAT, dAg, BT, BR, **kw):
    above, below, i = stack(lam, dE, dH, dHAT, dAg, **kw)
    ne_ = np.real(BEB[i]) if kw.get('n_eml') is None else kw['n_eml']
    e = (ne_, ne_)
    nb = kw.get('n_back', 1.0) if kw.get('mirror', 'Ag') == 'DBR' else 1.0
    St = cps2.Stack(lam, e, 25.0, kw.get('z0', 12.5), above=above, below=below, n_sub=NL, n_top=nb, h=THETA)
    P = np.real(cps2.solve(St)['P_tot'])
    sU = np.nan_to_num(np.real(cps2.sub_angular(St, TH))) / P
    lay = [(l[0], l[2]) for l in below[::-1]] + [(e[0], 25.0)] + [(l[0], l[2]) for l in above]
    Ru = R_from_lens(lay, lam, nb)
    a = np.linalg.solve(np.eye(90) - Ru[:, None] * BR, sU)
    return float(BT @ a), float(sU.sum())


def eqe(dE, dH, dHAT, dAg, BT, BR, lams=(580, 600, 615, 630, 650, 680), **kw):
    w = np.array([S[idx(l)] for l in lams]); w = w / w.sum()
    eo = np.array([eta_out(l, dE, dH, dHAT, dAg, BT, BR, **kw)[0] for l in lams])
    return ETA_RAD * float(w @ eo)

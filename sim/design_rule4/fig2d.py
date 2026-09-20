# -*- coding: utf-8 -*-
"""Fig. 2(d): angle- and wavelength-resolved round-trip loss of the candidate electrode
structures, for light arriving from the substrate.

1 - R is what the recycling process actually loses per round trip; for an opaque reflector it
is all absorption, for the dielectric mirror part of it is transmission, so the two are
reported separately.  Dispersive n,k from nk_JH_total.mat (Koenig ITO read from
nk_ITO_Konig2014.csv); the organic stack is the generic non-absorbing n = 1.8 slab used
throughout Fig. 2, so the panel isolates the electrodes."""
import numpy as np, scipy.io as sio, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
LAM = np.arange(400.0, 801.0)                      # the library grid, 1 nm steps
_m = sio.loadmat(os.path.join(ROOT, 'nk_JH_total.mat'), squeeze_me=True, struct_as_record=False)
_mat, _spec = _m['material'], _m['spectrum']
def lib(name): return np.atleast_1d(np.asarray(getattr(_mat, name))).astype(complex)

def konig():
    raw = open(os.path.join(HERE, 'nk_ITO_Konig2014.csv')).read().split('wl,k')
    n = np.array([[float(x) for x in l.split(',')] for l in raw[0].strip().splitlines()[1:] if l.strip()])
    k = np.array([[float(x) for x in l.split(',')] for l in raw[1].strip().splitlines()
                  if l.strip() and not l.startswith('wl')])
    return np.interp(LAM/1000, n[:, 0], n[:, 1]) + 1j*np.interp(LAM/1000, k[:, 0], k[:, 1])

AG, AL, IZO, ITO = lib('l_Ag_McPeak'), lib('l_Al_JO'), lib('l_IZO'), konig()
ZNS, LIF = lib('l_ZnS'), lib('l_LiF')
ORG = np.full_like(LAM, 1.8, dtype=complex)
AIR = np.ones_like(LAM, dtype=complex)
ORG_D = 420.0            # HTL 200 + EML 20 + ETL 200, as in Fig. 2(a)-(c)

def dbr(pairs=4.5, d_zns=70.0, d_lif=115.0):
    """ZnS/LiF stack, ZnS first (the high-index layer faces the organics).
    Defaults are the thicknesses deposited on the orange device."""
    out = []
    n_full = int(pairs)
    for _ in range(n_full):
        out += [(ZNS, d_zns), (LIF, d_lif)]
    if pairs - n_full >= 0.5:
        out += [(ZNS, d_zns)]
    return out

def dbr_qw(pairs=10, lam0=550.0):
    """Quarter-wave stack at lam0, using the measured indices there."""
    i = int(np.argmin(abs(LAM - lam0)))
    return dbr(pairs, lam0/(4*ZNS[i].real), lam0/(4*LIF[i].real))

def dbr_chirp(pairs=10, d_zns=56.0, d_lif=88.0, factor=1.42):
    """Linearly chirped ZnS/LiF stack: both thicknesses scale from 1 to `factor` across the
    stack, which widens the stopband instead of deepening it.  The defaults minimise the
    flux- and spectrum-weighted round-trip loss on glass with the green emitter."""
    s = np.linspace(1.0, factor, pairs)
    out = []
    for i in range(pairs):
        out += [(ZNS, d_zns*s[i]), (LIF, d_lif*s[i])]
    return out

QW = dbr_qw()
D_ZNS, D_LIF = QW[0][1], QW[1][1]
# thicknesses that minimise the flux- and spectrum-weighted loss for randomised light
# (grid search, 10 pairs): glass + green emitter, and n_sub = 1.8 + orange emitter
OPT15 = dbr(10, 69.0, 97.0)
OPT18 = dbr(10, 73.0, 110.0)

# ---- the structures, listed from the substrate inwards -----------------------
def stacks(d_ito=150.0):
    """Everything fixed except the reflector: the transparent electrode is 150 nm of ITO
    in all three.  The DBR replaces the metal as a purely optical mirror; the transparent
    cathode it would need in a real device is the variant below."""
    common = [(ITO, d_ito), (ORG, ORG_D)]
    return {
      'Al':                dict(layers=common + [(AL, 100.0)], exit=AIR,
                                note='conventional Al cathode'),
      'Ag':                dict(layers=common + [(AG, 100.0)], exit=AIR,
                                note='low-loss Ag cathode'),
      'DBR (10 pairs)':    dict(layers=common + QW, exit=AIR,
                                note=f'ZnS {D_ZNS:.0f} nm / LiF {D_LIF:.0f} nm, quarter-wave at 550 nm'),
      'DBR, chirped':      dict(layers=common + dbr_chirp(), exit=AIR,
                                note='10 pairs, ZnS 56->80 nm / LiF 88->125 nm; the chirp widens the stopband'),
      'DBR, chirped 20':   dict(layers=common + dbr_chirp(20, 67.0, 71.0, 1.42), exit=AIR,
                                note='the same idea with 20 pairs'),
      'DBR, re-optimised': dict(layers=common + OPT15, exit=AIR,
                                note='ZnS 69 nm / LiF 97 nm, minimising the loss for randomised light on glass'),
      'DBR, re-opt. 1.8':  dict(layers=common + OPT18, exit=AIR,
                                note='ZnS 73 nm / LiF 110 nm, the same for a high-index substrate'),
      'DBR + IZO cathode': dict(layers=common + [(IZO, 50.0)] + QW, exit=AIR,
                                note='the quarter-wave stack with the transparent cathode a real device needs'),
    }

def stacks_v1():
    """The earlier set, kept for reference: the transparent electrode varied too."""
    return {
      'Al / ITO 150 nm':   dict(layers=[(ITO, 150.0), (ORG, ORG_D), (AL, 100.0)], exit=AIR, note=''),
      'Ag / ITO 150 nm':   dict(layers=[(ITO, 150.0), (ORG, ORG_D), (AG, 100.0)], exit=AIR, note=''),
      'Ag / IZO 50 nm':    dict(layers=[(IZO, 50.0), (ORG, ORG_D), (AG, 100.0)], exit=AIR, note=''),
      'Ag / Ag 10 nm':     dict(layers=[(AG, 10.0), (ORG, ORG_D), (AG, 100.0)], exit=AIR, note=''),
      'DBR / IZO 50 nm':   dict(layers=[(IZO, 50.0), (ORG, ORG_D), (IZO, 50.0)] + dbr(), exit=AIR, note=''),
    }

# ---- TMM --------------------------------------------------------------------
def RT(layers, n_inc, n_exit, th, lam):
    """R and T for a plane wave from n_inc through the finite layers into n_exit.
    th: incidence angles in the incident medium (rad).  Arrays are (len(lam), len(th))."""
    ns = [n_inc] + [l[0] for l in layers] + [n_exit]
    ds = [0.0] + [l[1] for l in layers] + [0.0]
    kx = n_inc[:, None]*np.sin(th)[None, :]
    kz = [np.sqrt(n[:, None]**2 - kx**2 + 0j) for n in ns]
    out = {}
    for pol in ('s', 'p'):
        eta = kz if pol == 's' else [ns[i][:, None]**2/kz[i] for i in range(len(ns))]
        M = np.zeros(kx.shape + (2, 2), dtype=complex)
        M[..., 0, 0] = M[..., 1, 1] = 1.0
        for i in range(len(ns)-1):
            r = (eta[i] - eta[i+1])/(eta[i] + eta[i+1])
            t = 2*eta[i]/(eta[i] + eta[i+1])
            I = np.empty(kx.shape + (2, 2), dtype=complex)
            I[..., 0, 0] = 1/t; I[..., 0, 1] = r/t; I[..., 1, 0] = r/t; I[..., 1, 1] = 1/t
            M = M @ I
            if ds[i+1] > 0:
                ph = 2*np.pi/lam[:, None]*kz[i+1]*ds[i+1]
                L = np.zeros(kx.shape + (2, 2), dtype=complex)
                L[..., 0, 0] = np.exp(-1j*ph); L[..., 1, 1] = np.exp(1j*ph)
                M = M @ L
        r = M[..., 1, 0]/M[..., 0, 0]
        t = 1.0/M[..., 0, 0]
        R = np.abs(r)**2
        T = np.abs(t)**2*np.real(eta[-1]/eta[0])
        out[pol] = (R, np.where(np.isfinite(T), T, 0.0))
    R = 0.5*(out['s'][0] + out['p'][0]); T = 0.5*(out['s'][1] + out['p'][1])
    return R, T

def maps(name, n_sub=1.5, nth=451, lam=None, table=None):
    lam = LAM if lam is None else lam
    sel = np.isin(LAM, lam)
    S = (table or stacks())[name]
    layers = [(m[sel], d) for m, d in S['layers']]
    th = np.linspace(0, np.pi/2, nth + 1); th = 0.5*(th[1:] + th[:-1])
    nsub = np.full(lam.shape, n_sub, dtype=complex)
    R, T = RT(layers, nsub, S['exit'][sel], th, lam)
    return th, lam, R, T

def weighted(name, n_sub=1.5, spectrum=None, lam=None, table=None):
    """Flux-weighted (cos.sin, the ergodic weight an MLA enforces) and spectrum-weighted
    round-trip loss: absorbed, transmitted, and their sum = A' = 1 - <R>."""
    th, lam_, R, T = maps(name, n_sub=n_sub, lam=lam, table=table)
    w = np.cos(th)*np.sin(th); w /= w.sum()
    Ra, Ta = (R*w).sum(1), (T*w).sum(1)
    if spectrum is None:
        sw = np.ones_like(lam_)
    else:
        sw = np.interp(lam_, LAM, np.atleast_1d(np.asarray(spectrum)).real)
    sw = np.clip(sw, 0, None); sw /= sw.sum()
    return dict(absorbed=float(((1-Ra-Ta)*sw).sum()), transmitted=float((Ta*sw).sum()),
                loss=float(((1-Ra)*sw).sum()))

GREEN = np.atleast_1d(np.asarray(getattr(_spec, 'l_I_Irppy2acac'))).real
ORANGE = np.atleast_1d(np.asarray(getattr(_spec, 'l_I_Irdmppyph2tmd'))).real

if __name__ == '__main__':
    vis = np.arange(430.0, 701.0)
    for n_sub, spec, tag in ((1.5, GREEN, 'n_sub 1.5, green emitter'),
                             (1.8, ORANGE, 'n_sub 1.8, orange emitter')):
        print(f'--- {tag} ---')
        print(f"{'structure':20s} {'absorbed':>9} {'transmitted':>12} {'A (=1-<R>)':>11}")
        for k in stacks():
            w = weighted(k, n_sub=n_sub, spectrum=spec, lam=vis)
            print(f"{k:20s} {w['absorbed']:9.3f} {w['transmitted']:12.3f} {w['loss']:11.3f}")
        print()

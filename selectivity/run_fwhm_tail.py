"""Emitter requirement map for >1000 cd/A: FWHM x red-tail fraction.
Spectrum = (1-f) G(530 nm, w) + f G(566 nm, 40 nm)  (vibronic shoulder ~0.16 eV below the 0-0 band).
Device: low-loss n 1.8 stack + scattering slab (g 0.8, tau* 1) + front filter: best quarter-wave edge DBR for
on-axis luminance (TiO2/SiO2 or ZnS/LiF, 8/12/16 pairs, edge scanned), ideal 10/20 deg cones, no filter.
References per spectrum: hemispherical MLA, best thin-Ag strong cavity (Ag cathode, ETL/HTL scan).
Outputs: on-axis cd/A, total EQE, du'v' (<=30, <=60 deg), CIE 1931 (x, y) on axis."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import json, itertools, numpy as np
from multiprocessing import Pool
import sel, optfilter as O, run_display as RD, common

FW = (10, 15, 20, 25, 30, 40); TAIL = (0.0, 0.05, 0.10, 0.15, 0.20, 0.30)
DV = RD.DV; SL = RD.SL; MLA = dict(np.load('data/ref_mla.npz')); CAV = np.load('data/ref_cavity.npz')
dev = dict(eta_sub=DV['LL18_eta_sub'], P=DV['LL18_P'], R=DV['LL18_R'], n_sub=1.8)
SLAB = tuple(SL['g0.8_t1.0_w1'])
CIE = np.load('data/cie.npz'); XB, YB, ZB = (CIE[k][sel.IL] for k in 'xyz')


def spectrum(w, f):
    g = lambda c, s: np.exp(-4 * np.log(2) * ((sel.LAM - c) / s) ** 2)
    a, b = g(530.0, w), g(566.0, 40.0)
    s = (1 - f) * a / a.sum() + f * b / b.sum()
    return s / s.sum()


def k_cdA(spec):
    return 683 * (spec * YB * 1239.84 / sel.LAM).sum() / np.pi


def xy(S):
    X, Y, Z = ((S * c / sel.LAM).sum() for c in (XB, YB, ZB))
    return X / (X + Y + Z), Y / (X + Y + Z)


def one(args):
    w, f = args
    spec = spectrum(w, f); name = f'W{w}_T{f}'
    common.SPECTRA[name] = spec                       # optfilter.Problem looks the spectrum up by name
    best = None
    for mat in ('TiO2/SiO2', 'ZnS/LiF'):
        P = O.Problem('LL18', name, 5.0, mat, key='fwd')
        for N in (8, 12, 16):
            for lb in np.arange(540.0, 760.0, 5.0):
                v = P.eqe(O.edge_seq(mat, N, lb))
                if best is None or v > best[0]:
                    best = (v, mat, N, float(lb))
    out = dict(fwhm=w, tail=f, dbr_design=f'{best[1]} {best[2]}p {best[3]:.0f} nm')
    for lab, filt in (('nofilter', 'fresnel'), ('dbr', ('edge',) + best[1:]), ('ideal20', ('ideal', 20.0)), ('ideal10', ('ideal', 10.0))):
        Rf, Tf = sel.filter_RT(filt, 1.8)
        air = sel.solve(dev, Rf, Tf, *SLAB)
        m = RD.display_metrics(RD.to_air_angle(air, 1.8), dev['eta_sub'], spec)
        I = RD.to_air_angle(air, 1.8)[:, 0] * spec * dev['eta_sub']
        out[lab] = dict(cdA=m['cdA0'], EQE=m['EQE'], duv30=m['duv30'], duv60=m['duv60'], half=m['half_lum_angle'], xy=xy(I))
    k = k_cdA(spec)
    ml = spec @ (dev['eta_sub'][:, None] * MLA['LL18'])
    out['MLA'] = dict(cdA=float(ml[0] / np.sin(np.radians(5)) ** 2 * k), EQE=float(ml[-1]))
    E = np.einsum('l,kl,klb->kb', spec, CAV['eta_sub'], CAV['cones'])
    ks = [i for i, kk in enumerate(CAV['keys']) if kk.startswith('Ag|')]; i = max(ks, key=lambda i: E[i, 0])
    out['cavity'] = dict(cdA=float(E[i, 0] / np.sin(np.radians(5)) ** 2 * k), EQE=float(E[i, -1]), geom=str(CAV['keys'][i]))
    return out


if __name__ == '__main__':
    with Pool(4) as p:
        res = p.map(one, list(itertools.product(FW, TAIL)))
    json.dump(res, open('data/fwhm_tail.json', 'w'), indent=1)
    print('FWHM tail | DBR cd/A EQE duv30 | ideal20 cd/A | ideal10 cd/A | cavity cd/A | MLA cd/A | xy(DBR)')
    for r in res:
        print(f"{r['fwhm']:3d} {r['tail']:.2f} | {r['dbr']['cdA']:6.0f} {r['dbr']['EQE']:.3f} {r['dbr']['duv30']:.4f} | {r['ideal20']['cdA']:6.0f} | {r['ideal10']['cdA']:6.0f} | "
              f"{r['cavity']['cdA']:5.0f} | {r['MLA']['cdA']:4.0f} | ({r['dbr']['xy'][0]:.3f}, {r['dbr']['xy'][1]:.3f})")

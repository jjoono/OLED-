"""Joint fit of the green-device cavity (ETL = B3PyMPM, HTL = TAPC thickness) to bare, hemispherical and freeform
devices at once: tilt-corrected angular intensity (bare 10-80 deg, lens devices 10-60 deg) and Lambertian-equivalent
on-axis EQE with one common internal efficiency eta (least squares). 20 nm spectral grid, N = 40000 rays per wavelength.
Usage: python3 green_joint_fit.py targets.json out.json"""
import sys, json, itertools, numpy as np
from multiprocessing import Pool
from scipy.interpolate import CubicSpline
from scipy.integrate import trapezoid
sys.path.insert(0, '../selectivity'); import sel
from sel import cps2
import device as DV, simulate as SM, layout as LY
LAMS = np.arange(440, 641, 20).astype(float); J = (LAMS - 400).astype(int)
E = np.arange(91); dOm = 2 * np.pi * (np.cos(np.radians(E[:-1])) - np.cos(np.radians(E[1:])))
TH_H, Z0, Q = 0.77, 12.5, 0.98


def rel(prof):
    p = np.array(prof); i0 = (p[:5] * dOm[:5]).sum() / dOm[:5].sum()
    return np.array([1.0] + [(p[a - 3:a + 3] * dOm[a - 3:a + 3]).sum() / dOm[a - 3:a + 3].sum() / i0 for a in range(10, 90, 10)]), np.pi * i0


def stacks(etl, htl):
    """Per-wavelength source angular distribution P, back reflectances R (pixel + 4 regions) and spectral weights wq."""
    s = sel.D['s_Ir(ppy)2acac'][J] * LAMS; s = s / s.sum(); out = []; Fs = []; fs = []
    for jj, (j, lam) in enumerate(zip(J, LAMS)):
        org = lambda m, d: (sel.D['n_' + m + '_o'][j].real, sel.D['n_' + m + '_e'][j].real, d)
        met = lambda m, d: (sel.D['n_' + m][j], sel.D['n_' + m][j], d)
        eml = (sel.D['n_EML_o'][j].real, sel.D['n_EML_e'][j].real)
        cath = [org('B3', etl), met('LiF', 1.0), met('Al', 100.0)]
        S = cps2.Stack(lam, eml, 25.0, Z0, above=cath, below=[org('TCTA', 10.0), org('TAPC', htl), met('ITO', 150.0)], n_sub=1.51, h=TH_H)
        r = cps2.solve_pol(S, npts=3000); ref = cps2.solve_pol(cps2.Stack(lam, eml, 25.0, Z0, above=[], below=[], n_sub=eml[0], n_top=eml[0], h=TH_H), npts=3000)
        th_f = np.linspace(0, 90, 1801)[:-1] + 0.025; pf = np.clip(np.real(cps2.sub_angular(S, th_f)), 0, None).reshape(90, 20).sum(1); P = pf / pf.sum()
        R = {'pix': np.clip(np.real(cps2.stack_reflectance(S, DV.TH)), 0, 1)}
        for ito in (True, False):
            for al in (True, False):
                S2 = cps2.Stack(lam, eml, 25.0, Z0, above=cath if al else [org('B3', etl)], below=[org('TCTA', 10.0), org('TAPC', htl)] + ([met('ITO', 150.0)] if ito else []), n_sub=1.51, n_top=1.0, h=TH_H)
                R[('ITO' if ito else 'noITO') + '_' + ('Al' if al else 'noAl')] = np.clip(np.real(cps2.stack_reflectance(S2, DV.TH)), 0, 1)
        Fs.append(r['P_tot'] / ref['P_tot']); fs.append(r['air'] + r['sub']); out.append((P, R))
    Fs = np.array(Fs); fs = np.array(fs); wq = Q * s * Fs * fs / (Q * np.sum(s * Fs) + 1 - Q)
    return out, wq


def run(args):
    etl, htl = args[:2]; hk = args[2] if len(args) > 2 else 'hemisphere'
    PL = SM.planar_kernel(); KER = {'bare': PL, 'hemisphere': SM.lens_kernel(hk), 'freeform': SM.lens_kernel('green_freeform')}
    BTn = lambda K: K['BT'] / np.maximum(K['p'][:, None], 1e-15)
    ST, wq = stacks(etl, htl); air = {k: np.zeros((len(J), 90)) for k in KER}
    for jj, (P, R) in enumerate(ST):
        for nm, (t, hp) in {'bare': (0.7, 0.0), 'hemisphere': (1.4, 1.5), 'freeform': (1.4, 1.5)}.items():
            eS, eP, sd = LY.mc(KER[nm], PL, R, P, t, hp, True, N=40000, seed=jj)
            air[nm][jj] = eS @ BTn(KER[nm]) + eP @ BTn(PL)
    out = dict(etl=etl, htl=htl, hemi_kernel=hk)
    for nm in KER:
        A = wq[:, None] * air[nm]; rl, fr = rel(A.sum(0) / dOm); out[nm] = dict(EQE=float(A.sum()), front=float(fr), PE=float(A.sum() / fr), rel=rl.tolist())
    return out


if __name__ == '__main__':
    T = json.load(open(sys.argv[1]))
    grid = list(itertools.product((35., 40., 45., 50., 55., 60., 65., 70., 75., 80.), (80., 90., 100., 110., 120., 130.)))
    with Pool(4) as pool:
        res = []
        for r in pool.imap_unordered(run, grid):
            res.append(r); print(f"ETL {r['etl']:.0f} HTL {r['htl']:.0f}: " + ' | '.join(f"{k} front {r[k]['front']*100:.1f} PE {r[k]['PE']:.2f}" for k in ('bare', 'hemisphere', 'freeform')), flush=True)
            json.dump(res, open(sys.argv[2], 'w'))

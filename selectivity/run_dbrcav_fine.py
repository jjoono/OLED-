"""High-Q dielectric-mirror microcavity, resolved properly: 0.5 nm wavelength grid over the emission band,
HTL scanned in 2 nm steps (ETL from the coarse scan) so the resonance is placed on the emission, then
on-axis cd/A, luminance and u'v' versus angle and total EQE for the best HTL of each pair count."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import json, time, sys, numpy as np
from multiprocessing import Pool
import sel, common, run_display as RD
from sel import cps2

LF = np.arange(400.0, 700.01, 0.5)
L1 = np.arange(400, 701)
CIE = np.load('data/cie.npz')
XB, YB, ZB = (np.interp(LF, L1, CIE[k]) for k in 'xyz')
NKF = {k[2:]: np.interp(LF, L1, sel.D[k].real) + 1j * np.interp(LF, L1, sel.D[k].imag) for k in sel.D if k.startswith('n_')}


def spec_fine(name):
    if name.startswith('G'):
        f = float(name[1:]); s = np.exp(-4 * np.log(2) * ((LF - 530.0) / f) ** 2)
    else:
        s = np.interp(LF, L1, sel.D['s_' + name]) * LF
    return s / s.sum()


def stack(N, e, h, i, lam_d=530.0):
    nk = lambda m: NKF[m][i]; org = lambda m, d: (nk(m + '_o').real, nk(m + '_e').real, d)
    jd = int(lam_d) - 400; nH, nL = sel.nk('TiO2', jd).real, sel.nk('SiO2', jd).real
    dbr = []
    for _ in range(N):
        dbr += [(nk('SiO2'), nk('SiO2'), lam_d / 4 / nL), (nk('TiO2'), nk('TiO2'), lam_d / 4 / nH)]
    below = [org('TCTA', 10.0), org('TAPC', h), (nk('ITO'), nk('ITO'), 50.0)] + dbr
    return cps2.Stack(float(LF[i]), (nk('EML_o').real, nk('EML_e').real), 25.0, 12.5,
                      above=[org('B3', e), (nk('Ag'), nk('Ag'), 100.0)], below=below, n_sub=1.51, h=0.865)


def air_spectrum(args):
    """Per generated photon: power exiting into air per substrate bin (planar glass, angle preserving)."""
    N, e, h, idx = args
    out = np.zeros((len(idx), sel.NB)); Rf, Tf = sel.tmm_RT([], 530.0, sel.TH, 1.51)
    for k, i in enumerate(idx):
        S = stack(N, e, h, i)
        Pt = cps2.solve_pol(S, npts=800)['P_tot']
        th_f = np.linspace(0, 90, 1801)[:-1] + 0.025
        pf = np.clip(np.real(cps2.sub_angular(S, th_f)), 0, None).reshape(sel.NB, 10).sum(1) * 0.05 / Pt
        R = np.clip(np.real(cps2.stack_reflectance(S, sel.TH)), 0, 1)
        out[k] = Tf * pf / (1 - Rf * R)
    return args[:3], out


def metrics(air, idx, spec):
    I = RD.to_air_angle(air, 1.51)                         # (lam, angle) per generated photon per sr
    S = spec[idx][:, None] * I; lam = LF[idx]
    lum = 683 * (S * (YB[idx] * 1239.84 / lam)[:, None]).sum(0) / np.cos(np.radians(RD.AC))
    X, Y, Z = ((S * c[idx][:, None] / lam[:, None]).sum(0) for c in (XB, YB, ZB))
    den = X + 15 * Y + 3 * Z; u, v = 4 * X / den, 9 * Y / den; duv = np.hypot(u - u[0], v - v[0])
    half = float(np.interp(0.5, lum[::-1] / lum[0], RD.AC[::-1])) if lum.min() < 0.5 * lum[0] else 90.0
    return dict(cdA0=float(lum[0]), EQE=float((S * RD.SR[None, :]).sum()), half_lum_angle=half,
                duv30=float(duv[RD.AC <= 30].max()), duv60=float(duv[RD.AC <= 60].max()),
                lum_rel=(lum / lum[0]).tolist(), duv=duv.tolist())


if __name__ == '__main__':
    t0 = time.time(); z = np.load('data/display_scan.npz', allow_pickle=True); res = list(z['res'])
    out = {}
    for sp in sys.argv[1:] or ('G30', 'Ir(ppy)2acac', 'G10'):
        spec = spec_fine(sp); idx = np.where(spec > 1e-3 * spec.max())[0]
        s4 = common.SPECTRA[sp][::2]
        for N in (3, 6, 10):
            e = max((r for r in res if r[0][0] == N), key=lambda r: s4 @ (r[1] * r[2][:, 0]))[0][1]
            H = np.arange(20.0, 201.0, 2.0)
            onax = {}
            F5 = sel.cone_frac(1.51)[0]
            with Pool(4) as p:
                for (n_, e_, h), air in p.imap_unordered(air_spectrum, [(N, e, h, idx[::2]) for h in H]):
                    onax[h] = float((spec[idx[::2]] * (air @ F5)).sum() / spec[idx[::2]].sum())
            hb = max(onax, key=onax.get)
            hs = [hh for hh in H if abs(hh - hb) <= 4]
            with Pool(4) as p:
                fine = dict(((r[0][2], r[1]) for r in p.imap_unordered(air_spectrum, [(N, e, hh, idx) for hh in hs])))
            best = max(fine, key=lambda hh: metrics(fine[hh], idx, spec)['cdA0'])
            m = metrics(fine[best], idx, spec); m.update(N=N, ETL=e, HTL=float(best))
            out[f'{sp}|DBR cavity N{N}'] = m
            print(f"{sp:13s} N {N:2d} ETL {e:.0f} HTL {best:.0f}: cd/A {m['cdA0']:6.0f}  EQE {m['EQE']:.3f}  half-lum {m['half_lum_angle']:.1f}  du'v' <=30 {m['duv30']:.4f} <=60 {m['duv60']:.4f}  [{time.time() - t0:.0f}s]", flush=True)
    old = json.load(open('data/display_dbrcav.json')) if os.path.exists('data/display_dbrcav.json') else {}
    old.update(out); json.dump(old, open('data/display_dbrcav.json', 'w'))

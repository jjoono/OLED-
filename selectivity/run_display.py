"""Display-relevant comparison: on-axis cd/A is not the whole story.

(A) Strong dielectric-mirror microcavity: glass | (SiO2/TiO2)^N quarter-wave at the emission peak |
    ITO 50 | TAPC | TCTA 10 | EML 25 | B3PyMPM | Ag 100, planar exit, ETL/HTL scanned for maximum
    on-axis luminance (classic DBR-cavity OLED, many pairs).
(B) For each configuration: on-axis cd/A, luminance and colour (CIE 1976 u'v') versus viewing
    angle, total EQE.  Spectral radiant intensity per generated photon, I(lam, theta_air)."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import json, time, numpy as np
from multiprocessing import Pool
import sel, common, cps2

CIE = np.load('data/cie.npz'); XB, YB, ZB = (CIE[k][sel.IL] for k in 'xyz')
DV = dict(np.load('data/devices.npz')); SL = dict(np.load('data/slabs.npz'))
LAMC = 530.0
AE = np.arange(0, 91, 2.5); AC = 0.5 * (AE[:-1] + AE[1:]); SR = 2 * np.pi * (np.cos(np.radians(AE[:-1])) - np.cos(np.radians(AE[1:])))
SPECS = ('G30', 'G10', 'Ir(ppy)2acac')


def cav_device(N, d_etl, d_htl, lam_d=LAMC, lam_grid=None):
    lam_grid = sel.LAM if lam_grid is None else lam_grid
    es, P, R = [], [], []
    for lam in lam_grid:
        j = int(lam) - 400; nk = lambda m: sel.nk(m, j); org = lambda m, d: (nk(m + '_o').real, nk(m + '_e').real, d)
        jd = int(lam_d) - 400; nH, nL = sel.nk('TiO2', jd).real, sel.nk('SiO2', jd).real
        dbr = []
        for _ in range(N):                                   # listed outwards from the EML: ITO, then L H L H ... to glass
            dbr += [(nk('SiO2'), nk('SiO2'), lam_d / 4 / nL), (nk('TiO2'), nk('TiO2'), lam_d / 4 / nH)]
        below = [org('TCTA', 10.0), org('TAPC', d_htl), (nk('ITO'), nk('ITO'), 50.0)] + dbr
        S = cps2.Stack(float(lam), (nk('EML_o').real, nk('EML_e').real), 25.0, 12.5,
                       above=[org('B3', d_etl), (nk('Ag'), nk('Ag'), 100.0)], below=below, n_sub=1.51, h=0.865)
        r = cps2.solve_pol(S, npts=2000)
        th_f = np.linspace(0, 90, 1801)[:-1] + 0.025
        pf = np.clip(np.real(cps2.sub_angular(S, th_f)), 0, None).reshape(sel.NB, 10).sum(1)
        es.append(r['air'] + r['sub']); P.append(pf / max(pf.sum(), 1e-30))
        R.append(np.clip(np.real(cps2.stack_reflectance(S, sel.TH)), 0, 1))
    return dict(eta_sub=np.array(es), P=np.array(P), R=np.array(R), n_sub=1.51)


def planar_air(dev):
    Rf, Tf = sel.filter_RT('fresnel', 1.51)
    return Tf * dev['P'] / (1 - Rf * dev['R'])


def scan_one(args):
    N, e, h = args
    lg = sel.LAM[::2]                                        # 4 nm for the scan
    d = cav_device(N, e, h, lam_grid=lg)
    Rf, Tf = sel.tmm_RT([], 530.0, sel.TH, 1.51)
    air = Tf[None] * d['P'] / (1 - Rf[None] * d['R'])
    F = sel.cone_frac(1.51)
    return args, d['eta_sub'], air @ F[[0, -1]].T          # per lambda: (on-axis 5 deg, total)


def rebin(power, lo, hi):
    """Distribute power held in air-angle bins [lo, hi] (deg) onto AE, uniformly in sin^2 (Lambertian within a bin)."""
    s2 = lambda a: np.sin(np.radians(np.clip(a, 0, 90))) ** 2
    W = np.clip(np.minimum(s2(hi)[:, None], s2(AE[1:])[None, :]) - np.maximum(s2(lo)[:, None], s2(AE[:-1])[None, :]), 0, None)
    W = W / np.maximum(W.sum(1, keepdims=True), 1e-30)
    return power @ W


def to_air_angle(air_sub, n_sub):
    """(lam, substrate bin) power -> spectral radiant intensity (lam, air-angle bin) per sr."""
    lo = np.degrees(np.arcsin(np.clip(n_sub * np.sin(np.radians(sel.EDG[:-1])), 0, 1)))
    hi = np.degrees(np.arcsin(np.clip(n_sub * np.sin(np.radians(sel.EDG[1:])), 0, 1)))
    ok = n_sub * np.sin(np.radians(sel.EDG[:-1])) < 1
    return rebin(air_sub[:, ok], lo[ok], hi[ok]) / SR[None, :]


def display_metrics(I, eta_sub, spec):
    """I: (lam, angle) per substrate-delivered photon per sr.  Returns cd/A(theta), u'v'(theta), EQE."""
    S = spec[:, None] * eta_sub[:, None] * I                 # photons per generated photon per sr
    lum = 683 * (S * (YB * 1239.84 / sel.LAM)[:, None]).sum(0) / np.cos(np.radians(AC))
    X, Y, Z = ((S * c[:, None] / sel.LAM[:, None]).sum(0) for c in (XB, YB, ZB))     # energy-weighted tristimulus
    den = np.maximum(X + 15 * Y + 3 * Z, 1e-30); u = 4 * X / den; v = 9 * Y / den
    u = np.where(Y > 1e-6 * Y.max(), u, np.nan); v = np.where(Y > 1e-6 * Y.max(), v, np.nan)
    eqe = float((S * SR[None, :]).sum())
    duv = np.hypot(u - u[0], v - v[0])
    half = float(np.interp(0.5, lum[::-1] / lum[0], AC[::-1])) if lum.min() < 0.5 * lum[0] else 90.0
    return dict(cdA0=float(lum[0]), lum_rel=(lum / lum[0]).tolist(), duv=duv.tolist(), EQE=eqe, half_lum_angle=half,
                duv30=float(np.nanmax(duv[AC <= 30])), duv60=float(np.nanmax(duv[AC <= 60])),
                cone20=float((S * SR[None, :])[:, AC < 20].sum()))


if __name__ == '__main__':
    t0 = time.time()
    jobs = [(N, e, h) for N in (3, 6, 10) for e in np.arange(30., 301., 30.) for h in np.arange(30., 301., 30.)]
    if os.path.exists('data/display_scan.npz'):
        z = np.load('data/display_scan.npz', allow_pickle=True); res = list(z['res'])
    else:
        with Pool(2) as p:
            res = p.map(scan_one, jobs, chunksize=4)
        np.savez('data/display_scan.npz', res=np.array(res, dtype=object))
    print('scan done', time.time() - t0, flush=True)
    out = {'angles': AC.tolist(), 'configs': {}}
    lg = sel.LAM[::2]
    for sp in SPECS:
        spec4 = common.SPECTRA[sp][::2]; spec4 = spec4 / spec4.sum()
        for N in (3, 6, 10):
            best = max((r for r in res if r[0][0] == N), key=lambda r: spec4 @ (r[1] * r[2][:, 0]))
            (N_, e, h) = best[0]
            d = cav_device(N_, e, h)
            m = display_metrics(to_air_angle(planar_air(d), 1.51), d['eta_sub'], common.SPECTRA[sp])
            out['configs'][f'{sp}|DBR cavity N{N_} (ETL {e:.0f}, HTL {h:.0f})'] = m
            print(sp, N_, e, h, round(m['cdA0']), round(m['EQE'], 3), round(m['duv60'], 4), f'[{time.time() - t0:.0f}s]', flush=True)
    json.dump(out, open('data/display_cavity.json', 'w'))

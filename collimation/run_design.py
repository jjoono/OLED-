"""Design study: thin-film collimator (scatterer + DBR) on a weak-cavity top-emitting micro-OLED,
versus a conventional strong-cavity micro-OLED, evaluated as power delivered into the
acceptance cone (half-angle beta) of a pancake lens.

Outputs results.json and printed tables.  Green emitter centred at 530 nm, Gaussian spectra.
"""
import json, itertools, time, numpy as np
import model as Mo

LAM = np.arange(470.0, 601.0, 5.0)
FWHMS = [10.0, 20.0, 30.0, 50.0]
BETAS = [10.0, 15.0, 20.0, 25.0, 30.0]
A_SC = [0.0, 0.01, 0.02]                    # scatterer loss per traversal
MATS = {'TiO2/SiO2': (2.40, 1.46), 'Ta2O5/SiO2': (2.10, 1.46), 'SiNx/SiO2': (1.95, 1.46)}
t0 = time.time()


def spec_eval(dev, hist_fn):
    """dev: dict lam->(eta_sub,P,R); hist_fn(lam, P, R) -> air hist per unit delivered power.
    Returns per-FWHM air histogram per unit generated power."""
    H = {lam: dev[lam][0] * hist_fn(lam, dev[lam][1], dev[lam][2]) for lam in LAM}
    return {f: sum(Mo.spectrum(LAM, 530.0, f)[i] * H[lam] for i, lam in enumerate(LAM)) for f in FWHMS}


# ---------- 1. weak-cavity device: maximise eta_sub at 530 nm ----------
best = None
for de, dr in itertools.product(range(30, 301, 15), range(30, 301, 15)):
    es, P, R = Mo.device(530.0, de, dr, cathode='ito')
    if best is None or es > best[0]:
        best = (es, de, dr)
WEAK = dict(d_exit=best[1], d_refl=best[2])
devW = {lam: Mo.device(lam, WEAK['d_exit'], WEAK['d_refl'], cathode='ito') for lam in LAM}
Rbar = {lam: float((devW[lam][2] * Mo.WN).sum()) for lam in LAM}
print(f'weak cavity: ETL {best[1]} / HTL {best[2]} nm, eta_sub(530) {best[0]:.3f}, '
      f"A'(530) {1 - Rbar[530.0]:.4f}  [{time.time() - t0:.0f}s]", flush=True)

# ---------- 2. strong-cavity planar reference, optimised per beta ----------
TF = Mo.fresnel_T()
grid = list(itertools.product(range(20, 201, 10), range(20, 251, 10)))
single = {}
for de, dr in grid:
    es, P, R = Mo.device(530.0, de, dr, cathode='ag')
    single[(de, dr)] = es * Mo.planar(P, R, TF)
REF = {}
for b in BETAS:
    de, dr = max(grid, key=lambda g: Mo.cone(single[g], b))
    devS = {lam: Mo.device(lam, de, dr, cathode='ag') for lam in LAM}
    REF[b] = dict(d=(de, dr), hist=spec_eval(devS, lambda lam, P, R: Mo.planar(P, R, TF)))
print(f'strong-cavity references done [{time.time() - t0:.0f}s]', flush=True)

# weak cavity + scatterer only (no DBR) and weak planar (no scatterer)
SCAT_ONLY = spec_eval(devW, lambda lam, P, R: Mo.multipass(P, R, TF, 1.0, 0.0))
WEAK_PLANAR = spec_eval(devW, lambda lam, P, R: Mo.planar(P, R, TF))


# ---------- 3. DBR search (full randomisation -> closed form per wavelength) ----------
def closed_hist(T, Rb, a):
    p = float((T * Mo.WN).sum()); q = (1 - a) ** 2 * Rb
    shape = T * Mo.WN
    tot = p / (1 - (1 - p) * q)
    h = np.bincount(np.clip(np.searchsorted(Mo.AIR_EDGES, Mo.TH_AIR[Mo.ESC]) - 1, 0, 359),
                    weights=shape[Mo.ESC], minlength=360)
    return h / max(h.sum(), 1e-12) * tot


designs = []
for mat, (nH, nL) in MATS.items():
    for N in range(3, 17):
        for lb in np.arange(560.0, 801.0, 5.0):
            ns, ds = Mo.dbr_stack(lb, N, nH, nL)
            designs.append((mat, N, lb, {lam: Mo.tmm_T(ns, ds, lam, Mo.THR) for lam in LAM}))
print(f'{len(designs)} DBR designs traced [{time.time() - t0:.0f}s]', flush=True)

RES = {'weak': WEAK, 'beta': BETAS, 'fwhm': FWHMS, 'a_sc': A_SC, 'rows': []}
for f in FWHMS:
    w = Mo.spectrum(LAM, 530.0, f)
    for a in A_SC:
        # per design: spectrum-weighted air histogram per generated photon
        Hs = []
        for mat, N, lb, Tl in designs:
            h = sum(w[i] * devW[lam][0] * closed_hist(Tl[lam], Rbar[lam], a) for i, lam in enumerate(LAM))
            Hs.append(h)
        for b in BETAS:
            k = int(np.argmax([Mo.cone(h, b) for h in Hs]))
            h = Hs[k]; mat, N, lb, _ = designs[k]
            # ideal angular filter with cone = beta
            Tid = (Mo.TH_AIR < b) & Mo.ESC
            hid = sum(w[i] * devW[lam][0] * closed_hist(Tid.astype(float), Rbar[lam], a) for i, lam in enumerate(LAM))
            ref = REF[b]['hist'][f]; so = SCAT_ONLY[f]
            row = dict(fwhm=f, a_sc=a, beta=b, mat=mat, N=N, lam_b=float(lb),
                       eta_total=float(h.sum()), eta_beta=float(Mo.cone(h, b)), L0=float(Mo.axial_radiance(h)),
                       ideal_eta_beta=float(Mo.cone(hid, b)), ideal_total=float(hid.sum()),
                       ref_eta_beta=float(Mo.cone(ref, b)), ref_total=float(ref.sum()), ref_L0=float(Mo.axial_radiance(ref)),
                       ref_d=REF[b]['d'], scat_eta_beta=float(Mo.cone(so, b)), scat_total=float(so.sum()),
                       weakplanar_eta_beta=float(Mo.cone(WEAK_PLANAR[f], b)),
                       hist=h.tolist(), ref_hist=ref.tolist())
            RES['rows'].append(row)
            print(f'FWHM {f:4.0f} a_sc {a:.2f} beta {b:4.0f}: best {mat:10s} N={N:2d} lb={lb:.0f} | '
                  f'eta_beta {row["eta_beta"]:.3f} (ideal {row["ideal_eta_beta"]:.3f}, strong-cavity ref {row["ref_eta_beta"]:.3f}, '
                  f'scat-only {row["scat_eta_beta"]:.3f}) gain vs ref x{row["eta_beta"] / row["ref_eta_beta"]:.2f} | '
                  f'total {row["eta_total"]:.3f} vs ref {row["ref_total"]:.3f}', flush=True)
json.dump(RES, open('results.json', 'w'))
print(f'done [{time.time() - t0:.0f}s]')

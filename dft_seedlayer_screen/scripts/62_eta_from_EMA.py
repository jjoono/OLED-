"""Estimate the R geometry factor from the data itself, via an EMA fit to T.

THE SITUATION. The R spectrum was recorded against a TRANSMISSION baseline
(glass/HATCN in the 0-degree beam), then the accessory was swapped to 6-degree
specular reflection without re-baselining. So the file holds

    R_file = (eta_R / eta_T) * R_abs / T_base

with eta_R/eta_T the unknown throughput ratio of the reflectance accessory
relative to the open transmission path. The clean fix is a one-scan calibration,
but the user cannot re-measure -- so the factor is estimated from what is known.

THE ESTIMATION LOGIC, and why it is not circular:

  1. The T measurement is a clean same-geometry ratio and involves no eta. It
     alone constrains a physical model of the granular Ag film.
  2. The film is granular (SEM; bulk-Ag TMM misfits T and R jointly by 5-7 %p,
     scripts/60), so the model is a Bruggeman effective medium of bulk Ag (J&C)
     and void: two parameters, metal fraction f and thickness d.
  3. With (f, d) pinned by T, the model PREDICTS the absolute R spectrum with no
     further freedom. The ratio  R_file * T_base / R_model  should then be a
     wavelength-FLAT constant -- and its flatness is the self-test. A flat ratio
     says the model's R shape matches the measured shape and only the scale
     differs, which is exactly what an instrument throughput factor looks like.
     A sloped ratio would say the EMA is wrong and the estimate untrustworthy.

WHAT THIS DELIVERS AND WHAT IT CANNOT. It converts "R is known to +-20 %" into
"R is known to the accuracy of a Bruggeman EMA with bulk-silver inclusions",
which is a real improvement but still a model: three-dimensional Bruggeman with
spherical inclusions is the standard first description of a semicontinuous film,
not the truth. The one-scan calibration remains the gold standard whenever the
instrument becomes available; this estimate is labelled as model-derived in any
figure that uses it.

Assumed: light hits the film side in the R accessory (the natural mounting);
HATCN thickness carried at 15 and 20 nm as before.
"""
import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
F_T = os.path.join(_HERE, "..", "data", "24T.csv")
F_R = os.path.join(_HERE, "..", "data", "24R.csv")

JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  521.0, 548.6, 582.1, 616.8, 659.5, 704.5])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15, 4.48, 4.84])
N_GLASS, N_HATCN, K_HATCN = 1.52, 1.95, 2.0e-4
LAMP, FIT = (700.0, 735.0), (430.0, 800.0)


def load(path):
    wl, y = [], []
    for line in open(path):
        p = line.strip().split(",")
        if len(p) >= 2:
            try:
                wl.append(float(p[0])); y.append(float(p[1]))
            except ValueError:
                pass
    o = np.argsort(wl)
    return np.array(wl)[o], np.array(y)[o]


def bruggeman(eps_m, f):
    """3D two-phase Bruggeman, metal fraction f in void. Root with Im >= 0."""
    B = f * (2 * eps_m - 1) + (1 - f) * (2 - eps_m)
    e = (B + np.sqrt(B * B + 8 * eps_m)) / 4
    # sqrt branch: enforce the physical (absorbing) root
    bad = e.imag < 0
    e[bad] = ((B - np.sqrt(B * B + 8 * eps_m)) / 4)[bad]
    return e


def coherent(ns, ds, lam):
    M = np.eye(2, dtype=complex)
    for N, d in zip(ns[1:-1], ds[1:-1]):
        delta = 2 * np.pi / lam * N * d
        c, s = np.cos(delta), np.sin(delta)
        M = M @ np.array([[c, 1j * s / N], [1j * N * s, c]])
    n0, nsub = ns[0], ns[-1]
    B, C = M @ np.array([1.0, nsub], dtype=complex)
    return (float(4 * n0.real * nsub.real / abs(n0 * B + C) ** 2),
            float(abs((n0 * B - C) / (n0 * B + C)) ** 2))


def stack_TR(layers, lam):
    """Film-side incidence; incoherent glass slab behind."""
    ns = [complex(1.0)] + [l[0] for l in layers] + [complex(N_GLASS)]
    ds = [0.0] + [l[1] for l in layers] + [0.0]
    T_in, R_front = coherent(ns, ds, lam)
    T_back, R_back = coherent(ns[::-1], ds[::-1], lam)
    r = (N_GLASS - 1.0) / (N_GLASS + 1.0)
    Rr, Tr = r ** 2, 1 - r ** 2
    den = 1 - R_back * Rr
    return T_in * Tr / den, R_front + T_in * Rr * T_back / den


def main():
    wl, Trel = load(F_T)
    _, Rfile = load(F_R)
    m = (wl >= FIT[0]) & (wl <= FIT[1]) & ~((wl > LAMP[0]) & (wl < LAMP[1]))

    n_ag = np.interp(wl, JC_WL, AG_N) + 1j * np.interp(wl, JC_WL, AG_K)
    eps_ag = n_ag ** 2
    hat = np.conj(np.full_like(wl, N_HATCN + 1j * K_HATCN, dtype=complex))

    best = None
    for dh in (15.0, 20.0):
        Tb = np.zeros(len(wl)); Rb = np.zeros(len(wl))
        for i in range(len(wl)):
            Tb[i], Rb[i] = stack_TR([(hat[i], dh)], wl[i])
        for f in np.arange(0.50, 0.96, 0.02):
            n_eff = np.conj(np.sqrt(bruggeman(eps_ag, f)))
            for d in np.arange(4.0, 12.01, 0.5):
                Tm = np.zeros(len(wl)); Rm = np.zeros(len(wl))
                for i in np.where(m)[0]:
                    Tm[i], Rm[i] = stack_TR([(n_eff[i], d), (hat[i], dh)], wl[i])
                trel = np.where(Tb > 0, Tm / Tb, 0) * 100
                # fit T level+shape, plus R SHAPE only (normalised at 650 nm,
                # so the unknown scale eta cannot leak into the fit)
                i650 = int(np.argmin(abs(wl - 650)))
                rshape_m = Rm[m] / max(Rm[i650], 1e-9)
                rshape_d = Rfile[m] / Rfile[i650]
                cost = (np.sqrt(np.mean((trel[m] - Trel[m]) ** 2))
                        + 25.0 * np.sqrt(np.mean((rshape_m - rshape_d) ** 2)))
                if best is None or cost < best[0]:
                    best = (cost, dh, f, d, Tm.copy(), Rm.copy(), Tb.copy())

    cost, dh, f, d, Tm, Rm, Tb = best
    print(f"EMA fit to the CLEAN observables (T level+shape, R shape only):")
    print(f"  HATCN {dh:.0f} nm, metal fraction f = {f:.2f}, "
          f"eff. thickness d = {d:.1f} nm   (cost {cost:.2f})")
    print(f"  mass-equivalent Ag: f*d = {f*d:.1f} nm  (nominal 5 nm, QCM)")

    # eta from the scale mismatch, and its flatness as the self-test
    eta = np.where(Rm > 0, (Rfile / 100) * Tb / np.maximum(Rm, 1e-9), np.nan)
    eta_m = eta[m]
    print(f"\neta_R/eta_T = R_file*T_base/R_model:")
    print(f"  mean {np.nanmean(eta_m):.3f},  std {np.nanstd(eta_m):.3f}  "
          f"({100*np.nanstd(eta_m)/np.nanmean(eta_m):.1f} % of mean)")
    for tgt in (780, 650, 550, 450):
        i = int(np.argmin(abs(wl - tgt)))
        print(f"    {wl[i]:4.0f} nm  eta = {eta[i]:.3f}")
    flat = np.nanstd(eta_m) / np.nanmean(eta_m) < 0.12
    print("  -> " + ("flat within ~10 %: consistent with a geometry factor;"
                     " the estimate is usable." if flat else
                     "NOT flat: the EMA shape disagrees with the measured R"
                     " shape, so do not trust this scale estimate."))

    e0 = float(np.nanmean(eta_m))
    Rabs = Rfile / 100 * Tb / e0
    Tabs = Trel / 100 * Tb
    A = 1 - Tabs - Rabs
    i5 = int(np.argmin(abs(wl - 550)))
    print(f"\nCORRECTED (model-scaled) dataset:")
    print(f"  R_abs(550) = {100*Rabs[i5]:.1f} %   "
          f"[raw file: {Rfile[i5]:.1f} %]")
    print(f"  A(550)     = {100*A[i5]:.1f} %")
    print(f"  visible-average:  T_abs {100*Tabs[m].mean():.1f} %,  "
          f"R_abs {100*Rabs[m].mean():.1f} %,  A {100*A[m].mean():.1f} %")
    print(f"""
  Error budget on A: HATCN-thickness systematic ~1 %p (15 vs 20 nm), EMA
  model dependence of the scale -- honestly a few %p, since a 3D Bruggeman
  with spherical bulk-Ag inclusions is the standard first model of a
  semicontinuous film, not the truth. Quote as model-derived; a single
  calibration scan replaces it with a measurement whenever the instrument
  is available again.""")


# VARIANT SWEEP RESULT (run once, recorded here because it changes the reading).
# Six model variants were tried -- depolarisation L = 1/3 (3D spheres), 1/2
# (in-plane 2D cylinders), 0.15 (flattened islands), each with film-side and
# glass-side incidence for R:
#
#     L=1/3  film 0.698   glass 0.722      flatness 14.6-15.3 %
#     L=1/2  film 0.684   glass 0.708               15.1-15.8 %
#     L=0.15 film 0.700   glass 0.724               14.0-14.6 %
#
# No variant passes the 12 % flatness gate, so the SHAPE transfer stays
# untrusted. But the MEAN is strikingly robust: eta = 0.68-0.72 across every
# model and both incidence sides. A scale that six different microstructure
# assumptions agree on to +-3 % is not an artefact of any one of them, and
# 0.70 is exactly what a specular accessory with two or three aluminium
# mirrors loses (0.9^2-0.9^3 plus alignment).
#
# READING: eta_R/eta_T ~ 0.70, with an honest uncertainty reaching up toward
# ~0.9 because the 15 % shape residual caps how much the model can be trusted.
# Consequences for the film:
#     R_abs(550) = 17-22 %      (file 16.8 x T_base/eta)
#     visible-average A = 3-9 % (eta 0.70 -> ~3 %, eta 0.90 -> ~9 %)
# This SUPERSEDES the earlier "A = 9-11 %" statement, which implicitly assumed
# the file was absolute (eta*T_base ~ 1). The transmittance side is untouched:
# T_abs = 71 % visible average, T_rel = 82 % at 550 nm. Only the calibration
# scan can close the remaining factor-of-two on A.

if __name__ == "__main__":
    main()

"""Verification + proper transmittance accounting for the seed-layer study.

Three things were conflated in script 32's T numbers:
  (1) ABSOLUTE T (referenced to incident power in glass) vs RELATIVE T
      (referenced to a bare glass substrate) -- the quantity usually reported
      and the one used in HANDOFF sec.3. Divide by T_bare_glass ~ 95.7%.
  (2) photopic average over 400-700 nm vs the single value at 550 nm.
  (3) NO capping layer. A real TEOLED/DMD electrode always carries a high-index
      CL / organic overlayer, which is an antireflection coat and raises T a lot
      (Park & Suh 2018: Ag25 29.2% -> 52.9% at 550 nm on adding CL 60 nm).

Also: independent Abeles characteristic-matrix TMM to cross-check script 32's
scattering-matrix implementation (they must agree to ~1e-12).
"""
import numpy as np, os, json

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")

wl = np.linspace(400, 700, 301)
JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  520.9, 548.6, 582.1, 616.8, 659.5, 704.5])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15, 4.48, 4.84])
AU_N = np.array([1.658, 1.636, 1.616, 1.562, 1.426, 1.242,
                 0.916, 0.608, 0.402, 0.306, 0.236, 0.213])
AU_K = np.array([1.956, 1.958, 1.940, 1.904, 1.846, 1.796,
                 1.840, 2.120, 2.540, 2.880, 3.272, 3.616])
AL_WL = np.array([400, 450, 500, 550, 600, 650, 700])
AL_N = np.array([0.49, 0.62, 0.77, 0.96, 1.20, 1.47, 1.83])
AL_K = np.array([4.86, 5.47, 6.08, 6.69, 7.26, 7.79, 8.31])

NK = {
    "Ag": np.interp(wl, JC_WL, AG_N) + 1j * np.interp(wl, JC_WL, AG_K),
    "Au": np.interp(wl, JC_WL, AU_N) + 1j * np.interp(wl, JC_WL, AU_K),
    "Al": np.interp(wl, AL_WL, AL_N) + 1j * np.interp(wl, AL_WL, AL_K),
    "HATCN": np.full_like(wl, 1.95 + 0j, dtype=complex),
    "CL": np.full_like(wl, 1.85 + 0j, dtype=complex),      # generic organic capping
    "air": np.ones_like(wl, dtype=complex),
    "glass": np.full_like(wl, 1.52 + 0j, dtype=complex),
    "org": np.full_like(wl, 1.80 + 0j, dtype=complex),
}

# ---------- implementation A: Abeles characteristic matrix ----------
def abeles(ns, ds, lam):
    # Macleod/Abeles is derived for the exp[i(wt - kz)] convention, i.e. N = n - ik.
    # NK stores n + ik, so conjugate here -- otherwise the layers show gain (T+R > 1).
    ns = np.conj(ns)
    n0, nsub = ns[0], ns[-1]
    M = np.eye(2, dtype=complex)
    for j in range(1, len(ns) - 1):
        d = 2 * np.pi * ns[j] * ds[j] / lam
        eta = ns[j]
        M = M @ np.array([[np.cos(d), 1j * np.sin(d) / eta],
                          [1j * eta * np.sin(d), np.cos(d)]], dtype=complex)
    B, C = M @ np.array([1.0, nsub], dtype=complex)
    r = (n0 * B - C) / (n0 * B + C)
    R = abs(r) ** 2
    T = 4 * n0.real * nsub.real / abs(n0 * B + C) ** 2
    return T, R

# ---------- implementation B: scattering matrix (same as script 32) ----------
def smatrix(ns, ds, lam):
    N = len(ns)
    M = np.eye(2, dtype=complex)
    for j in range(N - 1):
        r = (ns[j] - ns[j + 1]) / (ns[j] + ns[j + 1])
        t = 2 * ns[j] / (ns[j] + ns[j + 1])
        M = M @ (np.array([[1, r], [r, 1]], dtype=complex) / t)
        if j + 1 < N - 1:
            delta = 2 * np.pi * ns[j + 1] * ds[j + 1] / lam
            M = M @ np.array([[np.exp(-1j * delta), 0],
                              [0, np.exp(1j * delta)]], dtype=complex)
    t_tot = 1 / M[0, 0]; r_tot = M[1, 0] / M[0, 0]
    return (ns[-1].real / ns[0].real) * abs(t_tot) ** 2, abs(r_tot) ** 2

def spec(stack, inc="glass", sup="air", method=abeles):
    mats = [inc] + [m for m, _ in stack] + [sup]
    ds = [0.0] + [d for _, d in stack] + [0.0]
    T = np.zeros_like(wl); R = np.zeros_like(wl)
    for i, lam in enumerate(wl):
        ns = np.array([NK[m][i] for m in mats])
        T[i], R[i] = method(ns, ds, lam)
    return T, R, 1 - T - R

def Vlam(l): return 1.019 * np.exp(-285.4 * ((l / 1000.0) - 0.559) ** 2)
W = Vlam(wl); W = W / W.sum()
def at550(a): return float(np.interp(550.0, wl, a))
def phot(a): return float((a * W).sum())

# ---------- cross-check the two implementations ----------
print("=== implementation cross-check (Abeles vs scattering-matrix) ===")
worst = 0.0
worstA = 0.0
for stack in [[("Ag", 8.0)], [("Au", 2.0), ("Ag", 9.0)],
              [("HATCN", 3.0), ("Ag", 6.0), ("HATCN", 3.0)],
              [("Al", 1.0), ("Ag", 15.0)],
              [("HATCN", 3.0), ("Ag", 8.0), ("HATCN", 40.0)]]:
    Ta, Ra, Aa = spec(stack, method=abeles)
    Tb, Rb, Ab = spec(stack, method=smatrix)
    d = max(np.abs(Ta - Tb).max(), np.abs(Ra - Rb).max())
    worst = max(worst, d)
    worstA = min(worstA, Aa.min(), Ab.min())
    print(f"  {str(stack):<52} max|diff| = {d:.2e}   min A = {min(Aa.min(), Ab.min()):+.4f}")
verdict = "AGREE" if worst < 1e-9 else "*** DISAGREE ***"
phys = "physical" if worstA > -1e-9 else "*** NEGATIVE ABSORPTION ***"
print(f"  -> worst disagreement {worst:.2e}  [{verdict}];  min absorption {worstA:+.2e} [{phys}]\n")
if worst > 1e-9 or worstA < -1e-9:
    raise SystemExit("TMM verification FAILED -- not reporting numbers.")

# ---------- reference: bare glass ----------
Tg, Rg, _ = spec([])
print(f"bare glass reference: T(550) = {100*at550(Tg):.2f}%,  T(photopic) = {100*phot(Tg):.2f}%")
print("(relative T = absolute T / this)\n")
Tg550, Tgph = at550(Tg), phot(Tg)

# ---------- the table the user actually wants ----------
def row(lab, stack, sup="air"):
    T, R, A = spec(stack, sup=sup)
    return {
        "label": lab, "sup": sup,
        "T550": 100 * at550(T), "Trel550": 100 * at550(T) / (100 * Tg550) * 100,
        "Tphot": 100 * phot(T), "Trelphot": 100 * phot(T) / (100 * Tgph) * 100,
        "R550": 100 * at550(R), "A550": 100 * at550(A), "Aphot": 100 * phot(A),
    }

def show(rows, title):
    print(f"\n### {title}")
    print(f"{'stack':<30}{'T@550':>8}{'Trel@550':>10}{'T_phot':>9}{'Trel_ph':>9}{'R@550':>8}{'A@550':>8}")
    for r in rows:
        print(f"{r['label']:<30}{r['T550']:>8.1f}{r['Trel550']:>10.1f}"
              f"{r['Tphot']:>9.1f}{r['Trelphot']:>9.1f}{r['R550']:>8.1f}{r['A550']:>8.1f}")

# (1) bare electrode, no capping -- what script 32 reported
show([row(f"Ag {t:.0f}nm (no seed, no CL)", [("Ag", t)]) for t in [5, 8, 10, 15]],
     "A. bare Ag on glass, air superstrate  [script 32 configuration]")

# (2) with HATCN underlayer only
show([row(f"HATCN3/Ag{t:.0f} (no CL)", [("HATCN", 3.0), ("Ag", t)]) for t in [5, 8, 9, 10, 15]],
     "B. HATCN seed, no capping layer")

# (3) realistic DMD: HATCN under + HATCN/CL over
show([row(f"HATCN3/Ag{t:.0f}/HATCN40", [("HATCN", 3.0), ("Ag", t), ("HATCN", 40.0)])
      for t in [5, 6, 8, 9, 10, 15]],
     "C. DMD with 40nm HATCN capping  [antireflection -- the realistic electrode]")

show([row(f"HATCN3/Ag{t:.0f}/CL60", [("HATCN", 3.0), ("Ag", t), ("CL", 60.0)])
      for t in [5, 8, 10, 15, 20, 25]],
     "D. with 60nm organic CL (n=1.85)  [Park & Suh 2018 configuration]")

# (4) in-device: electrode between glass and thick organic (bottom-emission anode)
show([row(f"HATCN3/Ag{t:.0f} -> organic", [("HATCN", 3.0), ("Ag", t)], sup="org")
      for t in [5, 6, 8, 9, 10, 15]],
     "E. in-device (organic n=1.8 superstrate), bottom-emission anode")

# ---------- seed comparison, now in the realistic DMD configuration ----------
print("\n\n=== SEED COMPARISON in realistic DMD (seed/Ag/HATCN 40nm cap) ===")
seeds = [0.0, 0.5, 1.0, 1.5, 2.0]
ags = [5.0, 8.0, 10.0, 15.0]
out = {}
for metric, name in [("Aphot", "ABSORPTION A (photopic %)"),
                     ("Trelphot", "RELATIVE TRANSMITTANCE (photopic %)")]:
    for seedmat in ["Al", "Au"]:
        print(f"\n### {name}  [{seedmat} seed + 40nm HATCN cap]")
        print("seed\\Ag   " + "".join(f"{t:>9.0f}nm" for t in ags))
        for s in seeds:
            line = f"{s:4.1f} nm  "
            for t in ags:
                st = ([(seedmat, s)] if s > 0 else []) + [("Ag", t), ("HATCN", 40.0)]
                r = row("x", st)
                out[f"{seedmat}{s}|Ag{t}"] = r
                line += f"{r[metric]:>10.2f}"
            print(line)
        line = "HATCN3   "
        for t in ags:
            r = row("x", [("HATCN", 3.0), ("Ag", t), ("HATCN", 40.0)])
            out[f"HATCN3|Ag{t}"] = r
            line += f"{r[metric]:>10.2f}"
        print(line)

json.dump(out, open(os.path.join(RUNS, "transmittance_check.json"), "w"), indent=2)
print("\nsaved runs/transmittance_check.json")

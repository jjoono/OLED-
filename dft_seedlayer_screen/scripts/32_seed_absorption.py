"""Absorption penalty of METAL seed layers (Al, Au) under thin Ag, vs HATCN seed.

Grid: seed = 0 / 0.5 / 1 / 1.5 / 2 nm  x  Ag = 5 / 8 / 10 / 15 nm.
TMM at normal incidence, 400-700 nm, with LAYER-RESOLVED absorption
(Poynting-vector difference across each layer) so the seed's own loss is separated
from the Ag loss.

nk sources (all standard bulk optical constants, 12-pt J&C grid interpolated):
  Ag : Johnson & Christy 1972      (same table as scripts 19-21 -> consistent)
  Au : Johnson & Christy 1972
  Al : Palik / Rakic
  HATCN : n = 1.95, k = 0 (transparent model; real HATCN absorbs only < ~450 nm)

Configurations:
  (a) measurable   : glass / seed / Ag / air
  (b) in-device    : glass / seed / Ag / organic(n=1.8)   [bottom-emission anode]

NOTE ON THE MODEL: bulk nk applied to a 0.5-2 nm layer is an idealization -- such a
layer is not a continuous film. Real sub-nm metal is discontinuous with extra
damping, so the numbers here are a LOWER BOUND on the metal-seed absorption penalty.
"""
import numpy as np, os, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RUNS = os.path.join(BASE, "runs")

wl = np.linspace(400, 700, 301)

# ---------------- optical constants ----------------
# Johnson & Christy 1972 (same wavelength grid for Ag and Au)
JC_WL = np.array([397.4, 413.3, 430.5, 450.9, 471.4, 495.9,
                  520.9, 548.6, 582.1, 616.8, 659.5, 704.5])
AG_N = np.array([0.05, 0.05, 0.04, 0.04, 0.05, 0.05, 0.05, 0.06, 0.05, 0.06, 0.05, 0.04])
AG_K = np.array([2.07, 2.21, 2.36, 2.66, 2.83, 3.09, 3.34, 3.59, 3.93, 4.15, 4.48, 4.84])
AU_N = np.array([1.658, 1.636, 1.616, 1.562, 1.426, 1.242,
                 0.916, 0.608, 0.402, 0.306, 0.236, 0.213])
AU_K = np.array([1.956, 1.958, 1.940, 1.904, 1.846, 1.796,
                 1.840, 2.120, 2.540, 2.880, 3.272, 3.616])
# Al (Palik/Rakic), coarser grid
AL_WL = np.array([400, 450, 500, 550, 600, 650, 700])
AL_N = np.array([0.49, 0.62, 0.77, 0.96, 1.20, 1.47, 1.83])
AL_K = np.array([4.86, 5.47, 6.08, 6.69, 7.26, 7.79, 8.31])

n_ag = np.interp(wl, JC_WL, AG_N) + 1j * np.interp(wl, JC_WL, AG_K)
n_au = np.interp(wl, JC_WL, AU_N) + 1j * np.interp(wl, JC_WL, AU_K)
n_al = np.interp(wl, AL_WL, AL_N) + 1j * np.interp(wl, AL_WL, AL_K)
n_hatcn = np.full_like(wl, 1.95 + 0j, dtype=complex)
n_air = np.ones_like(wl, dtype=complex)
n_glass = np.full_like(wl, 1.52 + 0j, dtype=complex)
n_org = np.full_like(wl, 1.80 + 0j, dtype=complex)

NK = {"Ag": n_ag, "Au": n_au, "Al": n_al, "HATCN": n_hatcn,
      "air": n_air, "glass": n_glass, "org": n_org}

# ---------------- TMM with layer-resolved absorption ----------------
def tmm_layers(ns, ds, lam):
    """ns: complex n per layer (incident..exit), ds: thickness nm (ends ignored).
    Returns T, R, and absorption fraction in each finite layer."""
    N = len(ns)
    kz = 2 * np.pi * ns / lam                      # 1/nm
    # interface + propagation transfer matrices; M relates (v0,w0) to (v_N, 0)
    Ms = []
    for j in range(N - 1):
        r = (ns[j] - ns[j + 1]) / (ns[j] + ns[j + 1])
        t = 2 * ns[j] / (ns[j] + ns[j + 1])
        I = np.array([[1, r], [r, 1]], dtype=complex) / t
        Ms.append(I)
    # full matrix
    M = np.eye(2, dtype=complex)
    for j in range(N - 1):
        M = M @ Ms[j]
        if j + 1 < N - 1:
            d = ds[j + 1]
            delta = kz[j + 1] * d
            P = np.array([[np.exp(-1j * delta), 0], [0, np.exp(1j * delta)]], dtype=complex)
            M = M @ P
    t_tot = 1 / M[0, 0]
    r_tot = M[1, 0] / M[0, 0]
    T = (ns[-1].real / ns[0].real) * abs(t_tot) ** 2
    R = abs(r_tot) ** 2

    # forward/backward amplitudes at TOP of each finite layer, by back-propagation
    # amplitudes in exit layer
    vw = np.array([t_tot, 0], dtype=complex)
    amps = [None] * N
    amps[N - 1] = vw.copy()
    for j in range(N - 2, 0, -1):
        # undo interface j->j+1 then propagation inside layer j
        vw = Ms[j] @ vw
        amps[j] = vw.copy()                        # at BOTTOM boundary of layer j...
        d = ds[j]
        delta = kz[j] * d
        P = np.array([[np.exp(-1j * delta), 0], [0, np.exp(1j * delta)]], dtype=complex)
        vw = P @ vw
        amps[j] = vw.copy()                        # at TOP boundary of layer j
    S0 = ns[0].real                                # incident Poynting (unit amplitude)

    def poyn(nj, v, w, kzj, z):
        f = v * np.exp(1j * kzj * z)
        b = w * np.exp(-1j * kzj * z)
        return np.real(np.conj(nj) * (f + b) * np.conj(f - b))

    A = np.zeros(N)
    for j in range(1, N - 1):
        v, w = amps[j]
        A[j] = (poyn(ns[j], v, w, kz[j], 0.0) - poyn(ns[j], v, w, kz[j], ds[j])) / S0
    return T, R, A

def spectrum(stack, sup="air"):
    """stack: list of (material, thickness_nm). Returns dict of arrays."""
    mats = ["glass"] + [m for m, _ in stack] + [sup]
    ds = [0.0] + [d for _, d in stack] + [0.0]
    T = np.zeros_like(wl); R = np.zeros_like(wl)
    Alay = np.zeros((len(wl), len(mats)))
    for i, lam in enumerate(wl):
        ns = np.array([NK[m][i] for m in mats])
        T[i], R[i], a = tmm_layers(ns, ds, lam)
        Alay[i] = a
    return {"T": T, "R": R, "A": 1 - T - R, "Alay": Alay, "mats": mats}

# photopic-ish weighting (CIE V(lambda) approximation, gaussian fit)
def Vlam(l):
    return 1.019 * np.exp(-285.4 * ((l / 1000.0) - 0.559) ** 2)
W = Vlam(wl); W = W / W.sum()

def at550(arr): return float(np.interp(550.0, wl, arr))
def photopic(arr): return float((arr * W).sum())

# ---------------- run the grid ----------------
seeds = [0.0, 0.5, 1.0, 1.5, 2.0]
ags = [5.0, 8.0, 10.0, 15.0]
results = {}
rows = []

for sup in ["air", "org"]:
    for seedmat in ["Al", "Au"]:
        for s in seeds:
            for t in ags:
                stack = ([(seedmat, s)] if s > 0 else []) + [("Ag", t)]
                r = spectrum(stack, sup)
                idx_seed = 1 if s > 0 else None
                idx_ag = 2 if s > 0 else 1
                a_seed = photopic(r["Alay"][:, idx_seed]) if idx_seed else 0.0
                a_ag = photopic(r["Alay"][:, idx_ag])
                key = f"{sup}|{seedmat}{s}|Ag{t}"
                results[key] = {
                    "A550": at550(r["A"]), "Aphot": photopic(r["A"]),
                    "T550": at550(r["T"]), "Tphot": photopic(r["T"]),
                    "R550": at550(r["R"]),
                    "A_seed_phot": a_seed, "A_Ag_phot": a_ag,
                }
                rows.append((sup, seedmat, s, t, results[key]))
        # HATCN reference for this superstrate
        for t in ags:
            r = spectrum([("HATCN", 3.0), ("Ag", t)], sup)
            key = f"{sup}|HATCN3|Ag{t}"
            results[key] = {
                "A550": at550(r["A"]), "Aphot": photopic(r["A"]),
                "T550": at550(r["T"]), "Tphot": photopic(r["T"]),
                "R550": at550(r["R"]),
                "A_seed_phot": photopic(r["Alay"][:, 1]),
                "A_Ag_phot": photopic(r["Alay"][:, 2]),
            }

json.dump(results, open(os.path.join(RUNS, "seed_absorption.json"), "w"), indent=2)

# ---------------- print tables ----------------
def table(sup, seedmat, metric, label):
    print(f"\n### {label}  [{seedmat} seed, superstrate={sup}]")
    print("seed\\Ag   " + "".join(f"{t:>9.0f}nm" for t in ags))
    for s in seeds:
        line = f"{s:4.1f} nm  "
        for t in ags:
            v = results[f"{sup}|{seedmat}{s}|Ag{t}"][metric]
            line += f"{100*v:>10.2f}"
        print(line)
    line = "HATCN3   "
    for t in ags:
        line += f"{100*results[f'{sup}|HATCN3|Ag{t}'][metric]:>10.2f}"
    print(line)

for sup in ["air", "org"]:
    for m in ["Al", "Au"]:
        table(sup, m, "Aphot", "ABSORPTION A (%, photopic-weighted 400-700nm)")
for sup in ["air"]:
    for m in ["Al", "Au"]:
        table(sup, m, "Tphot", "TRANSMITTANCE T (%, photopic)")

print("\n### Layer-resolved absorption (photopic %, superstrate=air, Ag 8nm)")
print(f"{'stack':<16}{'A_total':>10}{'A_seed':>10}{'A_Ag':>10}")
for m in ["Al", "Au"]:
    for s in seeds:
        r = results[f"air|{m}{s}|Ag8.0"]
        print(f"{m}{s}/Ag8{'':<6}{100*r['Aphot']:>10.2f}{100*r['A_seed_phot']:>10.2f}{100*r['A_Ag_phot']:>10.2f}")
r = results["air|HATCN3|Ag8.0"]
print(f"{'HATCN3/Ag8':<16}{100*r['Aphot']:>10.2f}{100*r['A_seed_phot']:>10.2f}{100*r['A_Ag_phot']:>10.2f}")

# ---------------- figure ----------------
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
for ax, (m, ttl) in zip(axes[:2], [("Al", "Al seed"), ("Au", "Au seed")]):
    for s in seeds:
        r = spectrum(([(m, s)] if s > 0 else []) + [("Ag", 8.0)], "air")
        ax.plot(wl, 100 * r["A"], label=f"{m} {s} nm" if s > 0 else "no seed (neat Ag)")
    rh = spectrum([("HATCN", 3.0), ("Ag", 8.0)], "air")
    ax.plot(wl, 100 * rh["A"], "k--", lw=2, label="HATCN 3 nm")
    ax.set_title(f"{ttl} / Ag 8 nm"); ax.set_xlabel("wavelength (nm)")
    ax.set_ylabel("absorption (%)"); ax.legend(fontsize=8); ax.grid(alpha=.3)

ax = axes[2]
for m, c in [("Al", "#c0392b"), ("Au", "#d4a017")]:
    for t, ls in zip(ags, ["-", "--", "-.", ":"]):
        ys = [100 * results[f"air|{m}{s}|Ag{t}"]["Aphot"] for s in seeds]
        ax.plot(seeds, ys, ls, color=c, marker="o", ms=3,
                label=f"{m}, Ag {t:.0f}nm")
ax.set_xlabel("seed thickness (nm)"); ax.set_ylabel("photopic absorption (%)")
ax.set_title("Absorption vs seed thickness"); ax.legend(fontsize=7, ncol=2); ax.grid(alpha=.3)
plt.tight_layout()
os.makedirs(os.path.join(BASE, "optics"), exist_ok=True)
plt.savefig(os.path.join(BASE, "optics", "seed_absorption.png"), dpi=140)
print("\nsaved seed_absorption.png + runs/seed_absorption.json")

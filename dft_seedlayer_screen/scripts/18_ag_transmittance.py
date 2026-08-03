"""Transfer-matrix transmittance of Ag films (8, 15 nm) on glass, 400-700 nm.
Ag n,k: Johnson & Christy (1972) bulk values, linearly interpolated.
Outputs: PNG plot + CSV. T_rel = T(film on glass) / T(bare glass).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
plt.rcParams["font.family"] = "Malgun Gothic"
plt.rcParams["axes.unicode_minus"] = False

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Johnson & Christy Ag (wavelength nm, n, k)
JC = np.array([
    [397.4, 0.05, 2.07], [413.3, 0.05, 2.21], [430.5, 0.04, 2.36],
    [450.9, 0.04, 2.66], [471.4, 0.05, 2.83], [495.9, 0.05, 3.09],
    [520.9, 0.05, 3.34], [548.6, 0.06, 3.59], [582.1, 0.05, 3.93],
    [616.8, 0.06, 4.15], [659.5, 0.05, 4.48], [704.5, 0.04, 4.84],
])

wl = np.linspace(400, 700, 301)
n_ag = np.interp(wl, JC[:, 0], JC[:, 1]) + 1j * np.interp(wl, JC[:, 0], JC[:, 2])
n_air, n_glass = 1.0, 1.52

def tmm_T(n_list, d_list, lam):
    """normal incidence; first/last media semi-infinite (d ignored)."""
    # interface + propagation matrices
    M = np.eye(2, dtype=complex)
    for j in range(len(n_list) - 1):
        n1, n2 = n_list[j], n_list[j + 1]
        r = (n1 - n2) / (n1 + n2)
        t = 2 * n1 / (n1 + n2)
        I = np.array([[1, r], [r, 1]], dtype=complex) / t
        if j + 1 < len(n_list) - 1:
            delta = 2 * np.pi * n_list[j + 1] * d_list[j + 1] / lam
            P = np.array([[np.exp(-1j * delta), 0], [0, np.exp(1j * delta)]])
            M = M @ I @ P
        else:
            M = M @ I
    t_tot = 1 / M[0, 0]
    return (np.real(n_list[-1]) / np.real(n_list[0])) * abs(t_tot) ** 2

def spectrum(d_ag):
    return np.array([tmm_T([n_air, na, n_glass], [0, d_ag, 0], l)
                     for na, l in zip(n_ag, wl)])

# bare glass reference (single interface air/glass)
T_glass = np.array([tmm_T([n_air, n_glass], [0, 0], l) for l in wl])

out = {"wl": wl}
plt.figure(figsize=(8, 5.2))
for d, c in [(8, "#2b7bba"), (15, "#c0392b")]:
    T = spectrum(d)
    out[f"T_abs_{d}nm"] = T
    out[f"T_rel_{d}nm"] = T / T_glass
    plt.plot(wl, 100 * T, c=c, lw=2, label=f"Ag {d} nm (절대 T)")
    plt.plot(wl, 100 * T / T_glass, c=c, lw=2, ls="--", label=f"Ag {d} nm (유리 대비 상대 T)")
    i550 = np.argmin(abs(wl - 550))
    print(f"Ag {d} nm @550nm: T_abs = {100*T[i550]:.1f}%, T_rel = {100*(T/T_glass)[i550]:.1f}%")
    print(f"  400-700 평균: T_abs = {100*T.mean():.1f}%, T_rel = {100*(T/T_glass).mean():.1f}%")

plt.xlabel("Wavelength (nm)"); plt.ylabel("Transmittance (%)")
plt.title("Ideal smooth Ag on glass — TMM, Johnson & Christy n,k")
plt.legend(); plt.grid(alpha=0.3); plt.ylim(0, 100)
png = os.path.join(BASE, "Ag_transmittance.png")
plt.tight_layout(); plt.savefig(png, dpi=160)
print("saved", png)

import csv
with open(os.path.join(BASE, "Ag_transmittance.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["wavelength_nm", "T_abs_8nm", "T_rel_8nm", "T_abs_15nm", "T_rel_15nm"])
    for i in range(len(wl)):
        w.writerow([f"{wl[i]:.0f}", f"{out['T_abs_8nm'][i]:.4f}", f"{out['T_rel_8nm'][i]:.4f}",
                    f"{out['T_abs_15nm'][i]:.4f}", f"{out['T_rel_15nm'][i]:.4f}"])
print("saved CSV")

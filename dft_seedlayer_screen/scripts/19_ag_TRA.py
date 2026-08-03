"""T / R / A of Ag films on glass, 400-700 nm.
(1) ideal smooth 8, 15 nm (TMM, J&C n,k)
(2) island-type film: Maxwell-Garnett effective medium (Ag in vacuum, fill f=0.5)
    with same equivalent mass thickness 8 nm -> physical layer 16 nm.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
plt.rcParams["font.family"] = "Malgun Gothic"
plt.rcParams["axes.unicode_minus"] = False

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

JC = np.array([
    [397.4, 0.05, 2.07], [413.3, 0.05, 2.21], [430.5, 0.04, 2.36],
    [450.9, 0.04, 2.66], [471.4, 0.05, 2.83], [495.9, 0.05, 3.09],
    [520.9, 0.05, 3.34], [548.6, 0.06, 3.59], [582.1, 0.05, 3.93],
    [616.8, 0.06, 4.15], [659.5, 0.05, 4.48], [704.5, 0.04, 4.84],
])
wl = np.linspace(400, 700, 301)
n_ag = np.interp(wl, JC[:, 0], JC[:, 1]) + 1j * np.interp(wl, JC[:, 0], JC[:, 2])
eps_ag = n_ag ** 2
n_air, n_glass = 1.0, 1.52

def tmm_TR(n_list, d_list, lam):
    M = np.eye(2, dtype=complex)
    for j in range(len(n_list) - 1):
        n1, n2 = n_list[j], n_list[j + 1]
        r = (n1 - n2) / (n1 + n2); t = 2 * n1 / (n1 + n2)
        I = np.array([[1, r], [r, 1]], dtype=complex) / t
        if j + 1 < len(n_list) - 1:
            delta = 2 * np.pi * n_list[j + 1] * d_list[j + 1] / lam
            P = np.array([[np.exp(-1j * delta), 0], [0, np.exp(1j * delta)]])
            M = M @ I @ P
        else:
            M = M @ I
    t_tot = 1 / M[0, 0]; r_tot = M[1, 0] / M[0, 0]
    T = (np.real(n_list[-1]) / np.real(n_list[0])) * abs(t_tot) ** 2
    R = abs(r_tot) ** 2
    return T, R

def spec(n_layer, d):
    T = np.empty_like(wl); R = np.empty_like(wl)
    for i, l in enumerate(wl):
        T[i], R[i] = tmm_TR([n_air, n_layer[i], n_glass], [0, d, 0], l)
    return T, R, 1 - T - R

# Maxwell-Garnett: Ag inclusions (f) in vacuum host
def mg_eps(eps_i, eps_h, f):
    return eps_h * (eps_i * (1 + 2 * f) + 2 * eps_h * (1 - f)) / \
           (eps_i * (1 - f) + eps_h * (2 + f))

f = 0.5
n_mg = np.sqrt(mg_eps(eps_ag, 1.0 + 0j, f))

cases = [
    ("Ag 8 nm 연속막", n_ag, 8.0, "#2b7bba"),
    ("Ag 15 nm 연속막", n_ag, 15.0, "#c0392b"),
    ("Ag 8 nm 상당 island막 (MG, f=0.5, 16 nm)", n_mg, 16.0, "#7d3c98"),
]

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
rows = []
for name, nl, d, c in cases:
    T, R, A = spec(nl, d)
    i550 = np.argmin(abs(wl - 550))
    print(f"{name}: @550nm T={100*T[i550]:.1f}% R={100*R[i550]:.1f}% A={100*A[i550]:.1f}% | "
          f"평균 T={100*T.mean():.1f}% R={100*R.mean():.1f}% A={100*A.mean():.1f}%")
    for ax, y, lab in zip(axes, [T, R, A], ["T", "R", "A"]):
        ax.plot(wl, 100 * y, c=c, lw=2, label=name)
    rows.append((name, T, R, A))

for ax, lab in zip(axes, ["투과 T", "반사 R", "흡수 A"]):
    ax.set_title(lab); ax.set_xlabel("Wavelength (nm)"); ax.grid(alpha=0.3)
axes[0].set_ylabel("(%)"); axes[0].set_ylim(0, 100)
axes[2].legend(fontsize=8.5, loc="upper right")
plt.suptitle("Ag 박막 T/R/A — 연속막(TMM) vs island막(Maxwell-Garnett 유효매질)", y=1.02)
png = os.path.join(BASE, "Ag_TRA.png")
plt.tight_layout(); plt.savefig(png, dpi=160, bbox_inches="tight")
print("saved", png)

import csv
with open(os.path.join(BASE, "Ag_TRA.csv"), "w", newline="", encoding="utf-8-sig") as fcsv:
    w = csv.writer(fcsv)
    hdr = ["wavelength_nm"]
    for name, *_ in cases:
        for q in ("T", "R", "A"):
            hdr.append(f"{q}_{name}")
    w.writerow(hdr)
    for i in range(len(wl)):
        row = [f"{wl[i]:.0f}"]
        for _, T, R, A in rows:
            row += [f"{T[i]:.4f}", f"{R[i]:.4f}", f"{A[i]:.4f}"]
        w.writerow(row)
print("saved CSV")

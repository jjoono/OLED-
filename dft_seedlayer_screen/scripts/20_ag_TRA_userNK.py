"""T/R/A with user's Ag n,k (Palik, McPeak from nk_JH_total.mat), 400-700 nm.
Continuous 8/15 nm + Maxwell-Garnett island films at f=0.3, 0.5 (mass-equiv 8 nm).
"""
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os, csv
plt.rcParams["font.family"] = "Malgun Gothic"
plt.rcParams["axes.unicode_minus"] = False

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
m = loadmat(r"C:\Users\Junho\Dropbox\Linkstation\Simulation\LosslessEML_single_distribution\nk_JH_total.mat")
mat = m["material"][0, 0]
wl = np.arange(400, 701, 1.0)
NK = {"Palik": np.asarray(mat["Ag_Palik"]).ravel(),
      "McPeak": np.asarray(mat["Ag_McPeak"]).ravel()}
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
    return (np.real(n_list[-1]) / np.real(n_list[0])) * abs(t_tot) ** 2, abs(r_tot) ** 2

def spec(n_layer, d):
    T = np.empty_like(wl); R = np.empty_like(wl)
    for i, l in enumerate(wl):
        T[i], R[i] = tmm_TR([n_air, n_layer[i], n_glass], [0, d, 0], l)
    return T, R, 1 - T - R

def mg_eps(eps_i, eps_h, f):
    return eps_h * (eps_i * (1 + 2 * f) + 2 * eps_h * (1 - f)) / \
           (eps_i * (1 - f) + eps_h * (2 + f))

fig, axes = plt.subplots(2, 3, figsize=(15, 9), sharex=True, sharey=True)
csv_rows = {"wavelength_nm": wl}
for row, (src, nk) in enumerate(NK.items()):
    eps = nk ** 2
    cases = [
        (f"연속 8 nm", nk, 8.0, "#2b7bba", "-"),
        (f"연속 15 nm", nk, 15.0, "#c0392b", "-"),
        (f"island f=0.5 (16 nm)", np.sqrt(mg_eps(eps, 1+0j, 0.5)), 16.0, "#7d3c98", "--"),
        (f"island f=0.3 (26.7 nm)", np.sqrt(mg_eps(eps, 1+0j, 0.3)), 26.7, "#1e8449", "--"),
    ]
    print(f"=== Ag {src} ===")
    for name, nl, d, c, ls in cases:
        T, R, A = spec(nl, d)
        tag = f"{src}_{name.split()[0]}_{d:g}nm" if "연속" in name else f"{src}_island_f{name.split('=')[1][:3]}"
        for q, y in zip(("T", "R", "A"), (T, R, A)):
            csv_rows[f"{q}_{tag}"] = y
        i550 = np.argmin(abs(wl - 550)); iA = int(np.argmax(A))
        print(f"  {name}: @550 T={100*T[i550]:.1f} R={100*R[i550]:.1f} A={100*A[i550]:.1f} | "
              f"A_max={100*A[iA]:.1f}% @{wl[iA]:.0f}nm | 평균 A={100*A.mean():.1f}%")
        for col, y in enumerate((T, R, A)):
            axes[row][col].plot(wl, 100 * y, c=c, ls=ls, lw=1.8,
                                label=name if col == 2 else None)
    for col, lab in enumerate(["투과 T", "반사 R", "흡수 A"]):
        axes[row][col].set_title(f"Ag {src} — {lab}", fontsize=11)
        axes[row][col].grid(alpha=0.3)
    axes[row][0].set_ylabel("(%)")
    axes[row][2].legend(fontsize=8)
for col in range(3):
    axes[1][col].set_xlabel("Wavelength (nm)")
plt.ylim(0, 100)
png = os.path.join(BASE, "Ag_TRA_userNK.png")
plt.tight_layout(); plt.savefig(png, dpi=150)
print("saved", png)

with open(os.path.join(BASE, "Ag_TRA_userNK.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    keys = list(csv_rows.keys())
    w.writerow(keys)
    for i in range(len(wl)):
        w.writerow([f"{csv_rows[k][i]:.4f}" if k != "wavelength_nm" else f"{wl[i]:.0f}" for k in keys])
print("saved CSV")

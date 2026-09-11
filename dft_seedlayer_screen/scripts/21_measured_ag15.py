"""T/R/A of measured 15 nm Ag films (#2, #3 on HATCN, ellipsometry nk from xlsx)
vs ideal Ag (McPeak). Air / Ag(15nm) / glass, 400-700 nm. T_rel = T / T(bare glass).
"""
import numpy as np, openpyxl, os, csv
from scipy.io import loadmat
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Malgun Gothic"
plt.rcParams["axes.unicode_minus"] = False

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
XLSX = r"C:\Users\Junho\.claude\uploads\298d98bb-4f41-4a12-af8c-f1c865056809\80fe6f6b-Ag15nm_MgAg_8nm_260703.xlsx"

wb = openpyxl.load_workbook(XLSX, data_only=True)
ws = wb["260703"]
lam, n2, n3, k2, k3 = [], [], [], [], []
for r in ws.iter_rows(min_row=5, values_only=True):
    if r[2] is None: continue
    lam.append(float(r[2]))
    n2.append(float(r[4])); n3.append(float(r[5]))
    k2.append(float(r[10])); k3.append(float(r[11]))
lam = np.array(lam)
sel = (lam >= 400) & (lam <= 700)
wl = lam[sel]
nk2 = np.array(n2)[sel] + 1j * np.array(k2)[sel]
nk3 = np.array(n3)[sel] + 1j * np.array(k3)[sel]
print(f"rows {len(wl)}, 550nm: #2 n+ik = {nk2[np.argmin(abs(wl-550))]:.3f}, "
      f"#3 = {nk3[np.argmin(abs(wl-550))]:.3f}")

# ideal reference
m = loadmat(r"C:\Users\Junho\Dropbox\Linkstation\Simulation\LosslessEML_single_distribution\nk_JH_total.mat")
wl_mc = np.arange(400, 701, 1.0)
nk_mc_raw = np.asarray(m["material"][0, 0]["Ag_McPeak"]).ravel()
nk_mc = np.interp(wl, wl_mc, nk_mc_raw.real) + 1j * np.interp(wl, wl_mc, nk_mc_raw.imag)

n_air, n_glass = 1.0, 1.52
def tmm_TR(n_list, d_list, l):
    M = np.eye(2, dtype=complex)
    for j in range(len(n_list) - 1):
        a, b = n_list[j], n_list[j + 1]
        r = (a - b) / (a + b); t = 2 * a / (a + b)
        I = np.array([[1, r], [r, 1]], dtype=complex) / t
        if j + 1 < len(n_list) - 1:
            d = 2 * np.pi * n_list[j + 1] * d_list[j + 1] / l
            P = np.array([[np.exp(-1j * d), 0], [0, np.exp(1j * d)]])
            M = M @ I @ P
        else:
            M = M @ I
    return (np.real(n_list[-1]) / np.real(n_list[0])) * abs(1 / M[0, 0]) ** 2, abs(M[1, 0] / M[0, 0]) ** 2

def spec(nk, d=15.0):
    T = np.empty_like(wl); R = np.empty_like(wl)
    for i, l in enumerate(wl):
        T[i], R[i] = tmm_TR([n_air, nk[i], n_glass], [0, d, 0], l)
    return T, R, 1 - T - R

T_glass = np.array([tmm_TR([n_air, n_glass], [0, 0], l)[0] for l in wl])

cases = [("#2 (측정 nk)", nk2, "#2b7bba"), ("#3 (측정 nk)", nk3, "#c0392b"),
         ("이상적 Ag (McPeak)", nk_mc, "#666666")]
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6), sharey=True)
out = {"wavelength_nm": wl}
for name, nk, c in cases:
    T, R, A = spec(nk)
    Trel = T / T_glass
    i550 = np.argmin(abs(wl - 550)); iA = int(np.argmax(A))
    print(f"{name}: @550 T={100*T[i550]:.1f}% T_rel={100*Trel[i550]:.1f}% R={100*R[i550]:.1f}% "
          f"A={100*A[i550]:.1f}% | A_max={100*A[iA]:.1f}%@{wl[iA]:.0f} | 평균 T={100*T.mean():.1f} "
          f"R={100*R.mean():.1f} A={100*A.mean():.1f}")
    key = name.split()[0]
    out[f"T_{key}"] = T; out[f"T_rel_{key}"] = Trel; out[f"R_{key}"] = R; out[f"A_{key}"] = A
    ls = "--" if "McPeak" in name else "-"
    ax[0].plot(wl, 100 * Trel, c=c, ls=ls, lw=2, label=name)
    ax[1].plot(wl, 100 * R, c=c, ls=ls, lw=2)
    ax[2].plot(wl, 100 * A, c=c, ls=ls, lw=2)
for a, t in zip(ax, ["상대투과 T/T(유리)", "반사 R", "흡수 A"]):
    a.set_title(f"Ag 15 nm — {t}"); a.set_xlabel("Wavelength (nm)"); a.grid(alpha=0.3)
ax[0].set_ylabel("(%)"); ax[0].set_ylim(0, 100); ax[0].legend(fontsize=9)
os.makedirs(os.path.join(BASE, "optics"), exist_ok=True)
png = os.path.join(BASE, "optics", "Ag15_measured_TRA.png")
plt.tight_layout(); plt.savefig(png, dpi=150)
print("saved", png)

with open(os.path.join(BASE, "Ag15_measured_TRA.csv"), "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f); keys = list(out.keys()); w.writerow(keys)
    for i in range(len(wl)):
        w.writerow([f"{out[k][i]:.4f}" if k != "wavelength_nm" else f"{wl[i]:.0f}" for k in keys])
print("saved CSV")

"""Fig 3(b) mock plot for the Ag reflector, from fig3b_Ag.csv (produced by fig3b_Ag.m). Recovered from the session plotting code."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
d = np.loadtxt("fig3b_Ag.csv", delimiter=","); nes = [1.8, 1.7, 1.6, 1.5]
col_b = {"Air":"#eb6834","Sub-confined":"#1baf7a","WG":"#4a3aa7","SPP (u > 1)":"#eda100","Abs":"#2a78d6"}
plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})
fig = plt.figure(figsize=(11, 6.6), facecolor="white")
gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.1], hspace=0.5, wspace=0.3, left=0.07, right=0.98, top=0.93, bottom=0.14)
for j, ne in enumerate([1.8, 1.6, 1.5]):
    ax = fig.add_subplot(gs[0, j]); r = d[np.isclose(d[:,0], ne)]; r = r[np.argsort(r[:,1])]; x = r[:,1]
    ax.stackplot(x, r[:,2], r[:,3], r[:,4], r[:,5], r[:,6], labels=list(col_b.keys()), colors=list(col_b.values()), alpha=0.85, lw=0)
    es = r[:,2]+r[:,3]; ax.plot(x, es, color="#0b0b0b", lw=1.4); ax.text(x[-1]-5, es[-1]+0.02, "η_sub", ha="right", va="bottom", fontsize=8)
    ax.set_xlim(x[0], x[-1]); ax.set_ylim(0, 1); ax.set_xlabel("d_ETL (nm)")
    ttl = {1.8: "isotropic ETL (n_e = n_o = 1.8): bound SPP", 1.6: "n_e,ETL = 1.6: n_SPP ≈ n_org, marginal", 1.5: "n_e,ETL = 1.5: leaky SPP → substrate"}[ne]
    ax.set_title(f"(b{j+1}) {ttl}", loc="left", fontsize=9)
    if ne == 1.5: ax.annotate("SPP power re-appears\nas substrate light", xy=(85, 0.7), xytext=(140, 0.45), fontsize=7.5, arrowprops=dict(arrowstyle="->", color="#0b0b0b", lw=0.8))
    if j == 0: ax.set_ylabel("Power ratio"); ax.legend(loc="lower right", fontsize=7, frameon=True, facecolor="white", framealpha=0.9)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
ax = fig.add_subplot(gs[1, :])
for ne, c in zip(nes, ["#0b0b0b", "#2a78d6", "#eb6834", "#1baf7a"]):
    r = d[np.isclose(d[:,0], ne)]; r = r[np.argsort(r[:,1])]
    ax.plot(r[:,1], r[:,2]+r[:,3], color=c, lw=2, marker="o", ms=3.5, label=("isotropic ETL, n_e = 1.8" if ne == 1.8 else f"n_e,ETL = {ne}"))
ax.axvspan(200, 400, color="#e6e6e3", lw=0); ax.text(300, 0.79, "driving-voltage penalty (undoped)", ha="center", fontsize=8, color="#52514e")
ax.axhline(0.9, color="#9a9a96", lw=0.8, ls=":"); ax.text(4, 0.905, "90 %", fontsize=8, color="#52514e")
ax.set_xlim(0, 400); ax.set_ylim(0.2, 1.0); ax.set_xlabel("d_ETL (nm)   (dipole at the centre of a 20-nm EML; distance to the metal = d_ETL + 10 nm)"); ax.set_ylabel("Substrate-delivered power η_sub")
ax.set_title("(b4)  Substrate-delivered power versus metal-adjacent layer thickness, Ag reflector", loc="left", fontsize=9.5)
ax.legend(loc="lower right", fontsize=8, frameon=False)
for s in ("top", "right"): ax.spines[s].set_visible(False)
ax.grid(axis="y", color="#e6e6e3", lw=0.6)
ins = ax.inset_axes([0.40, 0.22, 0.2, 0.5])
ne_grid = np.linspace(1.4, 1.8, 41); eo = 1.8**2
for lab, n, c in [("Ag", 0.044+3.819j, "#0b0b0b"), ("Al", 0.958+6.687j, "#9a9a96")]:
    em = n**2; ee = ne_grid**2; nspp = np.sqrt(em*ee*(em-eo)/(em**2 - eo*ee)).real
    ins.plot(ne_grid, nspp, color=c, lw=1.6, label=lab)
ins.axhline(1.8, color="#eb6834", lw=0.9, ls="--"); ins.text(1.41, 1.815, "n_org = n_sub = 1.8", fontsize=6.5, color="#eb6834")
ins.set_xlabel("n_e,ETL", fontsize=7, labelpad=1); ins.set_ylabel("n_SPP", fontsize=7); ins.tick_params(labelsize=6.5)
ins.set_title("SPP index: leaky when n_SPP < n_sub", fontsize=7); ins.legend(fontsize=6.5, frameon=False, loc="upper left")
fig.text(0.07, 0.058, "Generic model, 550 nm, isotropic dipole, PLQY = 1: Ag (McPeak n,k, 100 nm) / ETL (n_o = 1.8, n_e varied) / EML 1.8 (20 nm) / HTL 1.8 (50 nm) /", fontsize=7.5, color="#52514e")
fig.text(0.07, 0.035, "TCO 1.8+0.02i (50 nm, index-matched) / substrate 1.8. Planar_sweep22_preprint_v2.m with the branch-fixed TMF.", fontsize=7.5, color="#52514e")
fig.text(0.07, 0.012, "Inset: TM surface-plasmon dispersion for a uniaxial dielectric on a metal, k_SPP² = k0² ε_m ε_e (ε_m − ε_o)/(ε_m² − ε_o ε_e).", fontsize=7.5, color="#52514e")
fig.savefig("fig3b_Ag_mock.png", dpi=170); print("saved b")

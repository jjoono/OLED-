"""Fig 3(a) mock plot for the Ag reflector, from fig3a_Ag.csv (produced by fig3a_Ag.m). Recovered from the session plotting code."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
col = {"Air":"#eb6834","Sub-confined":"#1baf7a","WG":"#4a3aa7","Evanescent (SPP + TCO-guided)":"#eda100","Abs":"#2a78d6"}
plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})
a = np.loadtxt("fig3a_Ag.csv", delimiter=",")
fig, axs = plt.subplots(1, 3, figsize=(11, 3.8), facecolor="white")
plt.subplots_adjust(left=0.06, right=0.99, top=0.87, bottom=0.24, wspace=0.28)
for ax, ns in zip(axs, [1.5, 1.65, 1.8]):
    r = a[np.isclose(a[:,0], ns)]; r = r[np.argsort(r[:,1])]; x = r[:,1]
    ax.stackplot(x, r[:,2], r[:,3], r[:,4], r[:,5], r[:,6], labels=list(col.keys()), colors=list(col.values()), alpha=0.85, lw=0)
    es = r[:,2]+r[:,3]; ax.plot(x, es, color="#0b0b0b", lw=1.4); ax.text(x[-1]-6, es[-1]+0.02, "η_sub", ha="right", va="bottom", fontsize=8)
    ax.set_xlim(x[0], x[-1]); ax.set_ylim(0, 1); ax.set_xlabel("d_org (nm)   (dipole at d_org/2)")
    ax.set_title(f"(a{[1.5,1.65,1.8].index(ns)+1})  n_sub = {ns}", loc="left", fontsize=9.5)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axs[0].set_ylabel("Power ratio"); axs[2].legend(loc="center right", fontsize=6.8, frameon=True, facecolor="white", framealpha=0.9)
fig.text(0.06, 0.06, "Generic model, 550 nm, isotropic dipole, PLQY = 1: Ag (McPeak n,k, 100 nm) / organic n = 1.8 (d_org, dipole at centre) / ITO 2+0.02i (50 nm) / substrate n_sub.", fontsize=7.5, color="#52514e")
fig.text(0.06, 0.02, "For n_sub = 1.8 the u > 1 bin also contains ITO-guided light (n_ITO > n_sub). Computed with Planar_sweep22_preprint_v2.m and the branch-fixed TMF.", fontsize=7.5, color="#52514e")
fig.savefig("fig3a_Ag_mock.png", dpi=170); print("saved a")

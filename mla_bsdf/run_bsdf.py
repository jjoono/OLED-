"""Run the full incident-angle sweep and produce the BSDF maps + summary curves.

Outputs (written next to this script):
    mla_bsdf_result.npz   raw data
    bsdf_maps.png         BSDF_T(theta_i, theta_t) and BSDF_R(theta_i, theta_r)
    total_RT.png          integrated reflectance / transmittance vs incident angle
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from raytrace import sweep

HERE = os.path.dirname(os.path.abspath(__file__))


def main(N=300_000, seed=0):
    theta_i = np.arange(0.5, 90.0, 1.0)          # 1 deg bin centres, 0..90
    res = sweep(theta_i, N=N, n_glass=1.5, n_air=1.0, R=1.0, pitch=2.0,
                max_bounce=200, nbins=90, seed=seed)

    np.savez(os.path.join(HERE, "mla_bsdf_result.npz"), **{
        k: v for k, v in res.items() if k != "params"})

    _plot_maps(res)
    _plot_totals(res)

    # console summary
    print(" thi(deg)   T       R      T+R")
    for ti, T, R in zip(res["theta_i"], res["T_total"], res["R_total"]):
        if int(round(ti - 0.5)) % 10 == 0:
            print("  %5.1f  %6.3f  %6.3f  %6.3f" % (ti, T, R, T + R))
    print("max lost fraction: %.4g" % res["lost"].max())
    return res


def _plot_maps(res):
    ti = res["theta_i"]
    to = res["theta_out_centers"]
    ext = [to.min(), to.max(), ti.min(), ti.max()]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
    vmax = 5.0     # match the LightTools colour scale ( >5 % clipped )

    im0 = axes[0].imshow(res["bsdf_T"], origin="lower", aspect="auto",
                         extent=ext, cmap="jet", vmin=0, vmax=vmax)
    axes[0].set_title(r"$\mathrm{BSDF_T}(\theta_i,\theta_t)$")
    axes[0].set_xlabel("Transmitted angle (deg)")
    axes[0].set_ylabel("Incident angle (deg)")

    # reflected-angle axis reversed to match the LightTools layout
    im1 = axes[1].imshow(res["bsdf_R"][:, ::-1], origin="lower", aspect="auto",
                         extent=[to.max(), to.min(), ti.min(), ti.max()],
                         cmap="jet", vmin=0, vmax=vmax)
    axes[1].set_title(r"$\mathrm{BSDF_R}(\theta_i,\theta_r)$")
    axes[1].set_xlabel("Reflected angle (deg)")
    axes[1].set_ylabel("Incident angle (deg)")

    cb = fig.colorbar(im1, ax=axes, fraction=0.046, pad=0.04)
    cb.set_label("Power ratio (%)  [per 1$\\degree$ bin]")
    fig.suptitle("Hemispherical MLA (R=20 um, hex close-packed, n=1.5 on n=1.5 glass) "
                 "- geometric Monte-Carlo ray trace", fontsize=10)
    fig.savefig(os.path.join(HERE, "bsdf_maps.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)


def _plot_totals(res):
    ti = res["theta_i"]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(ti, res["T_total"] * 100, "-o", ms=3, label="Transmittance T")
    ax.plot(ti, res["R_total"] * 100, "-s", ms=3, label="Reflectance R")
    ax.axvline(np.rad2deg(np.arcsin(1 / 1.5)), color="grey", ls="--", lw=1,
               label=r"critical angle 41.8$\degree$")
    ax.set_xlabel("Incident angle in glass (deg)")
    ax.set_ylabel("Power fraction (%)")
    ax.set_title("Angle-integrated reflectance / transmittance")
    ax.set_xlim(0, 90)
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.savefig(os.path.join(HERE, "total_RT.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

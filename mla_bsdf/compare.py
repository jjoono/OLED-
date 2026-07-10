"""Compare convex vs concave MLA: overlay integrated R/T and export all CSVs."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(prefix):
    return np.load(os.path.join(HERE, prefix + "result.npz"))


def export_csv(prefix):
    d = _load(prefix)
    ti, to = d["theta_i"], d["theta_out_centers"]
    hdr = "incident_angle_deg\\output_angle_deg," + ",".join(f"{a:.1f}" for a in to)
    for tag, mat in (("T", d["bsdf_T"]), ("R", d["bsdf_R"])):
        np.savetxt(os.path.join(HERE, f"{prefix}bsdf_{tag}.csv"),
                   np.column_stack([ti, mat]), delimiter=",", header=hdr,
                   comments="", fmt=["%.1f"] + ["%.6g"] * mat.shape[1])
    np.savetxt(os.path.join(HERE, prefix + "totals.csv"),
               np.column_stack([ti, d["T_total"] * 100, d["R_total"] * 100,
                                d["lost"] * 100]),
               delimiter=",",
               header="incident_angle_deg,T_percent,R_percent,lost_percent",
               comments="", fmt="%.4f")


def compare_totals():
    cv = _load("convex_")
    cc = _load("concave_")
    ti = cv["theta_i"]
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    ax.plot(ti, cv["T_total"] * 100, "-", color="C0", label="convex  T")
    ax.plot(ti, cc["T_total"] * 100, "--", color="C0", label="concave T")
    ax.plot(ti, cv["R_total"] * 100, "-", color="C3", label="convex  R")
    ax.plot(ti, cc["R_total"] * 100, "--", color="C3", label="concave R")
    ax.axvline(np.rad2deg(np.arcsin(1 / 1.5)), color="grey", ls=":", lw=1)
    ax.text(42.5, 5, "41.8$\\degree$", color="grey")
    ax.set_xlabel("Incident angle in glass (deg)")
    ax.set_ylabel("Power fraction (%)")
    ax.set_title("Convex vs concave hemispherical MLA — integrated R / T")
    ax.set_xlim(0, 90)
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.3)
    ax.legend(ncol=2)
    fig.savefig(os.path.join(HERE, "compare_total_RT.png"), dpi=140,
                bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    for p in ("convex_", "concave_"):
        export_csv(p)
    compare_totals()
    print("wrote CSVs and compare_total_RT.png")

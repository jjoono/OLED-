"""Export the BSDF result (.npz) to CSV files for spreadsheet / LightTools comparison."""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, "mla_bsdf_result.npz"))
ti = d["theta_i"]
to = d["theta_out_centers"]


def save_matrix(name, mat):
    """First column = incident angle, header row = output-angle bin centres."""
    path = os.path.join(HERE, name)
    header = "incident_angle_deg\\output_angle_deg," + ",".join(f"{a:.1f}" for a in to)
    rows = np.column_stack([ti, mat])
    np.savetxt(path, rows, delimiter=",", header=header, comments="",
               fmt=["%.1f"] + ["%.6g"] * mat.shape[1])
    return path


# 1) angle-integrated totals
p1 = os.path.join(HERE, "mla_bsdf_totals.csv")
np.savetxt(p1, np.column_stack([ti, d["T_total"] * 100, d["R_total"] * 100,
                                d["lost"] * 100]),
           delimiter=",", header="incident_angle_deg,T_percent,R_percent,lost_percent",
           comments="", fmt="%.4f")

# 2) full BSDF maps (percent of incident power per 1 deg output bin)
p2 = save_matrix("mla_bsdf_T.csv", d["bsdf_T"])
p3 = save_matrix("mla_bsdf_R.csv", d["bsdf_R"])

for p in (p1, p2, p3):
    print(p)

# -*- coding: utf-8 -*-
"""afterAgMLA (2026-10-23) pixel d44_p3: JVL with the two available angle-resolved spectra.
Writes out/<case>/{*.xlsx,*.png} and out/afterAgMLA_summary.csv.  Run from this folder: python run_afterAgMLA.py"""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pixel_report as P
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "..", "data"); OUT = os.path.join(HERE, "out")
JVL = os.path.join(DATA, "2026-10-23_afterAgMLA_d44_p3_jvl.csv")
cases = [("d44_p3__angles_d44_p1", "2026-10-23_afterAgMLA_d44_p1_gon-spec-deg.csv", "2026-10-23 afterAgMLA d44_p3 (angles d44_p1)"),
         ("d44_p3__angles_d45_p3", "2026-10-23_afterAgMLA_d45_p3_gon-spec-deg.csv", "2026-10-23 afterAgMLA d44_p3 (angles d45_p3)")]
rows = []
for key, deg, name in cases:
    summ, xl, png, df = P.build(JVL, os.path.join(DATA, deg), os.path.join(OUT, key), name, tag=key)
    s = dict(zip(summ["항목"], summ["값"]))
    r = P.evaluate(JVL, os.path.join(DATA, deg)); c = r["corr"]
    rows.append(dict(case=key, EQE_max=s["최대 EQE"], EQE_100=s["EQE @ 100 cd/m²"], EQE_1000=s["EQE @ 1000 cd/m²"],
                     EQE_max_GUI_0to90=df["EQE 0~+90 only (%)"].max(), CE_max=s["최대 CE"], PE_max=s["최대 PE"],
                     PE_100=s["PE @ 100 cd/m²"], PE_1000=s["PE @ 1000 cd/m²"], V_100=s["전압 @ 100 cd/m²"], V_1000=s["전압 @ 1000 cd/m²"],
                     L_max=s["최대 휘도"], CIE=s["CIE (x, y), 0°"], peak_nm=s["피크 파장, 0°"], fwhm_nm=s["FWHM, 0°"],
                     e_corr=c["e"], e_pos=c["e_pos"], e_neg=c["e_neg"], v_corr=c["v"]))
t = pd.DataFrame(rows); t.to_csv(os.path.join(OUT, "afterAgMLA_summary.csv"), index=False, float_format="%.4f")
pd.set_option("display.width", 250); print(t.round(3).T.to_string())

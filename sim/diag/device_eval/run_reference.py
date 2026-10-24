# -*- coding: utf-8 -*-
"""Reference (flat glass) pixel 2026-10-22 glass d32_p5, evaluated exactly like run_afterAgMLA.py.
The uploaded JVL was re-saved with every row in double quotes; a de-quoted copy goes to out/ref_d32_p5/.
Run from this folder: python run_reference.py"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pixel_report as P
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "..", "data"); OUT = os.path.join(HERE, "out", "ref_d32_p5")
os.makedirs(OUT, exist_ok=True)
jvl = os.path.join(OUT, "2026-10-22_glass_d32_p5_jvl_dequoted.csv")
open(jvl, "w").write(open(os.path.join(DATA, "2026-10-22_glass_d32_p5_jvl.csv")).read().replace('"', ""))
summ, xl, png, df = P.build(jvl, os.path.join(DATA, "2026-10-22_glass_d32_p5_gon-spec-deg.csv"), OUT, "2026-10-22 glass d32_p5", tag="ref_d32_p5")
print(summ.to_string()); print("EQE max, GUI 0..+90: %.2f %%" % df["EQE 0~+90 only (%)"].max())

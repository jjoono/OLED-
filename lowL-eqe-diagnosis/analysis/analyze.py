# -*- coding: utf-8 -*-
"""Classify pixels and analyze low-luminance EQE anomaly."""
import numpy as np, pandas as pd, os

d = os.path.dirname(os.path.abspath(__file__))
df = pd.read_csv(os.path.join(d, "metrics.csv"))
n = len(df)
print(f"total jvl files: {n}")

# --- classification ---
dead = df["Lmax"] < 1.0                      # never emits
nir = (df["lum_per_mW"] < 20) & (df["Lmax"] >= 1) & (df["PDmax"] > 0.01)  # emits (PD sees light) but few cd per mW -> NIR-ish
dim = (df["Lmax"] >= 1) & (df["Lmax"] < 100) & ~nir   # never reaches 100 cd/m2
leaky = df["J_leak"] > 1e-2                  # >0.01 mA/cm2 at ~1V forward
good = ~dead & ~nir & ~dim
reliable = good & ~leaky & (df["Lmax"] >= 500)

print(f"dead (Lmax<1):            {dead.sum():6d} ({100*dead.mean():.1f}%)")
print(f"NIR-like (lum/mW<20):     {nir.sum():6d}")
print(f"dim (1<=Lmax<100):        {dim.sum():6d}")
print(f"leaky (J@1V>1e-2 mA/cm2): {leaky.sum():6d}")
print(f"good (Lmax>=100, vis):    {good.sum():6d}")
print(f"reliable (good, not leaky, Lmax>=500): {reliable.sum():6d}")

r = df[reliable].copy()

# noise floor stats
print("\n--- PD noise (off-state), luminance-equivalent [cd/m2] ---")
print(r["noise_L_equiv"].describe(percentiles=[.1,.25,.5,.75,.9,.99]))

# EQE ratio stats where all defined
for c in ["EQE_1","EQE_10","EQE_100","EQE_1000"]:
    r[c] = pd.to_numeric(r[c], errors="coerce")
m = r["EQE_100"] > 0.01
sub = r[m & r["EQE_1"].notna() & r["EQE_10"].notna()].copy()
sub["r1"] = sub["EQE_1"]/sub["EQE_100"]
sub["r10"] = sub["EQE_10"]/sub["EQE_100"]
sub["r1000"] = sub["EQE_1000"]/sub["EQE_100"]
print(f"\npixels with EQE defined at 1,10,100 cd/m2: {len(sub)}")
print("\n--- EQE(1)/EQE(100) ---")
print(sub["r1"].describe(percentiles=[.05,.25,.5,.75,.95]))
print("\n--- EQE(10)/EQE(100) ---")
print(sub["r10"].describe(percentiles=[.05,.25,.5,.75,.95]))
print("\n--- EQE(1000)/EQE(100) ---")
print(sub["r1000"].describe(percentiles=[.05,.25,.5,.75,.95]))

# how many are anomalous (>1.5x or <0.67x at 1 cd/m2)
print(f"\nEQE(1) > 1.5x EQE(100): {(sub['r1']>1.5).mean()*100:.1f}%")
print(f"EQE(1) < 0.67x EQE(100): {(sub['r1']<0.67).mean()*100:.1f}%")
print(f"EQE(1) within +-20% of EQE(100): {((sub['r1']>0.8)&(sub['r1']<1.2)).mean()*100:.1f}%")
print(f"EQE(10) within +-20% of EQE(100): {((sub['r10']>0.8)&(sub['r10']<1.2)).mean()*100:.1f}%")

# correlation of anomaly with noise floor
sub["absdev1"] = np.abs(np.log(sub["r1"]))
sub["noise_ratio_at_1"] = sub["noise_L_equiv"] / 1.0   # noise vs signal at 1 cd/m2
valid = sub[np.isfinite(sub["noise_ratio_at_1"]) & (sub["noise_ratio_at_1"]>0)]
print("\n--- correlation: |log EQE ratio at 1cd| vs log(noise floor) ---")
print(np.corrcoef(np.log10(valid["noise_ratio_at_1"]), valid["absdev1"])[0,1])

# bucket by noise floor
valid = valid.copy()
valid["nbucket"] = pd.cut(valid["noise_L_equiv"], [0, 0.03, 0.1, 0.3, 1, 3, 100])
print("\n--- anomaly at 1 cd/m2 by PD-noise bucket ---")
print(valid.groupby("nbucket", observed=True).agg(
    n=("r1","size"),
    median_r1=("r1","median"),
    p95_r1=("r1", lambda x: x.quantile(.95)),
    p05_r1=("r1", lambda x: x.quantile(.05)),
    frac_bad=("r1", lambda x: ((x>1.5)|(x<0.67)).mean())
))
print("\n--- anomaly at 10 cd/m2 by PD-noise bucket ---")
print(valid.groupby("nbucket", observed=True).agg(
    n=("r10","size"),
    median_r10=("r10","median"),
    frac_bad=("r10", lambda x: ((x>1.5)|(x<0.67)).mean())
))

# where does EQEmax occur?
lm = r[np.isfinite(r["L_at_EQEmax"]) & (r["EQEmax"]>0)]
print("\n--- luminance at which EQEmax occurs (reliable pixels) ---")
print(pd.cut(lm["L_at_EQEmax"], [0,1,10,100,1000,1e6]).value_counts().sort_index())

# leaky comparison: r1 for leaky-but-bright pixels
lk = df[good & leaky & (df["EQE_100"]>0.01)].copy()
lk["r1"] = lk["EQE_1"]/lk["EQE_100"]
print(f"\nleaky pixels r1 median: {lk['r1'].median():.2f} (n={lk['r1'].notna().sum()})")
print(f"leaky pixels: EQE(1)<0.67x EQE(100): {(lk['r1']<0.67).mean()*100:.1f}%  >1.5x: {(lk['r1']>1.5).mean()*100:.1f}%")

sub.to_csv(os.path.join(d, "reliable_ratios.csv"), index=False)
df.assign(cls=np.select(
    [dead, nir, dim, good & leaky, reliable],
    ["dead","nir","dim","leaky","reliable"], default="marginal"
)).to_csv(os.path.join(d, "classified.csv"), index=False)
print("\nwrote reliable_ratios.csv, classified.csv")

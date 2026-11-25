"""Quantitative SEM texture analysis for the HATCN / Ag thickness series.

Run locally on the original files (TIF is fine -- no need to convert):

    python 35_sem_texture.py  HATCN30.tif  HATCN30_Ag12.tif  HATCN30_Ag25.tif

Outputs sem_texture_results.csv + sem_texture.png with, per image:
  - correlation length from the 2D autocorrelation (characteristic feature size)
  - dominant length scale from the radially-averaged power spectral density
  - area fraction / count density / equivalent-diameter distribution of the
    bright protruding features (hillocks / abnormal grains)
  - gray-level RMS contrast (texture strength; NOT height -- use AFM for height)

Requires: numpy, scipy, pillow  (scikit-image optional, improves segmentation)
    pip install numpy scipy pillow scikit-image

Scale: the FEI/Thermo databar is auto-cropped. Set NM_PER_PX below, or leave it
None and pass --scalebar <px> <nm> if the automatic guess is wrong.
"""
import sys, os, csv
import numpy as np
from PIL import Image
from scipy import ndimage, signal

# ---- calibration -----------------------------------------------------------
# 50,000x on a 758-px-wide export -> ~3.4 nm/px. MEASURE the scale bar in your
# own export and set this; everything downstream scales linearly with it.
NM_PER_PX = 3.36
DATABAR_FRAC = 0.10          # bottom fraction of the frame occupied by the databar

def load(path):
    a = np.array(Image.open(path).convert("I")).astype(float)
    h = a.shape[0]
    a = a[: int(h * (1 - DATABAR_FRAC)), :]          # drop databar
    a = (a - a.min()) / (a.max() - a.min() + 1e-12)  # normalise 16/8-bit alike
    return a

def flatten(a, sigma=40):
    """Remove slow illumination/tilt gradient so texture stats are unbiased."""
    return a - ndimage.gaussian_filter(a, sigma)

def autocorr_length(a):
    """Correlation length: radius where the normalised autocorrelation drops to 1/e."""
    f = np.fft.fft2(a - a.mean())
    ac = np.real(np.fft.ifft2(f * np.conj(f)))
    ac = np.fft.fftshift(ac); ac /= ac.max()
    cy, cx = np.array(ac.shape) // 2
    yy, xx = np.indices(ac.shape)
    r = np.hypot(yy - cy, xx - cx).astype(int)
    prof = ndimage.mean(ac, labels=r, index=np.arange(r.max() + 1))
    below = np.where(prof < 1 / np.e)[0]
    return (below[0] if len(below) else np.nan) * NM_PER_PX, prof

def psd_peak(a):
    """Dominant lateral length scale from the radially averaged PSD."""
    f = np.fft.fftshift(np.fft.fft2(a - a.mean()))
    p = np.abs(f) ** 2
    cy, cx = np.array(p.shape) // 2
    yy, xx = np.indices(p.shape)
    r = np.hypot(yy - cy, xx - cx).astype(int)
    prof = ndimage.mean(p, labels=r, index=np.arange(1, r.max() + 1))
    k = np.arange(1, r.max() + 1)
    # weight by k to find the scale carrying the most variance per octave
    i = int(np.argmax(prof * k))
    lam = a.shape[0] / k[i] * NM_PER_PX
    return lam, k, prof

def bright_features(a):
    """Segment the bright protruding features above the fine-grained matrix."""
    fl = flatten(a, sigma=25)
    fl = ndimage.gaussian_filter(fl, 1.5)             # kill detector noise
    try:
        from skimage.filters import threshold_otsu
        thr = threshold_otsu(fl)
    except Exception:
        thr = fl.mean() + 1.0 * fl.std()
    mask = fl > thr
    mask = ndimage.binary_opening(mask, np.ones((3, 3)))
    lab, n = ndimage.label(mask)
    if n == 0:
        return mask, 0.0, 0.0, np.array([])
    sizes = np.bincount(lab.ravel())[1:]              # px per feature
    sizes = sizes[sizes >= 9]                         # drop noise specks
    d_eq = 2 * np.sqrt(sizes / np.pi) * NM_PER_PX     # equivalent diameter, nm
    area_frac = mask.sum() / mask.size
    area_um2 = mask.size * (NM_PER_PX / 1000.0) ** 2
    return mask, area_frac, len(sizes) / area_um2, d_eq

def main(paths):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = []
    fig, axes = plt.subplots(3, len(paths), figsize=(4.2 * len(paths), 11))
    if len(paths) == 1:
        axes = axes.reshape(3, 1)

    for j, p in enumerate(paths):
        a = load(p)
        lc, acprof = autocorr_length(flatten(a))
        lam, k, psd = psd_peak(flatten(a))
        mask, af, dens, d_eq = bright_features(a)
        rms = flatten(a, sigma=40).std()

        rows.append({
            "file": os.path.basename(p),
            "corr_length_nm": round(lc, 1),
            "psd_scale_nm": round(lam, 1),
            "bright_area_frac_%": round(100 * af, 2),
            "bright_count_per_um2": round(dens, 1),
            "bright_d_eq_median_nm": round(float(np.median(d_eq)), 1) if len(d_eq) else np.nan,
            "bright_d_eq_p90_nm": round(float(np.percentile(d_eq, 90)), 1) if len(d_eq) else np.nan,
            "graylevel_rms": round(float(rms), 4),
        })

        axes[0, j].imshow(a, cmap="gray"); axes[0, j].set_title(os.path.basename(p), fontsize=9)
        axes[1, j].imshow(mask, cmap="gray")
        axes[1, j].set_title(f"bright features: {100*af:.1f}% area, {dens:.0f}/um2", fontsize=8)
        if len(d_eq):
            axes[2, j].hist(d_eq, bins=30, color="#2b7bba")
            axes[2, j].set_xlabel("equivalent diameter (nm)")
            axes[2, j].set_title(f"median {np.median(d_eq):.0f} nm", fontsize=8)
        for ax in axes[:2, j]:
            ax.set_xticks([]); ax.set_yticks([])

    plt.tight_layout(); plt.savefig("sem_texture.png", dpi=140)
    with open("sem_texture_results.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    print(f"\n{'file':<28}{'corr_len':>10}{'psd_scale':>11}{'bright%':>9}"
          f"{'count/um2':>11}{'d50':>8}{'d90':>8}{'rms':>8}")
    for r in rows:
        print(f"{r['file']:<28}{r['corr_length_nm']:>10}{r['psd_scale_nm']:>11}"
              f"{r['bright_area_frac_%']:>9}{r['bright_count_per_um2']:>11}"
              f"{r['bright_d_eq_median_nm']:>8}{r['bright_d_eq_p90_nm']:>8}{r['graylevel_rms']:>8}")
    print(f"\n(NM_PER_PX = {NM_PER_PX}; check it against your scale bar)")
    print("wrote sem_texture_results.csv, sem_texture.png")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])

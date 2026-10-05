# -*- coding: utf-8 -*-
"""Compare backside light leakage between device photos in greyscale.

The DBR device's backside looks black to the eye, so a linear picture cannot
show how much light it lets through.  This converts every photo to *linear*
luminance and puts them all on one shared logarithmic grey scale, so a faint
glow becomes visible without any photo getting its own private gain.  It then
measures the central (lit) finger in each photo and reports how bright it is
relative to the reference photo.

    python3 figures/backside_grayscale.py \
        "w/o reflector=wo.jpg" "w/ DBR=dbr.jpg" \
        --meter "w/o reflector=998" --meter "w/ DBR=52" --out figures/backside

What the numbers mean, and when they mean nothing
-------------------------------------------------
* Pixels are sRGB-decoded to linear light, then combined into luminance
  Y = 0.2126 R + 0.7152 G + 0.0722 B -- the photopic weighting a luminance meter
  uses, so a photo ratio is directly comparable with an L_back ratio.
* The metric is the *peak surface brightness* of the lit finger (99.8th
  percentile of a 5x5-smoothed map), not an integral, so it does not depend on
  how far the camera was or how tightly the photo was cropped.
* A ratio between two photos is only physical if both were taken at the same
  exposure (shutter, ISO, aperture, manual white balance).  If the files carry
  EXIF, each photo is normalised by t*ISO/N^2 automatically.  Screenshots and
  crops pasted through PowerPoint usually carry none.
* Pass the luminance-meter readings with --meter and every photo ratio is
  checked against the meter ratio.  A mismatch beyond 1.5x is printed in red on
  the figure itself: it almost always means the exposures differed, and the
  photo numbers should not be quoted.
* A photo whose lit finger spans fewer than 32 code values is flagged as too
  dark or too quantised to measure.  Lossy exports (palette PNG, re-saved JPEG)
  are the usual cause; use the camera's original file.
"""
import argparse
import json
import os
import sys

import numpy as np
from PIL import Image

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

SMOOTH = 5            # box size (px) used to find the peak and the ROI
PEAK_PCT = 99.8       # "peak" = this percentile of the smoothed luminance
BLACK_PCT = 2.0       # camera black level = this percentile of the raw luminance
ROI_FRAC = 0.5        # ROI = pixels >= this fraction of the peak ...
ROI_RADIUS = 0.15     # ... within this fraction of the image width of the peak
BAND = 3              # line profile = mean of rows peak_y +- BAND
MIN_LEVELS = 32       # fewer code values than this in the finger -> unreliable
MISMATCH = 1.5        # photo/meter disagreement beyond this factor -> flagged
LOG_MIN = 1e-3        # bottom of the shared log grey scale (relative to reference)


def srgb_to_linear(c):
    c = c / 255.0
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def box_blur(a, k=SMOOTH):
    p = k // 2
    b = np.pad(a, p, mode="edge")
    c = np.pad(np.cumsum(np.cumsum(b, 0), 1), ((1, 0), (1, 0)))
    return (c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]) / float(k * k)


def exposure_of(img):
    """t * ISO / N^2 from EXIF, or None.  Brightness on the sensor scales with it."""
    try:
        ex = img.getexif().get_ifd(0x8769)
        t, n, iso = ex.get(33434), ex.get(33437), ex.get(34855)
        if isinstance(iso, (tuple, list)):
            iso = iso[0]
        if t and n and iso:
            return float(t) * float(iso) / float(n) ** 2
    except Exception:
        pass
    return None


def analyse(label, path, exposure=None):
    img = Image.open(path)
    exif_exp = exposure_of(img)
    rgb = np.asarray(img.convert("RGB")).astype(float)
    lin = srgb_to_linear(rgb)
    y = 0.2126 * lin[..., 0] + 0.7152 * lin[..., 1] + 0.0722 * lin[..., 2]

    exp = exposure if exposure is not None else exif_exp
    if exp:
        y = y / exp
    y = y - np.percentile(y, BLACK_PCT)

    sm = box_blur(y)
    peak = float(np.percentile(sm, PEAK_PCT))
    py, px = np.unravel_index(np.argmax(sm), sm.shape)
    yy, xx = np.mgrid[0:y.shape[0], 0:y.shape[1]]
    near = (yy - py) ** 2 + (xx - px) ** 2 <= (ROI_RADIUS * y.shape[1]) ** 2
    roi = near & (sm >= ROI_FRAC * peak) if peak > 0 else near & False
    code = rgb.max(2)
    levels = int(code[roi].max() - code[roi].min() + 1) if roi.any() else 0
    clipped = float((code[roi] >= 254).mean()) if roi.any() else 0.0

    lo, hi = max(0, py - BAND), min(y.shape[0], py + BAND + 1)
    return {
        "label": label, "path": path, "shape": list(y.shape),
        "exposure": exp, "exposure_source": "manual" if exposure else
                    ("exif" if exif_exp else None),
        "peak": peak, "roi_mean": float(y[roi].mean()) if roi.any() else 0.0,
        "roi_px": int(roi.sum()), "peak_xy": [int(px), int(py)],
        "code_levels": levels, "max_code": int(code.max()), "clipped": clipped,
        "Y": y, "profile": y[lo:hi].mean(0), "roi": roi,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("photos", nargs="+", help='"label=path", reference first')
    ap.add_argument("--meter", action="append", default=[],
                    help='"label=L_back" from the luminance meter, for the cross-check')
    ap.add_argument("--exposure", action="append", default=[],
                    help='"label=t*ISO/N^2" when the file has no EXIF')
    ap.add_argument("--px-per-mm", action="append", default=[],
                    help='"label=px/mm" to put the profile in millimetres')
    ap.add_argument("--out", default="backside_grayscale", help="output prefix")
    a = ap.parse_args(argv)

    kv = lambda items: {k: float(v) for k, v in (s.rsplit("=", 1) for s in items)}
    meter, expo, scale = kv(a.meter), kv(a.exposure), kv(a.px_per_mm)
    res = [analyse(lbl, p, expo.get(lbl)) for lbl, p in
           (s.split("=", 1) for s in a.photos)]
    ref = res[0]
    if ref["peak"] <= 0:
        sys.exit("reference photo has no measurable signal")

    notes = []
    for r in res:
        r["rel_peak"] = r["peak"] / ref["peak"]
        r["rel_roi"] = r["roi_mean"] / ref["roi_mean"] if ref["roi_mean"] else 0.0
        if r["code_levels"] < MIN_LEVELS:
            notes.append(("warn", "%s: lit finger spans only %d grey levels (max code %d) "
                          "-- too dark or quantised to measure" %
                          (r["label"], r["code_levels"], r["max_code"])))
        if r["clipped"] > 0.01:
            notes.append(("warn", "%s: %.0f%% of the finger is clipped at 255 "
                          "-- its brightness is underestimated" %
                          (r["label"], 100 * r["clipped"])))
        if r["label"] in meter and ref["label"] in meter and r is not ref:
            m = meter[r["label"]] / meter[ref["label"]]
            f = r["rel_peak"] / m if m else float("inf")
            r["meter_ratio"], r["photo_vs_meter"] = m, f
            bad = not (1 / MISMATCH <= f <= MISMATCH)
            notes.append(("bad" if bad else "ok",
                          "%s: photo %.3g vs meter %.3g (x%.2g)%s" %
                          (r["label"], r["rel_peak"], m, f,
                           " -- exposures likely differ" if bad else "")))
    if not any(r["exposure"] for r in res):
        notes.append(("warn", "no exposure data (no EXIF, no --exposure): ratios assume "
                      "identical camera settings"))

    # ------------------------------------------------------------ figure
    n = len(res)
    fig = plt.figure(figsize=(3.0 * n + 1.0, 6.2), dpi=200)
    gs = fig.add_gridspec(2, n + 1, width_ratios=[1] * n + [0.06],
                          height_ratios=[1.0, 0.9], hspace=0.42, wspace=0.12)
    norm = LogNorm(vmin=LOG_MIN, vmax=1.0)
    for i, r in enumerate(res):
        ax = fig.add_subplot(gs[0, i])
        im = ax.imshow(np.clip(r["Y"] / ref["peak"], LOG_MIN, None), cmap="gray",
                       norm=norm, interpolation="nearest")
        ax.contour(r["roi"], levels=[0.5], colors=["#e64b35"], linewidths=0.6)
        ax.set_title("%s\n%s of reference" % (r["label"], fmt_pct(r["rel_peak"])),
                     fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
    cax = fig.add_subplot(gs[0, n])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Relative luminance", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    axp = fig.add_subplot(gs[1, :n])
    for r in res:
        x = np.arange(r["profile"].size) - r["peak_xy"][0]
        unit = "px"
        if r["label"] in scale:
            x = x / scale[r["label"]]; unit = "mm"
        axp.plot(x, np.clip(r["profile"] / ref["peak"], LOG_MIN / 10, None),
                 lw=1.2, label=r["label"])
    axp.set_yscale("log")
    axp.set_ylim(LOG_MIN / 3, 2.0)
    axp.set_xlabel("Position from the lit finger (%s)" % unit, fontsize=8)
    axp.set_ylabel("Relative luminance", fontsize=8)
    axp.tick_params(labelsize=7)
    axp.legend(fontsize=7, frameon=False)
    axp.grid(True, which="major", lw=0.3, alpha=0.5)

    colour = {"ok": "#2e7d32", "bad": "#c62828", "warn": "#8a6d00"}
    for j, (kind, s) in enumerate(notes):
        fig.text(0.01, 0.005 + 0.028 * (len(notes) - 1 - j), s, fontsize=6.5,
                 color=colour[kind], ha="left", va="bottom")
    fig.subplots_adjust(bottom=0.08 + 0.028 * len(notes), top=0.92, left=0.09, right=0.95)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    fig.savefig(a.out + ".png")
    fig.savefig(a.out + ".svg")
    plt.close(fig)

    # one greyscale image per photo on the shared log scale, for the paper layout
    for r in res:
        v = np.log10(np.clip(r["Y"] / ref["peak"], LOG_MIN, 1.0))
        g = ((v - np.log10(LOG_MIN)) / -np.log10(LOG_MIN) * 255).round().astype("uint8")
        Image.fromarray(g).save("%s_%s.png" % (a.out, slug(r["label"])))

    keep = ("label", "path", "shape", "exposure", "exposure_source", "peak", "roi_mean",
            "roi_px", "peak_xy", "code_levels", "max_code", "clipped", "rel_peak",
            "rel_roi", "meter_ratio", "photo_vs_meter")
    with open(a.out + ".json", "w") as fh:
        json.dump({"photos": [{k: r[k] for k in keep if k in r} for r in res],
                   "notes": [s for _, s in notes]}, fh, indent=1)
    for r in res:
        print("%-22s peak %.4g  rel %-8s levels %3d  clipped %.0f%%" %
              (r["label"], r["peak"], fmt_pct(r["rel_peak"]), r["code_levels"],
               100 * r["clipped"]))
    for kind, s in notes:
        print("[%s] %s" % (kind, s))


def fmt_pct(v):
    return "%.2g%%" % (100 * v) if v < 0.1 else "%.0f%%" % (100 * v)


def slug(s):
    return "".join(c if c.isalnum() else "_" for c in s).strip("_").lower()


if __name__ == "__main__":
    main()

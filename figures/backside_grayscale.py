# -*- coding: utf-8 -*-
"""Compare backside light leakage between device photos in greyscale.

The DBR device's backside looks black to the eye, so a linear picture cannot
show how much light it lets through.  This converts every photo to *linear*
luminance and puts them all on one shared logarithmic grey scale, so a faint
glow becomes visible without any photo getting its own private gain, and it
measures how bright each device is relative to the first (reference) photo.

    python3 figures/backside_grayscale.py \
        "w/o reflector=wo.jpg" "w/ DBR=dbr.jpg" \
        --meter "w/o reflector=998" --meter "w/ DBR=52" --out figures/backside

Method
------
* Pixels are sRGB-decoded to linear light and combined into luminance
  Y = 0.2126 R + 0.7152 G + 0.0722 B -- the photopic weighting a luminance meter
  uses, so a photo ratio is comparable with an L_back ratio.
* The photos are of different substrates on the same holder, so they are
  aligned on the device itself: normalised cross-correlation of the
  log-luminance around the lit finger.  Every photo is then measured in the
  *same physical region* as the reference.
* The lit-finger region is where the reference is at least half its peak
  (5x5-smoothed), as one connected patch.  Its mean is the pixel's own backside
  brightness.  The same ratio is also taken over square windows of growing size,
  because the residual leakage of a reflector device need not come out where the
  pixel is: the curve shows where it goes.
* A luminance meter averages over its spot, so its L_back ratio should fall
  somewhere on that curve -- between the finger-only and the whole-device value.
  If it falls outside (beyond 1.5x), the photos were almost certainly taken at
  different exposures, and the figure says so in red.
* A ratio between photos is only physical at identical exposure (shutter, ISO,
  aperture, manual white balance).  EXIF exposure t*ISO/N^2 is divided out when
  present; uploads and slide exports usually strip it, so keep the originals.
"""
import argparse
import json
import os

import numpy as np
from PIL import Image

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter, FixedLocator, FuncFormatter

SMOOTH = 5             # box size (px) for the finger map
PEAK_SMOOTH = 9        # box size (px) used to locate the finger
ROI_FRAC = 0.5         # finger region = reference >= this fraction of its peak
LOCAL = 220            # half-size (px) of the patch the finger region may grow in
NCC_HALF = 120         # half-size (px) of the alignment template
NCC_SEARCH = 60        # search radius (px) for the alignment
BAND = 3               # line profile = mean of rows centre +- BAND
WINDOWS = (10, 20, 40, 80, 120, 160, 240, 320)   # half-sizes for the area curve
MIN_CODE = 32          # finger brighter than this code value is well resolved
MISMATCH = 1.5         # meter outside the curve by more than this -> flagged
LOG_MIN = 1e-3         # bottom of the shared log grey scale (relative to reference)
LINE_MARK = "#e64b35"  # the profile line drawn on the photos
LINE_COLOURS = ("#1a1a1a", "#2b62d9", "#d62728", "#3f9e4d")   # profile curves, in order


def srgb_to_linear(c):
    c = c / 255.0
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def box_blur(a, k):
    p = k // 2
    b = np.pad(a, p, mode="edge")
    c = np.pad(np.cumsum(np.cumsum(b, 0), 1), ((1, 0), (1, 0)))
    return (c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]) / float(k * k)


def exposure_of(img):
    """t * ISO / N^2 from EXIF, or None.  Signal on the sensor scales with it."""
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


def load(label, path, exposure=None):
    img = Image.open(path)
    exif_exp = exposure_of(img)
    rgb = np.asarray(img.convert("RGB")).astype(float)
    lin = srgb_to_linear(rgb)
    y = 0.2126 * lin[..., 0] + 0.7152 * lin[..., 1] + 0.0722 * lin[..., 2]
    exp = exposure if exposure is not None else exif_exp
    if exp:
        y = y / exp
    y = y - np.median(y)          # these scenes are mostly black: median = black level
    sm = box_blur(y, PEAK_SMOOTH)
    c = np.unravel_index(np.argmax(sm), sm.shape)
    return {"label": label, "path": path, "Y": y, "code": rgb.max(2),
            "gray8": np.asarray(img.convert("L")).astype(float),
            "exposure": exp, "exposure_source": "manual" if exposure else
                        ("exif" if exif_exp else None),
            "argmax": (int(c[0]), int(c[1]))}


def finger_region(y, centre):
    """Connected patch around `centre` where the smoothed map >= ROI_FRAC * its peak."""
    sm = box_blur(y, SMOOTH)
    y0, x0 = max(0, centre[0] - LOCAL), max(0, centre[1] - LOCAL)
    sub = sm[y0:centre[0] + LOCAL, x0:centre[1] + LOCAL]
    cy, cx = centre[0] - y0, centre[1] - x0
    mask = sub >= ROI_FRAC * sub[cy, cx]
    comp = np.zeros_like(mask)
    comp[cy, cx] = True
    while True:                    # grow the seed inside the mask until it stops
        d = comp.copy()
        d[1:] |= comp[:-1]; d[:-1] |= comp[1:]; d[:, 1:] |= comp[:, :-1]; d[:, :-1] |= comp[:, 1:]
        d &= mask
        if d.sum() == comp.sum():
            break
        comp = d
    full = np.zeros(y.shape, bool)
    full[y0:y0 + comp.shape[0], x0:x0 + comp.shape[1]] = comp
    return full


def align(ref, r):
    """Offset (dy, dx) that puts r's device onto the reference's, by NCC of log-luminance."""
    def logmap(p):
        return np.log10(np.clip(box_blur(p["Y"], SMOOTH), 0, None) + 3e-5)
    cy, cx = ref["centre"]
    t = logmap(ref)[cy - NCC_HALF:cy + NCC_HALF, cx - NCC_HALF:cx + NCC_HALF]
    t = (t - t.mean()) / (t.std() + 1e-12)
    L = logmap(r)
    H, W = L.shape

    def score(py, px):
        if py - NCC_HALF < 0 or px - NCC_HALF < 0 or py + NCC_HALF > H or px + NCC_HALF > W:
            return -2.0
        w = L[py - NCC_HALF:py + NCC_HALF, px - NCC_HALF:px + NCC_HALF]
        return float((t * (w - w.mean()) / (w.std() + 1e-12)).mean())

    gy, gx = r["argmax"]
    best = max(((score(gy + dy, gx + dx), gy + dy, gx + dx)
                for dy in range(-NCC_SEARCH, NCC_SEARCH + 1, 3)
                for dx in range(-NCC_SEARCH, NCC_SEARCH + 1, 3)))
    _, by, bx = best
    best = max(((score(by + dy, bx + dx), by + dy, bx + dx)
                for dy in range(-3, 4) for dx in range(-3, 4)))
    s, by, bx = best
    return (by - cy, bx - cx), s


def bilinear(img, ys, xs):
    y0 = np.clip(np.floor(ys).astype(int), 0, img.shape[0] - 2)
    x0 = np.clip(np.floor(xs).astype(int), 0, img.shape[1] - 2)
    fy, fx = ys - y0, xs - x0
    return (img[y0, x0] * (1 - fy) * (1 - fx) + img[y0 + 1, x0] * fy * (1 - fx)
            + img[y0, x0 + 1] * (1 - fy) * fx + img[y0 + 1, x0 + 1] * fy * fx)


def sample_line(img, centre, d0, d1, half, step=0.5):
    """Mean of `img` across a (2*half+1)-px band along the segment centre+d0 -> centre+d1.

    d0, d1 are (dx, dy) in px.  Positions are measured along the line from the foot
    of the perpendicular dropped from `centre`, so 0 is level with the lit finger.
    """
    p0 = np.array([centre[0] + d0[1], centre[1] + d0[0]], float)       # (y, x)
    p1 = np.array([centre[0] + d1[1], centre[1] + d1[0]], float)
    L = float(np.hypot(*(p1 - p0)))
    u = (p1 - p0) / L
    nrm = np.array([u[1], -u[0]])
    t = np.arange(0.0, L + 1e-9, step)
    acc = 0.0
    for w in range(-half, half + 1):
        acc = acc + bilinear(img, p0[0] + t * u[0] + w * nrm[0], p0[1] + t * u[1] + w * nrm[1])
    foot = -float(np.dot(np.array([d0[1], d0[0]], float), u))
    return t - foot, acc / (2 * half + 1)


def shifted(mask, off):
    return np.roll(np.roll(mask, off[0], 0), off[1], 1)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("photos", nargs="+", help='"label=path", reference first')
    ap.add_argument("--meter", action="append", default=[],
                    help='"label=L_back" from the luminance meter, for the cross-check')
    ap.add_argument("--exposure", action="append", default=[],
                    help='"label=t*ISO/N^2" when the file has no EXIF')
    ap.add_argument("--crop", type=int, default=150,
                    help="half-size (px) of the region shown in the figure")
    ap.add_argument("--line", default=None,
                    help='"dx0,dy0,dx1,dy1": a profile line in px relative to the lit '
                         'finger of the reference; every photo is sampled along the same '
                         'physical line')
    ap.add_argument("--line-half", type=int, default=2,
                    help="half-width (px) of the band averaged across the line")
    ap.add_argument("--out", default="backside_grayscale", help="output prefix")
    a = ap.parse_args(argv)

    kv = lambda items: {k: float(v) for k, v in (s.rsplit("=", 1) for s in items)}
    meter, expo = kv(a.meter), kv(a.exposure)
    res = [load(lbl, p, expo.get(lbl)) for lbl, p in (s.split("=", 1) for s in a.photos)]
    ref = res[0]
    ref["centre"] = ref["argmax"]
    ref["offset"], ref["ncc"] = (0, 0), 1.0
    roi = finger_region(ref["Y"], ref["centre"])
    yy, xx = np.mgrid[0:ref["Y"].shape[0], 0:ref["Y"].shape[1]]
    cy, cx = ref["centre"]
    wins = [(abs(yy - cy) <= h) & (abs(xx - cx) <= h) for h in WINDOWS]

    notes = []
    for r in res[1:]:
        r["offset"], r["ncc"] = align(ref, r)
        if r["ncc"] < 0.5:
            notes.append(("warn", "%s: alignment is weak (NCC %.2f) -- check the region "
                          "outlines in the figure" % (r["label"], r["ncc"])))
    for r in res:
        o = r["offset"]
        r["roi"] = shifted(roi, o)
        r["finger"] = float(r["Y"][r["roi"]].mean())
        r["windows"] = [float(r["Y"][shifted(w, o)].mean()) for w in wins]
        r["finger_code_max"] = int(r["code"][r["roi"]].max())
        r["clipped"] = float((r["code"][r["roi"]] >= 254).mean())
        r["centre"] = (cy + o[0], cx + o[1])
    for r in res:
        r["rel_finger"] = r["finger"] / ref["finger"]
        r["rel_windows"] = [v / w for v, w in zip(r["windows"], ref["windows"])]
        if r["finger_code_max"] < MIN_CODE:
            notes.append(("warn", "%s: finger reaches only code %d of 255 -- a mean over "
                          "%d px, so usable, but treat it as approximate (~+-20%%)" %
                          (r["label"], r["finger_code_max"], int(r["roi"].sum()))))
        if r["clipped"] > 0.01:
            notes.append(("warn", "%s: %.0f%% of the finger is clipped at 255 -- its "
                          "brightness is underestimated" % (r["label"], 100 * r["clipped"])))
        if r is not ref and r["label"] in meter and ref["label"] in meter:
            m = meter[r["label"]] / meter[ref["label"]]
            lo, hi = min(r["rel_windows"] + [r["rel_finger"]]), max(r["rel_windows"] + [r["rel_finger"]])
            r["meter_ratio"] = m
            ok = lo / MISMATCH <= m <= hi * MISMATCH
            notes.append(("ok" if ok else "bad",
                          "%s: meter %.3g lies %s the photo range %.3g (finger) - %.3g "
                          "(whole device)%s" % (r["label"], m,
                                                "within" if ok else "OUTSIDE", lo, hi,
                                                "" if ok else " -- exposures likely differ")))
    if not any(r["exposure"] for r in res):
        notes.append(("warn", "no exposure data (no EXIF, no --exposure): ratios assume "
                      "identical camera settings"))

    # ------------------------------------------------------------ line profile
    line = None
    if a.line:
        v = [float(t) for t in a.line.split(",")]
        line = ((v[0], v[1]), (v[2], v[3]))
        for r in res:
            r["line_pos"], r["line_Y"] = sample_line(r["Y"], r["centre"], *line, a.line_half)
            _, r["line_gray8"] = sample_line(r["gray8"], r["centre"], *line, a.line_half)
        top = float(ref["line_Y"].max())
        for r in res:
            r["line_rel"] = r["line_Y"] / top
            r["line_peak_rel"] = float(r["line_Y"].max()) / top

    # ------------------------------------------------------------ crops
    h = a.crop
    for r in res:
        y0, x0 = r["centre"][0] - h, r["centre"][1] - h
        r["crop"] = np.clip(r["Y"][max(0, y0):y0 + 2 * h, max(0, x0):x0 + 2 * h]
                            / ref["finger"], 0, None)
        r["crop_roi"] = r["roi"][max(0, y0):y0 + 2 * h, max(0, x0):x0 + 2 * h]
        r["profile"] = r["Y"][r["centre"][0] - BAND:r["centre"][0] + BAND + 1,
                              max(0, x0):x0 + 2 * h].mean(0) / ref["finger"]

    # ------------------------------------------------------------ analysis figure
    n = len(res)
    norm = LogNorm(vmin=LOG_MIN, vmax=1.0)
    fig = plt.figure(figsize=(3.0 * n + 4.6, 6.4), dpi=200)
    gs = fig.add_gridspec(2, n + 3, width_ratios=[1] * n + [0.06, 0.35, 1.25],
                          height_ratios=[1.0, 0.8], hspace=0.38, wspace=0.18)
    for i, r in enumerate(res):
        ax = fig.add_subplot(gs[0, i])
        im = ax.imshow(np.clip(r["crop"], LOG_MIN, None), cmap="gray", norm=norm,
                       interpolation="nearest")
        ax.contour(r["crop_roi"], levels=[0.5], colors=["#e64b35"], linewidths=0.6)
        if line:
            draw_line(ax, line, h, lw=0.9)
        ax.set_title("%s\nfinger: %s of reference" % (r["label"], pct(r["rel_finger"])),
                     fontsize=8.5)
        ax.set_xticks([]); ax.set_yticks([])
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[0, n]))
    cb.set_label("Relative luminance", fontsize=8); cb.ax.tick_params(labelsize=7)

    axw = fig.add_subplot(gs[0, n + 2])
    for r in res[1:]:
        axw.plot([2 * w + 1 for w in WINDOWS], r["rel_windows"], "o-", ms=3, lw=1.1,
                 label=r["label"])
        if "meter_ratio" in r:
            axw.axhline(r["meter_ratio"], ls="--", lw=0.9, color="0.4")
            axw.text(2 * WINDOWS[-1] + 1, r["meter_ratio"], " meter", fontsize=7,
                     va="center", color="0.3")
    axw.set_xscale("log"); axw.set_yscale("log")
    # explicit ticks: the ratios usually span less than a decade, which leaves a
    # log axis with no labelled tick at all, and its minor labels pile up
    lo_, hi_ = axw.get_ylim()
    yt = [v for v in (1e-3, 2e-3, 5e-3, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0) if lo_ <= v <= hi_]
    axw.yaxis.set_major_locator(FixedLocator(yt))
    axw.yaxis.set_major_formatter(FuncFormatter(lambda v, _: pct(v)))
    axw.xaxis.set_major_locator(FixedLocator([20, 50, 100, 200, 500]))
    axw.xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%d" % v))
    for ax_ in (axw.xaxis, axw.yaxis):
        ax_.set_minor_formatter(NullFormatter())
    axw.set_xlabel("Averaging window (px)", fontsize=8)
    axw.set_ylabel("Ratio to reference", fontsize=8)
    axw.tick_params(labelsize=7); axw.legend(fontsize=7, frameon=False)
    axw.grid(True, which="both", lw=0.3, alpha=0.5)

    axp = fig.add_subplot(gs[1, :n + 3])
    for r in res:
        x = np.arange(r["profile"].size) - h
        axp.plot(x, np.clip(r["profile"], LOG_MIN / 10, None), lw=1.1, label=r["label"])
    axp.set_yscale("log"); axp.set_ylim(LOG_MIN / 3, 2.0)
    axp.set_xlabel("Position across the fingers, from the lit one (px)", fontsize=8)
    axp.set_ylabel("Relative luminance", fontsize=8)
    axp.tick_params(labelsize=7); axp.legend(fontsize=7, frameon=False)
    axp.grid(True, which="major", lw=0.3, alpha=0.5)

    colour = {"ok": "#2e7d32", "bad": "#c62828", "warn": "#8a6d00"}
    for j, (kind, s) in enumerate(notes):
        fig.text(0.01, 0.005 + 0.026 * (len(notes) - 1 - j), s, fontsize=6.5,
                 color=colour[kind], ha="left", va="bottom")
    fig.subplots_adjust(bottom=0.09 + 0.026 * len(notes), top=0.92, left=0.06, right=0.97)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    fig.savefig(a.out + ".png"); fig.savefig(a.out + ".svg")
    plt.close(fig)

    # ------------------------------------------------------------ paper panel
    # just the greyscale crops on the shared scale and one colour bar, for panel c
    plt.rcParams["font.family"] = "sans-serif"     # first one installed wins
    plt.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
    fig = plt.figure(figsize=(1.55 * n + 0.45, 1.75), dpi=300)
    gs = fig.add_gridspec(1, n + 1, width_ratios=[1] * n + [0.07], wspace=0.06)
    for i, r in enumerate(res):
        ax = fig.add_subplot(gs[0, i])
        im = ax.imshow(np.clip(r["crop"], LOG_MIN, None), cmap="gray", norm=norm,
                       interpolation="bilinear")
        if line:
            draw_line(ax, line, h, lw=0.8)
        ax.set_title(r["label"], fontsize=7, pad=3)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[0, n]))
    cb.set_label("Relative luminance", fontsize=6.5, labelpad=2)
    cb.ax.tick_params(labelsize=6, length=2, pad=1)
    cb.outline.set_linewidth(0.5)
    fig.subplots_adjust(left=0.02, right=0.86, top=0.86, bottom=0.04)
    fig.savefig(a.out + "_panel.png"); fig.savefig(a.out + "_panel.svg")
    plt.close(fig)
    if line:
        fig, ax = plt.subplots(figsize=(2.6, 1.85), dpi=300)
        for r, col in zip(res, LINE_COLOURS):
            ax.plot(r["line_pos"], np.clip(r["line_rel"], LOG_MIN / 10, None), lw=1.0,
                    color=col, label=r["label"])
        ax.set_yscale("log")
        ax.set_ylim(LOG_MIN, 2.0)
        ax.set_xlim(r["line_pos"][0], r["line_pos"][-1])
        ax.set_xlabel("Position along the line (px)", fontsize=7, labelpad=2)
        ax.set_ylabel("Relative luminance", fontsize=7, labelpad=2)
        ax.tick_params(labelsize=6.5, length=2.5, width=0.6, pad=2, direction="in",
                       which="both", top=True, right=True)
        for sp in ax.spines.values():
            sp.set_linewidth(0.7)
        ax.legend(fontsize=6.5, frameon=False, loc="upper right", handlelength=1.4)
        fig.subplots_adjust(left=0.2, right=0.96, top=0.95, bottom=0.2)
        fig.savefig(a.out + "_line.png"); fig.savefig(a.out + "_line.svg")
        plt.close(fig)
        with open(a.out + "_line.csv", "w") as fh:
            cols = ["position_px"]
            for r in res:
                cols += ["%s rel. luminance" % r["label"], "%s gray (0-255)" % r["label"]]
            fh.write(",".join('"%s"' % c for c in cols) + "\n")
            for i in range(ref["line_pos"].size):
                row = ["%.2f" % ref["line_pos"][i]]
                for r in res:
                    row += ["%.6g" % r["line_rel"][i], "%.2f" % r["line_gray8"][i]]
                fh.write(",".join(row) + "\n")

    for r in res:
        v = np.log10(np.clip(r["crop"], LOG_MIN, 1.0))
        g = ((v - np.log10(LOG_MIN)) / -np.log10(LOG_MIN) * 255).round().astype("uint8")
        Image.fromarray(g).save("%s_%s.png" % (a.out, slug(r["label"])))

    keep = ("label", "path", "exposure", "exposure_source", "argmax", "centre", "offset",
            "ncc", "finger", "rel_finger", "rel_windows", "finger_code_max", "clipped",
            "meter_ratio", "line_peak_rel")
    with open(a.out + ".json", "w") as fh:
        json.dump({"finger_px": int(roi.sum()), "windows_px": [2 * w + 1 for w in WINDOWS],
                   "line": line, "line_half_px": a.line_half if line else None,
                   "photos": [{k: r[k] for k in keep if k in r} for r in res],
                   "notes": [s for _, s in notes]}, fh, indent=1, default=float)
    for r in res:
        print("%-20s offset %-12s NCC %.2f  finger %-7s whole-device %-7s" %
              (r["label"], tuple(int(v) for v in r["offset"]), r["ncc"],
               pct(r["rel_finger"]), pct(r["rel_windows"][-1])))
    if line:
        for r in res:
            print("%-20s line maximum %s of reference" % (r["label"], pct(r["line_peak_rel"])))
    for kind, s in notes:
        print("[%s] %s" % (kind, s))


def draw_line(ax, line, h, lw):
    """The profile line on a crop centred on the lit finger, as an arrow like the sketch."""
    (dx0, dy0), (dx1, dy1) = line
    ax.annotate("", xy=(h + dx1, h + dy1), xytext=(h + dx0, h + dy0),
                arrowprops=dict(arrowstyle="-|>", color=LINE_MARK, lw=lw,
                                mutation_scale=6, shrinkA=0, shrinkB=0))


def pct(v):
    return "%.2g%%" % (100 * v) if v < 0.1 else "%.0f%%" % (100 * v)


def slug(s):
    return "".join(c if c.isalnum() else "_" for c in s).strip("_").lower()


if __name__ == "__main__":
    main()

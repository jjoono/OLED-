# -*- coding: utf-8 -*-
"""Brighten a near-black device photo by an exact factor, for display.

    python3 figures/boost_photo.py dbr_back.jpg --gain 100 \
        --crop 854,602,1063,721 --size 401x228 --out dbr_back_x100.png

What is done, so a caption can say it in one line:

1. sRGB decode to linear light, multiply by --gain, re-encode.  "x100" means one
   hundred times the light, not a curve or contrast stretch; nothing is clipped
   unless the boosted pixel exceeds white, which the script reports.
2. Unless --raw: colour noise reduction.  A photo this dark lives in the bottom
   few code values, where JPEG colour subsampling leaves green/red/purple blocks
   that the gain makes obvious.  Luminance is kept from the photo (lightly
   smoothed, Gaussian sigma ~1.2 px, to soften the 8x8 JPEG steps); colour is
   the local average colour over ~9 px.
3. Crop and resize (bicubic) to match the photo it will sit next to.

Caption wording: "Back view of the DBR device, brightness x100 (linear gain;
colour noise reduced)."
"""
import argparse

import numpy as np
from PIL import Image


def srgb_decode(c):
    c = c / 255.0
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def srgb_encode(l):
    l = np.clip(l, 0, 1)
    return np.where(l <= 0.0031308, l * 12.92, 1.055 * l ** (1 / 2.4) - 0.055) * 255


def blur(a, k, passes=3):
    """Repeated box blur, ~Gaussian with sigma = k * sqrt(passes / 12)."""
    for _ in range(passes):
        p = k // 2
        b = np.pad(a, p, mode="edge")
        c = np.pad(np.cumsum(np.cumsum(b, 0), 1), ((1, 0), (1, 0)))
        a = (c[k:, k:] - c[:-k, k:] - c[k:, :-k] + c[:-k, :-k]) / float(k * k)
    return a


def boost(img, gain, raw=False):
    lin = srgb_decode(np.asarray(img.convert("RGB")).astype(float)) * gain
    clipped = float((lin.max(2) > 1.0).mean())
    if not raw:
        y = 0.2126 * lin[..., 0] + 0.7152 * lin[..., 1] + 0.0722 * lin[..., 2]
        ys = blur(y, 3, 2)
        yb = blur(y, 9) + 1e-12
        lin = np.stack([ys * blur(lin[..., i], 9) / yb for i in range(3)], -1)
    return Image.fromarray(srgb_encode(lin).round().astype("uint8")), clipped


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("photo")
    ap.add_argument("--gain", type=float, default=100.0)
    ap.add_argument("--crop", help="x0,y0,x1,y1 in the photo's pixels")
    ap.add_argument("--size", help="WxH of the output, e.g. 401x228")
    ap.add_argument("--raw", action="store_true", help="gain only, no noise reduction")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    img = Image.open(a.photo)
    if a.crop:
        img = img.crop(tuple(int(v) for v in a.crop.split(",")))
    out, clipped = boost(img, a.gain, a.raw)
    if a.size:
        out = out.resize(tuple(int(v) for v in a.size.lower().split("x")), Image.BICUBIC)
    out.save(a.out)
    print("wrote %s  gain x%g%s  clipped %.2f%% of pixels"
          % (a.out, a.gain, " (raw)" if a.raw else " (colour noise reduced)", 100 * clipped))


if __name__ == "__main__":
    main()

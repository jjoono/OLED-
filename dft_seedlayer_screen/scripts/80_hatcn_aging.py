"""Were the HATCN dendrites there at deposition, or did a week in air make them?

THE OBSERVATION. Bare HATCN films imaged 2026-08-12 (SEM metadata) show
dendritic islands, not continuous layers: Otsu coverage 30.5 % (2 nm nominal),
59.4 % (4 nm), 82.0 % (6 nm), continuous only at 30 nm. Dendritic branching is
the signature of diffusion-limited aggregation, i.e. crystalline growth.

WHY IT MATTERS. The two possible histories demand different fixes:
  as-deposited islands -> HATCN never wets the substrate; the seed layer is
      patchy from the start and needs more material or a different molecule.
  post-deposition dewetting -> HATCN deposits continuous and reorganises in
      storage; the fix is capping speed and handling, not the seed design.

THE TEST, using a control that already exists in the sample set. The samples
that received silver had their HATCN capped in the same vacuum run, so their
HATCN could not dewet afterwards. If the organic HAD been dendritic when the
silver landed, the silver would have to print that pattern:

  1. nucleation contrast is extreme -- our DFT puts Ag on the HATCN nitrile at
     1.03 eV and on the bare substrate near 0, so silver would decorate the
     islands preferentially;
  2. topography alone would do it -- 2 nm of material at 30 % coverage means
     islands ~6.7 nm tall, so a 3 nm silver film draped over them carries a
     6.7 nm height modulation at the 100-200 nm island scale, far stronger
     contrast than the 22-34 nm silver texture that dominates those images.

So: compute the radially-averaged PSD of each InLens frame and compare the
power in the DENDRITE band (100-250 nm, the island-spacing scale measured on
the bare films) against the SILVER-TEXTURE band (15-50 nm). Same detector,
same magnification, so the comparison is like-for-like.
"""
import os
import re

import numpy as np
from PIL import Image
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, ".."))
DENDRITE_BAND = (100.0, 250.0)      # nm, island spacing on the bare films
SILVER_BAND = (15.0, 50.0)          # nm, Ag grain/void texture (scripts/73)
FULL_BAND = (10.0, 400.0)

BARE = ["H2_1.tif", "H4_1.tif", "H6_1.tif", "H30_1.tif"]
CAPPED = ["H2_Ag3_1.tif", "H4_Ag3_1.tif", "H6_Ag_3_1.tif", "H30_Ag_3_1.tif",
          "H2_Ag5_1.tif", "H30_Ag5_1.tif", "H2_Ag25_1.tif"]


def band_powers(path):
    im = Image.open(path)
    blob = [v for v in im.tag_v2.values() if isinstance(v, str) and len(v) > 500][0]
    px = float(re.search(r"Image Pixel Size = ([0-9.]+)", blob).group(1))
    det = re.search(r"Signal A = (\w+)", blob).group(1)
    date = re.search(r"Date :([0-9]{1,2} [A-Za-z]{3} [0-9]{4})", blob)
    a = np.asarray(im.convert("L"), dtype=float)[:690, :]      # crop databar
    a = a - ndimage.gaussian_filter(a, 80)                     # flatten
    P = np.abs(np.fft.fftshift(np.fft.fft2(a - a.mean()))) ** 2
    cy, cx = np.array(P.shape) // 2
    yy, xx = np.ogrid[:P.shape[0], :P.shape[1]]
    rr = np.hypot(yy - cy, xx - cx).astype(int)
    k = np.arange(1, 340)
    pw = ndimage.mean(P, labels=rr, index=k) * k               # area weighted
    lam = a.shape[1] * px / k
    sel = lambda b: pw[(lam >= b[0]) & (lam <= b[1])].sum()
    tot = sel(FULL_BAND)
    return det, date.group(1) if date else "?", sel(DENDRITE_BAND) / tot, \
        sel(SILVER_BAND) / tot


def main():
    print(f"dendrite band {DENDRITE_BAND[0]:.0f}-{DENDRITE_BAND[1]:.0f} nm, "
          f"silver band {SILVER_BAND[0]:.0f}-{SILVER_BAND[1]:.0f} nm\n")
    print(f"{'file':<17}{'det':<8}{'imaged':<13}{'dendrite':>10}{'silver':>9}")
    for group, files in (("BARE HATCN (aged in air)", BARE),
                         ("HATCN CAPPED BY Ag", CAPPED)):
        print(f"\n-- {group}")
        for f in files:
            p = os.path.join(BASE, f)
            if not os.path.exists(p):
                print(f"{f:<17}  (missing)")
                continue
            det, date, d, s = band_powers(p)
            print(f"{f:<17}{det:<8}{date:<13}{100*d:>9.1f}%{100*s:>8.1f}%")


if __name__ == "__main__":
    main()


# RESULT (run 2026-08-16; SEM acquired 12 Aug 2026, ~1 week after deposition,
# samples stored in air throughout -- per the user):
#
#   bare HATCN   H2 54.4 %  H4 50.8 %  H6 39.3 %  H30 34.3 %   in the
#                                                              dendrite band
#   Ag-capped    Ag3 1.0 / 1.9 / 4.8 / 1.2 %,  Ag5 3.4 / 1.9 %,  Ag25 7.5 %
#
# A 10-50x difference. The silver films carry essentially NO power at the
# island-spacing scale, while the silver-texture band holds 61-72 % of theirs.
#
# VERDICT: the HATCN was continuous when the silver was deposited, and the
# dendrites grew during the week in air. The capped samples were protected;
# the bare ones were not.
#
# WHAT THIS RETRACTS. An earlier reading of scripts/73 -- that "Ag on H2 grew
# on 30 % HATCN plus 70 % bare substrate", and hence that crystallised HATCN
# is a worse silver surface than the bare substrate -- does NOT survive this
# test. Every Ag-covered sample sat on a continuous organic film. The
# scripts/73 trend (thinner HATCN -> smoother Ag, H2 best, H30 worst) still
# stands as an observation, but its MECHANISM is now open again, because the
# bare-film morphology we would use to explain it is an aged state, not the
# as-deposited one.
#
# WHAT SURVIVES, AND IS STRENGTHENED. HATCN's crystallisation risk (score 92
# in scripts/53) is now experimentally demonstrated: one week at room
# temperature in air is enough to take a continuous 2 nm film to 30 % coverage
# of dendrites. The mixed-seed strategy (HATCN:TPBi, scripts/68) keeps its
# rationale but changes its target: it is protection against reorganisation
# during storage and processing, not against crystallisation at deposition.
#
# PROTOCOL CONSEQUENCE. Every bare-HATCN reference measurement in this project
# -- SEM, and the ellipsometry that fitted 11-12 nm for 2-6 nm nominal -- was
# taken on aged films and describes the dendritic state, not the deposited
# one. Bare organic references must be imaged the same day or capped.
#
# STILL AN INFERENCE, NOT A PROOF. The decisive experiment is direct: deposit
# HATCN 2 nm, image the same spot on day 0, day 1 and day 7. AFM is the better
# instrument for it -- no vacuum, no charging, and the same area can be
# revisited.

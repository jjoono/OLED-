"""How large can HATCN's extinction coefficient be before the transparency claim fails?

WHY THIS EXISTS. scripts/32 reports "+0.00 %p absorption" for the HATCN seed. That
number is not a measurement and not really a result: the script sets k = 0 for
HATCN, so zero absorption is the input handed back. Quoting it as evidence that
HATCN is optically free is circular, and the manuscript draft did exactly that.

The honest question is not "is k zero" -- it is not, no organic film has k = 0 --
but "how large would k have to be before the seed's absorption matters?" That has
an answer that does not depend on knowing k, and it is the form the claim should
take until n,k are measured.

WHAT IS ACTUALLY KNOWN. HATCN's optical edge sits near 380-400 nm, so the visible
range is below the gap. Below-gap k in an amorphous organic film is not zero: an
Urbach tail carries it, typically 1e-3 to 1e-2, rising steeply at the edge. On top
of that, a HIL in contact with a donor (or with Ag) forms interfacial
charge-transfer states that absorb in the visible -- that is the mechanism that
makes HATCN work as a HIL. The CT absorption scales with INTERFACE AREA, not with
layer thickness, so it is not captured by a bulk k at all and is treated separately
in the discussion, not here.

This script sweeps k over four decades and reports photopic-weighted absorption of
the HATCN layer in the real stack, so the claim can be stated as a bound.
"""
import os, sys
import numpy as np
import importlib.util

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.abspath(os.path.join(HERE, ".."))

# Reuse the verified TMM from scripts/32 rather than reimplementing it. That
# implementation was cross-checked against an independent scattering-matrix code
# to 1.9e-15 (scripts/33), so the optics here are not the uncertain part.
_s = importlib.util.spec_from_file_location(
    "s32", os.path.join(HERE, "32_seed_absorption.py"))
S = importlib.util.module_from_spec(_s)
# scripts/32 runs its whole grid at import time and prints it. Swallow that so the
# sensitivity table is the only thing on stdout.
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    _s.loader.exec_module(S)

K_GRID = [0.0, 1e-3, 3e-3, 1e-2, 3e-2, 0.05, 0.1, 0.2, 0.5]
HATCN_T = [3.0, 5.0, 10.0, 30.0]      # nm; 30 nm is the film already deposited
AG_T = 8.0

# Threshold at which the "optically free" claim stops being defensible. 0.5 %p is
# roughly the reproducibility of a transmittance measurement on these stacks, so
# below it the seed's loss cannot be distinguished from run-to-run scatter.
CLAIM_LIMIT = 0.5


def absorption(k, d_hatcn, d_ag, superstrate="air"):
    """Photopic-weighted absorption of the HATCN layer alone, in the real stack."""
    S.NK["HATCN"] = np.full_like(S.wl, 1.95 + 1j * k, dtype=complex)
    r = S.spectrum([("HATCN", d_hatcn), ("Ag", d_ag)], superstrate)
    # spectrum() stacks as ["glass"] + stack materials + [superstrate], so the
    # seed is column 1, not 0. Taking column 0 would report the (lossless) glass
    # and would show 0.000 for every k -- the same circularity being fixed here.
    i_seed = r["mats"].index("HATCN")
    return 100.0 * S.photopic(r["Alay"][:, i_seed])


def main():
    print(f"Photopic-weighted absorption of the HATCN layer (%), Ag = {AG_T} nm\n")
    print(f"{'k':>8}" + "".join(f"{f'{d:.0f} nm':>11}" for d in HATCN_T))
    print("-" * (8 + 11 * len(HATCN_T)))
    table = {}
    for k in K_GRID:
        row = []
        for d in HATCN_T:
            try:
                row.append(absorption(k, d, AG_T))
            except Exception as exc:
                print(f"  failed at k={k}, d={d}: {type(exc).__name__}: {exc}")
                row.append(float("nan"))
        table[k] = row
        print(f"{k:>8.3f}" + "".join(f"{v:>11.3f}" for v in row))

    print("\n" + "=" * 70)
    print(f"Largest k for which each thickness stays under {CLAIM_LIMIT} %p:")
    for i, d in enumerate(HATCN_T):
        ok = [k for k in K_GRID if table[k][i] < CLAIM_LIMIT]
        best = max(ok) if ok else None
        print(f"  {d:>5.0f} nm   k <= {best if best is not None else 'none'}")

    print("\nHOW TO STATE THIS IN THE PAPER")
    print("  Not: 'HATCN contributes 0.00 %p' -- that is the k = 0 input returned.")
    print("  Instead: 'for any k below <bound>, the seed contributes less than")
    print("  0.5 %p, which is under the measurement reproducibility; the")
    print("  transparency claim therefore does not rest on the precise value.'")
    print("\n  And measure it: a 30 nm HATCN film -- already deposited for the SEM")
    print("  series -- gives n,k directly by ellipsometry, or alpha directly by")
    print("  transmission on quartz. That converts the bound into a number.")
    print("\n  Separately: interfacial CT absorption (HATCN/donor, HATCN/Ag) is NOT")
    print("  described by a bulk k and is not bounded by this sweep. It scales with")
    print("  interface area, so it does not vanish as the layer is made thinner.")


if __name__ == "__main__":
    main()

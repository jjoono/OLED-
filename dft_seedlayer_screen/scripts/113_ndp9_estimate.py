"""Predict E_b for NDP-9 from the correlation the screen already established.

The direct calculation is a counterpoise-corrected rigid scan on a 49-atom
complex, three SCF per scan point. It does not fit on four cores, and it belongs
on the workstation (scripts/112). Meanwhile the screen contains a usable
regularity: for a nitrile binder, E_b tracks how strongly the attached ring
pulls charge off the silver, and electron affinity measures that.

Five nitrile binders have both quantities on record. The fit is over five
points, so it is an estimate with an honest error bar, not a result -- and the
question it has to answer is narrow enough that an estimate can answer it. NDP-9
is already 1.35x behind HATCN on site density, which by Venables needs 27 meV of
extra E_d, about 0.10 eV of extra E_b, to overturn. The fit only has to say
whether NDP-9 clears 1.13 eV.

EA has to be computed at the SAME level as the fitted points or the new one is
not on their axis. Those came from PBE-D3BJ/def2-SVP vertical differences, so
this does too. A first attempt used GFN2-xTB for speed and returned 9.3 eV for
NDP-9 against 2.7 for HATCN -- an impossible number that survived only because
it was never compared against the method that produced the rest of the column.
"""
import json, os, sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathgeom import RUNS, STRUCT, read_xyz

H2EV = 27.211386
KT, CHI = 0.025852, 0.284
ED_FROM_EB = 0.28

# nitrile binders with both E_b (PBE-D3BJ/def2-SVP, CP-corrected) and EA on file
FIT = [
    ("DMABN",  -1.042, 0.287),
    ("PhCN",   -0.462, 0.292),
    ("oDCNB",   0.549, 0.379),     # E_b from the ortho isomer scan
    ("pDCNB",   0.742, 0.379),
    ("HATCN",   2.711, 1.029),
]


def dft_ea(fn):
    """Vertical electron affinity, PBE-D3BJ/def2-SVP -- the level the fit uses."""
    import psi4
    psi4.set_memory(f"{os.environ.get('PSI4_MEM_GB', '9')} GB")
    try:
        psi4.set_num_threads(len(os.sched_getaffinity(0)))
    except AttributeError:
        psi4.set_num_threads(os.cpu_count() or 4)
    psi4.core.set_output_file(os.path.join(RUNS, "psi4_ndp9_ea.out"), False)
    syms, xyz = read_xyz(os.path.join(STRUCT, fn))
    keep = [i for i, q in enumerate(syms) if q != "Ag"]
    syms, xyz = [syms[i] for i in keep], xyz[keep]

    def sp(charge, mult):
        st = f"{charge} {mult}\n"
        for a, c in zip(syms, xyz):
            st += f"{a} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}\n"
        st += "symmetry c1\nno_reorient\nno_com\n"
        for opts in (dict(BASE, maxiter=150, guess="sad", level_shift=1.0,
                          level_shift_cutoff=1e-3, damping_percentage=15.0),
                     dict(BASE, maxiter=200, guess="sad", level_shift=2.0,
                          level_shift_cutoff=1e-2, damping_percentage=40.0),
                     dict(BASE, maxiter=200, guess="sad", soscf=True,
                          soscf_start_convergence=1e-2)):
            try:
                psi4.set_options(opts)
                e = psi4.energy("pbe-d3bj", molecule=psi4.geometry(st))
                psi4.core.clean()
                return e
            except Exception:
                psi4.core.clean()
        return None

    print("    neutral...", flush=True)
    e0 = sp(0, 1)
    print(f"    neutral {e0}", flush=True)
    print("    anion...", flush=True)
    e1 = sp(-1, 2)
    print(f"    anion   {e1}", flush=True)
    if None in (e0, e1):
        return None
    return (e0 - e1) * H2EV


BASE = {"basis": "def2-svp", "scf_type": "df", "reference": "uks"}


def main():
    ea = np.array([p[1] for p in FIT])
    eb = np.array([p[2] for p in FIT])
    A = np.vstack([np.ones_like(ea), ea]).T
    coef, *_ = np.linalg.lstsq(A, eb, rcond=None)
    pred = A @ coef
    resid = eb - pred
    sd = float(np.sqrt(np.sum(resid**2) / (len(ea) - 2)))
    r2 = 1 - np.sum(resid**2) / np.sum((eb - eb.mean())**2)

    print("nitrile binders on record\n")
    print(f"{'':<8} {'EA (eV)':>9} {'E_b':>7} {'fit':>7} {'resid':>7}")
    for (n, a, b), p, r in zip(FIT, pred, resid):
        print(f"{n:<8} {a:>9.3f} {b:>7.3f} {p:>7.3f} {r:>+7.3f}")
    print(f"\nE_b = {coef[0]:.3f} + {coef[1]:.3f} x EA     R2 = {r2:.3f}, "
          f"residual sd = {sd:.3f} eV")

    print("\ncomputing NDP-9 electron affinity...", flush=True)
    a9 = dft_ea("NDP9.xyz")
    if a9 is None:
        print("  SCF failed")
        return
    e9 = coef[0] + coef[1] * a9
    print(f"  EA(NDP-9) = {a9:.3f} eV   (HATCN {FIT[-1][1]:.3f})")
    print(f"  predicted E_b = {e9:.3f} +/- {sd:.3f} eV\n")

    # what it would take to beat HATCN, from the site-density gap
    gap = 13.7 / 10.2                       # sites per nm2, script 111
    d_ed = KT / CHI * np.log(gap)
    need_eb = (0.286 + d_ed) / ED_FROM_EB
    print(f"NDP-9 is behind on site density by {gap:.2f}x, which Venables converts")
    print(f"to {d_ed*1000:.0f} meV of E_d, about {need_eb:.2f} eV of E_b. HATCN is at 1.029,")
    print(f"so NDP-9 has to reach roughly {need_eb:.2f} eV to draw level.\n")

    z = (need_eb - e9) / sd
    print(f"  predicted {e9:.3f}, needs {need_eb:.2f} -- short by "
          f"{need_eb-e9:+.3f} eV, {z:.1f} sd")
    print("\nCAVEATS, and they are not small. Five points, two of which (oDCNB and")
    print("pDCNB) share an E_b. The relation is being used above the range it was")
    print("fitted on if NDP-9's EA exceeds HATCN's, and extrapolation is where a")
    print("five-point line is least trustworthy. A vertical EA at def2-SVP is also")
    print("not an absolute electron affinity; it is consistent with the other")
    print("points here, which is what the fit needs, and nothing more.")
    print("Treat this as a prediction to be checked by scripts/112, not a result.")


if __name__ == "__main__":
    main()

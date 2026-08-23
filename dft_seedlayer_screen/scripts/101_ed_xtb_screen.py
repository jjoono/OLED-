"""E_d for the whole candidate set at GFN2-xTB, anchored to the DFT points.

The DFT protocol in script 97 is correct and unaffordable here: one SCF on a
31-atom Ag complex runs about ten minutes on four cores, and the campaign is
twenty-six systems times fifteen single points. So the barrier is screened with
GFN2-xTB, which is roughly three orders of magnitude cheaper, and the screen is
made trustworthy the only way a semiempirical screen can be -- by reproducing
the DFT numbers on the systems where both exist.

An earlier xTB attempt in this project returned 250 eV for B3PyMPM, which is why
xTB was set aside. That was a runaway geometry, not a method failure: nothing
constrained the adatom, and it was reported without a sanity check. Two things
are different here. The path geometry is identical to the DFT one -- same
structures, same drag, same z-scan -- so xTB is only ever asked for an energy
difference along a fixed path. And every result is range-checked before it is
kept: a barrier outside 0 to 3 eV is recorded as a failure, not as a number.

The anchor set decides whether the screen may be used at all. If xTB reproduces
the DFT ordering and rough magnitude on the systems that have both, the rest of
the column can be read as a screen; if it does not, this file reports that and
the DFT campaign is the only route.
"""
import json, os, shutil, subprocess, sys, tempfile
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import importlib
d97 = importlib.import_module("97_diffusion_barriers")

BASE = d97.BASE
STRUCT, RUNS = d97.STRUCT, d97.RUNS
OUT = os.path.join(RUNS, "diffusion_barriers_xtb.json")
XTB = "/root/miniforge3/envs/qc/bin/xtb"
H2EV = d97.H2EV
SANE = (-0.2, 3.0)          # eV; anything outside is a failed point, not a value

# DFT values already in hand, from scripts 24 and 24b
DFT_ANCHOR = {"HATCN": 0.2857, "AlOx": 0.1696, "cleanAl": 0.0, "TAPC": -0.0072}
ANCHOR_FILE = {"HATCN": "HATCN_Ag_CN.xyz", "AlOx": "Al4O6_Ag.xyz"}


def xtb_energy(syms, xyz, chrg=0, uhf=1):
    """GFN2-xTB single-point energy in hartree, or None if it does not converge."""
    d = tempfile.mkdtemp(prefix="xtb_")
    try:
        with open(os.path.join(d, "in.xyz"), "w") as f:
            f.write(f"{len(syms)}\n\n")
            for a, c in zip(syms, xyz):
                f.write(f"{a} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}\n")
        r = subprocess.run([XTB, "in.xyz", "--gfn", "2", "--chrg", str(chrg),
                            "--uhf", str(uhf), "--sp", "--norestart"],
                           cwd=d, capture_output=True, text=True, timeout=300)
        for line in r.stdout.splitlines():
            if "TOTAL ENERGY" in line:
                return float(line.split()[3])
        return None
    except (subprocess.TimeoutExpired, ValueError, IndexError):
        return None
    finally:
        shutil.rmtree(d, ignore_errors=True)


def barrier(tag, fn, rule):
    syms, xyz = d97.read_xyz(os.path.join(STRUCT, fn))
    sub_s, sub_x, ag, anchor, nrm = d97.geometry(syms, xyz)
    dest, cls = d97.destination(sub_s, sub_x, ag, anchor, rule)
    span = float(np.linalg.norm(dest - ag))
    npath = int(np.clip(round(span / d97.SPACING) + 1, d97.NPATH_MIN, d97.NPATH_MAX))

    E = []
    for t in np.linspace(0.0, 1.0, npath):
        pos = d97.place(sub_x, (1 - t) * ag + t * dest, nrm)
        best = None
        for dz in d97.ZSCAN:
            e = xtb_energy(sub_s + ["Ag"], np.vstack([sub_x, pos + dz * nrm]))
            if e is not None and (best is None or e < best):
                best = e
        E.append(best)

    ok = [e for e in E if e is not None]
    if E[0] is None or len(ok) < 3:
        return {"E_d_eV": None, "class": cls, "note": "too few converged points"}
    ed = (max(ok) - E[0]) * H2EV
    if not (SANE[0] <= ed <= SANE[1]):
        return {"E_d_eV": None, "class": cls,
                "note": f"out of range ({ed:.2f} eV) -- rejected"}
    return {"E_d_eV": round(ed, 4), "class": cls, "anchor": sub_s[anchor],
            "n_atoms": len(sub_s), "path_A": round(span, 3), "n_points": npath}


def main():
    res = {}
    for tag, fn, rule, _ in d97.CANDIDATES:
        p = os.path.join(STRUCT, fn)
        if not os.path.exists(p):
            print(f"{tag:<12} missing {fn}", flush=True)
            continue
        r = barrier(tag, fn, rule)
        res[tag] = r
        v = r["E_d_eV"]
        print(f"{tag:<12} {'--' if v is None else f'{v:6.3f} eV'}  {r['class']:<10}"
              f"{'  ' + r.get('note', '')}", flush=True)
        json.dump(res, open(OUT, "w"), indent=1)

    print("\nanchor check -- xTB against the DFT values already in hand")
    print(f"{'':<12} {'DFT':>8} {'xTB':>8} {'diff':>8}")
    pairs = []
    for tag, dft in DFT_ANCHOR.items():
        key = tag if tag in res else None
        if key is None and tag == "AlOx":
            key = "Al4O6"
        if key is None or res.get(key, {}).get("E_d_eV") is None:
            print(f"{tag:<12} {dft:>8.3f} {'--':>8}")
            continue
        x = res[key]["E_d_eV"]
        pairs.append((dft, x))
        print(f"{tag:<12} {dft:>8.3f} {x:>8.3f} {x-dft:>+8.3f}")
    if len(pairs) >= 2:
        a = np.array(pairs)
        mae = float(np.mean(np.abs(a[:, 1] - a[:, 0])))
        print(f"\nMAE over {len(pairs)} anchors: {mae:.3f} eV")
        print("A screen is usable if it keeps the ordering and lands within a few")
        print("tenths; it is not a substitute for DFT on any single number.")
    else:
        print("\nToo few anchors converged to judge the screen. Do not use this")
        print("column until the DFT campaign supplies more.")
    print(f"\nwrote {os.path.relpath(OUT, BASE)}")


if __name__ == "__main__":
    main()

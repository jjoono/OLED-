"""E_b for silver on NDP-9, by the same protocol as the rest of the screen.

Counterpoise-corrected rigid scan of a single Ag atom along the C-N axis of one
nitrile, PBE-D3BJ/def2-SVP, doublet. The molecule is frozen at its GFN2-xTB
geometry, matching how HATCN's 1.029 eV and F4TCNQ's 0.966 were obtained, so the
number lands on the same scale rather than needing a correction.

Forty-nine atoms against HATCN's thirty-one, and the SCF ladder from script 97
is reused because these silver-on-acceptor complexes oscillate under plain DIIS.
Checkpoints after every scan point.
"""
import json, os, sys, time

import numpy as np
import psi4

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathgeom import RUNS, STRUCT, H2EV, read_xyz

OUT = os.path.join(RUNS, "ndp9_binding_eV.json")
SCAN = [2.10, 2.20, 2.30, 2.40, 2.55, 2.75, 3.10]
psi4.set_memory(f"{os.environ.get('PSI4_MEM_GB', '9')} GB")
try:
    psi4.set_num_threads(len(os.sched_getaffinity(0)))
except AttributeError:
    psi4.set_num_threads(os.cpu_count() or 4)
psi4.core.set_output_file(os.path.join(RUNS, "psi4_ndp9.out"), False)

BASE = {"basis": "def2-svp", "scf_type": "df", "reference": "uks"}
LADDER = [
    dict(BASE, maxiter=150, guess="sad", level_shift=1.0,
         level_shift_cutoff=1e-3, damping_percentage=15.0),
    dict(BASE, maxiter=200, guess="sad", level_shift=2.0,
         level_shift_cutoff=1e-2, damping_percentage=40.0),
    dict(BASE, maxiter=200, guess="sad", soscf=True, soscf_start_convergence=1e-2),
]


def geom(syms, xyz, mult, ghost=()):
    st = f"0 {mult}\n"
    for i, (a, c) in enumerate(zip(syms, xyz)):
        tag = f"@{a}" if i in ghost else a
        st += f"{tag} {c[0]:.8f} {c[1]:.8f} {c[2]:.8f}\n"
    return st + "symmetry c1\nno_reorient\nno_com\n"


def energy(syms, xyz, mult, ghost=()):
    g = geom(syms, xyz, mult, ghost)
    for opts in LADDER:
        try:
            psi4.set_options(opts)
            e = psi4.energy("pbe-d3bj", molecule=psi4.geometry(g))
            psi4.core.clean()
            return e
        except Exception:
            psi4.core.clean()
    return None


def main():
    syms, xyz = read_xyz(os.path.join(STRUCT, "NDP9.xyz"))
    # pick a nitrile: N with exactly one heavy neighbour
    nit = []
    for i, q in enumerate(syms):
        if q != "N":
            continue
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        d[i] = 9
        if sum(1 for j, dd in enumerate(d) if dd < 1.75 and syms[j] != "H") == 1:
            nit.append(i)
    n_i = nit[0]
    d = np.linalg.norm(xyz - xyz[n_i], axis=1)
    d[n_i] = 9
    c_i = int(np.argmin(d))
    axis = xyz[n_i] - xyz[c_i]
    axis /= np.linalg.norm(axis)
    print(f"NDP-9: {len(syms)} atoms, {len(nit)} nitriles, "
          f"binding at N{n_i} along C{c_i}-N{n_i}\n", flush=True)

    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for r in SCAN:
        key = f"{r:.2f}"
        if key in res:
            print(f"  r = {key} A  cached", flush=True)
            continue
        t0 = time.time()
        ag = xyz[n_i] + r * axis
        full = np.vstack([xyz, ag])
        s_full = syms + ["Ag"]
        n = len(syms)
        e_c = energy(s_full, full, 2)
        # counterpoise: each fragment in the full basis of the complex
        e_m = energy(s_full, full, 1, ghost=(n,))
        e_a = energy(s_full, full, 2, ghost=tuple(range(n)))
        if None in (e_c, e_m, e_a):
            print(f"  r = {key} A  SCF failed on a fragment", flush=True)
            res[key] = None
        else:
            eb = (e_m + e_a - e_c) * H2EV
            res[key] = dict(Eb_eV=round(eb, 4), E_complex=e_c, E_mol=e_m, E_Ag=e_a)
            print(f"  r = {key} A  E_b = {eb:+.4f} eV   [{time.time()-t0:.0f} s]",
                  flush=True)
        json.dump(res, open(OUT, "w"), indent=1)

    ok = {k: v["Eb_eV"] for k, v in res.items() if isinstance(v, dict)}
    if ok:
        best = max(ok, key=lambda k: ok[k])
        print(f"\nE_b(NDP-9) = {ok[best]:.3f} eV at r = {best} A")
        print(f"  HATCN 1.029, F4TCNQ 0.966 by the same protocol")
        if float(best) in (SCAN[0], SCAN[-1]):
            print("  WARNING: the minimum sits at the edge of the scan; extend it")


if __name__ == "__main__":
    main()

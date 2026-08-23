"""E_d with Gaussian instead of psi4 -- same path, same definition, better SCF.

The barrier campaign has been blocked on convergence, not on cost: an open-shell
Ag on a strong acceptor has near-degenerate states, and plain DIIS oscillates
indefinitely. Gaussian's quadratically convergent SCF exists for exactly that,
so a workstation with Gaussian is a better machine for this than one with psi4.

The path is imported from pathgeom, unchanged, so a barrier from this driver and
one from script 97 describe the same physical process and can sit in the same
table. What differs is the program and therefore the absolute energies, which is
harmless -- E_d is a difference between two points computed the same way. Say
which program produced which row when reporting.

  method   PBEPBE/Def2SVP, EmpiricalDispersion=GD3BJ, charge 0 multiplicity 2
           Gaussian folds the dispersion correction into the nuclear repulsion,
           so the SCF Done energy already includes it.
  SCF      SCF=(XQC,MaxCycle=128): ordinary DIIS first, quadratic convergence
           automatically if that stalls.
  needs    G16. The Def2SVP keyword and its ECP for Ag are built in from G16
           onward; G09 has no def2 basis sets and would need an external
           basis-set block, which this driver does not write.

Run it exactly like the psi4 one:

    python scripts/102_ed_gaussian.py                # everything
    python scripts/102_ed_gaussian.py pyridine       # one system

Set GAUSS_EXE if the binary is not called g16 (g16.exe on Windows is found
automatically if it is on PATH).

A many-core workstation finishes sooner running several instances of this on
different systems than one instance with every thread, since shared-memory
scaling on a thirty-atom molecule flattens out well before sixty-four cores.
Give each instance its own output file and merge afterwards:

    set GAUSS_NPROC=16
    set GAUSS_OUT=part1.json
    python scripts\102_ed_gaussian.py HATCN F4TCNQ TPBi p-bPPhenB Cs2CO3 B3PyMPM Bphen

    python scripts\103_merge_Ed.py           # combines every part*.json

    python scripts/102_ed_gaussian.py --selftest     # no Gaussian needed

checks the input generation and the log parser against a canned log, which is
the part that can be verified without a license.
"""
import json, os, re, shutil, subprocess, sys, tempfile, time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathgeom import (CANDIDATES, H2EV, NPATH_MAX, NPATH_MIN, RUNS, SPACING,
                      STRUCT, ZSCAN, destination, geometry, place, read_xyz)

# Per-run output, so several instances can share a machine without overwriting
# each other. Each writes the whole dictionary it knows about, so two processes
# pointed at one file would each erase the other's systems.
OUT = os.environ.get("GAUSS_OUT") or os.path.join(
    RUNS, "diffusion_barriers_gaussian.json")
if not os.path.isabs(OUT):
    OUT = os.path.join(RUNS, OUT)
EXE = os.environ.get("GAUSS_EXE", "g16")
NPROC = int(os.environ.get("GAUSS_NPROC", os.cpu_count() or 4))
MEM_GB = int(os.environ.get("GAUSS_MEM_GB", "16"))
ROUTE = ("#P PBEPBE/Def2SVP EmpiricalDispersion=GD3BJ SP "
         "SCF=(XQC,MaxCycle=128) NoSymm")

# "SCF Done:  E(UPBEPBE) =  -1487.1234567     A.U. after  23 cycles"
E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")


def gjf(syms, xyz, title, chrg=0, mult=2):
    lines = [f"%NProcShared={NPROC}", f"%Mem={MEM_GB}GB", ROUTE, "",
             title, "", f"{chrg} {mult}"]
    for a, c in zip(syms, xyz):
        lines.append(f" {a:<2s} {c[0]:14.8f} {c[1]:14.8f} {c[2]:14.8f}")
    return "\n".join(lines) + "\n\n"


def parse(log):
    """Last converged SCF energy in the log, or None."""
    hits = E_RE.findall(log)
    return float(hits[-1]) if hits else None


def energy(syms, xyz, title):
    d = tempfile.mkdtemp(prefix="g16_")
    try:
        inp = os.path.join(d, "j.gjf")
        with open(inp, "w") as f:
            f.write(gjf(syms, xyz, title))
        try:
            r = subprocess.run([EXE, "j.gjf"], cwd=d, capture_output=True,
                               text=True, timeout=7200)
        except FileNotFoundError:
            print(f"    {EXE!r} not found on PATH -- set GAUSS_EXE", flush=True)
            raise SystemExit(2)
        except subprocess.TimeoutExpired:
            return None
        logp = os.path.join(d, "j.log")
        if not os.path.exists(logp):
            return parse(r.stdout)
        log = open(logp, errors="replace").read()
        if "Normal termination" not in log:
            return None
        return parse(log)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def barrier(tag, fn, rule, mult):
    syms, xyz = read_xyz(os.path.join(STRUCT, fn))
    sub_s, sub_x, ag, anchor, nrm = geometry(syms, xyz)
    dest, cls = destination(sub_s, sub_x, ag, anchor, rule)
    span = float(np.linalg.norm(dest - ag))
    npath = int(np.clip(round(span / SPACING) + 1, NPATH_MIN, NPATH_MAX))
    print(f"[{tag}] {len(sub_s)} substrate atoms, anchor {sub_s[anchor]}{anchor}, "
          f"path {cls}, span {span:.2f} A, {npath} points", flush=True)

    E = []
    for t in np.linspace(0.0, 1.0, npath):
        pos = place(sub_x, (1 - t) * ag + t * dest, nrm)
        best = None
        for dz in ZSCAN:
            e = energy(sub_s + ["Ag"], np.vstack([sub_x, pos + dz * nrm]),
                       f"{tag} t={t:.2f} dz={dz:+.2f}")
            if e is not None and (best is None or e < best):
                best = e
        E.append(best)
        print(f"    t={t:.2f}  E={'FAIL' if best is None else f'{best:.6f}'}",
              flush=True)

    ok = [e for e in E if e is not None]
    if E[0] is None or len(ok) < 3:
        print(f"[{tag}] insufficient converged points", flush=True)
        return None
    ed = (max(ok) - E[0]) * H2EV
    print(f"[{tag}] E_d = {ed:.3f} eV  ({cls})", flush=True)
    return {"E_d_eV": round(ed, 4), "class": cls, "anchor": sub_s[anchor],
            "n_atoms": len(sub_s), "path_A": round(span, 3), "n_points": npath,
            "program": "gaussian", "route": ROUTE}


def selftest():
    syms, xyz = read_xyz(os.path.join(STRUCT, "pyridine_Ag.xyz"))
    text = gjf(syms, xyz, "selftest")
    problems = []
    if not text.startswith("%NProcShared"):
        problems.append("link0 must come first")
    if text.count("\n\n") < 3:
        problems.append("missing a blank separator line")
    body = text.split("\n\n")
    if len(body) < 4:
        problems.append("route/title/molecule blocks not separated")
    if f"{len(syms)}" and len([l for l in text.splitlines()
                               if re.match(r"^ [A-Z][a-z]?\s+-?\d", l)]) != len(syms):
        problems.append("atom count in the input does not match the structure")
    canned = ("blah\n SCF Done:  E(UPBEPBE) =  -1487.12345678     A.U. after 23 cycles\n"
              " SCF Done:  E(UPBEPBE) =  -1487.22222222     A.U. after 31 cycles\n"
              " Normal termination of Gaussian 16\n")
    got = parse(canned)
    if got != -1487.22222222:
        problems.append(f"parser returned {got}, expected the LAST SCF Done")
    if parse("no energy here") is not None:
        problems.append("parser should return None when nothing matched")
    print(text.split("\n\n")[0])
    print(ROUTE)
    print(f"...{len(syms)} atoms written")
    if problems:
        print("\nSELFTEST FAILED:")
        for p in problems:
            print("  -", p)
        return 1
    print("\nselftest passed: input generation and log parsing are consistent.")
    print("The SCF itself cannot be checked without a Gaussian license.")
    return 0


def main():
    args = [a for a in sys.argv[1:]]
    if "--selftest" in args:
        raise SystemExit(selftest())
    only = args or None
    if only:
        known = {t for t, _, _, _ in CANDIDATES}
        unknown = [a for a in only if a not in known]
        if unknown:
            print(f"not candidate names, ignored: {' '.join(unknown)}")
            print(f"known: {' '.join(sorted(known))}\n", flush=True)
    print(f"writing to {OUT}\n", flush=True)
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for tag, fn, rule, mult in CANDIDATES:
        if only and tag not in only:
            continue
        if res.get(tag) is not None:
            print(f"[{tag}] already done, skipping", flush=True)
            continue
        if not os.path.exists(os.path.join(STRUCT, fn)):
            print(f"[{tag}] missing {fn}", flush=True)
            continue
        t0 = time.time()
        try:
            res[tag] = barrier(tag, fn, rule, mult)
        except SystemExit:
            raise
        except Exception as exc:
            print(f"[{tag}] aborted: {type(exc).__name__}: {exc}", flush=True)
            res[tag] = None
        print(f"[{tag}] {time.time()-t0:.0f} s\n", flush=True)
        json.dump(res, open(OUT, "w"), indent=1)
    print(json.dumps({k: (v or {}).get("E_d_eV") for k, v in res.items()}, indent=1))


if __name__ == "__main__":
    main()

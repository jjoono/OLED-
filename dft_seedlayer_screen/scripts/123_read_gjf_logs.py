"""Turn the Gaussian logs back into barriers.

Accepts either extension: Gaussian 16W on Windows writes .out, the Linux build
writes .log, and the same folder may contain both if it was run in two places.

Runs here, not on the workstation, so nothing has to be installed there beyond
Gaussian. Point it at the returned gaussian_jobs folder.

    python scripts/123_read_gjf_logs.py [folder]

For each candidate it takes the lowest converged energy over the height scan at
every path point -- the relaxed-adatom approximation -- and reports
max(path) - E(start) in eV, with the path class and a note on any point that
failed to converge. A candidate missing its starting point, or with fewer than
three converged points, is reported as incomplete rather than given a number.
"""
import json, os, re, sys

H2EV = 27.211386
E_RE = re.compile(r"SCF Done:\s+E\(\S+\)\s*=\s*(-?\d+\.\d+)")
S2_RE = re.compile(r"S\*\*2 *= *(-?\d+\.\d+)")
NAME_RE = re.compile(r"_t(\d+)_z(\d+)(_ref)?$")

# A doublet has S**2 = 0.75 exactly. UKS values drift a little above it for any
# real open-shell complex, but in the v5 run the points that had collapsed onto
# a wrong SCF solution were the contaminated ones every time -- benzene 0.7567
# at +3.41 eV, DMABN 0.7534 at +2.20 eV, Bphen 0.7551 -- while every point that
# agreed with its neighbours sat between 0.7512 and 0.7519. The separation is
# clean enough to reject on, and rejecting is right: a point on a different
# electronic state contributes the gap between states, not a barrier.
S2_MAX = 0.7525


def read_log(p):
    """(energy, S**2) of the last converged SCF, or (None, None) if unusable."""
    try:
        txt = open(p, errors="replace").read()
    except OSError:
        return None, None
    if "Normal termination" not in txt:
        return None, None
    hits = E_RE.findall(txt)
    if not hits:
        return None, None
    s2 = S2_RE.findall(txt)
    return float(hits[-1]), (float(s2[-1]) if s2 else None)


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "gaussian_jobs")
    if not os.path.isdir(root):
        print(f"no such folder: {root}")
        raise SystemExit(1)

    cls = {}
    idx = os.path.join(root, "INDEX.csv")
    if os.path.exists(idx):
        for line in open(idx).read().splitlines()[1:]:
            f = line.split(",")
            if len(f) >= 3:
                cls[f[1]] = (f[0], f[2])

    results, missing = {}, []
    for folder in sorted(os.listdir(root)):
        d = os.path.join(root, folder)
        if not os.path.isdir(d):
            continue
        pts, failed, total, contaminated = {}, 0, 0, []
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".gjf"):
                continue
            total += 1
            stem = fn[:-4]
            mo = NAME_RE.search(stem)
            if not mo:
                continue
            t = int(mo.group(1))
            e = s2 = None
            for ext in (".out", ".log"):      # G16W writes .out, g16 on Linux .log
                cand = os.path.join(d, stem + ext)
                if os.path.exists(cand):
                    e, s2 = read_log(cand)
                    break
            if e is None:
                failed += 1
                continue
            if s2 is not None and s2 > S2_MAX:
                contaminated.append((stem, s2))
                continue
            # The re-run of t=0 shares its path index with the original, so the
            # usual lowest-of-the-height-scan rule already keeps whichever of the
            # two landed on the lower solution.
            if t not in pts or e < pts[t]:
                pts[t] = e

        tag, path_class = cls.get(folder, (folder, "?"))
        if 0 not in pts or len(pts) < 3:
            missing.append((tag, len(pts), total, failed))
            continue
        order = sorted(pts)
        ed = (max(pts.values()) - pts[0]) * H2EV
        # A path that only goes downhill has no saddle on it, and the barrier it
        # reports is zero by construction rather than by physics. That happens
        # when the two endpoints are not the same site -- an asymmetric cluster
        # whose "nearest equivalent atom" is in a different environment -- or
        # when the grid is too coarse to resolve the bump between them.
        mono = all(pts[order[i + 1]] <= pts[order[i]] for i in range(len(order) - 1))
        lowest = min(pts, key=pts.get)
        results[tag] = {"E_d_eV": round(ed, 4), "class": path_class,
                        "points": len(pts), "failed_jobs": failed,
                        "program": "gaussian",
                        "monotonic_downhill": mono,
                        "rejected_spin_contaminated": [n for n, _ in contaminated],
                        "lowest_point": int(lowest),
                        "endpoint_gap_eV": round((pts[order[-1]] - pts[0]) * H2EV, 4),
                        "usable": not mono}

    if results:
        print(f"{'candidate':<12} {'E_d (eV)':>9} {'end-start':>10} "
              f"{'class':<10} {'pts':>4} {'failed':>7}")
        print("-" * 60)
        for tag, v in sorted(results.items(), key=lambda kv: -kv[1]["E_d_eV"]):
            flag = "" if v["usable"] else "  <- downhill, no saddle on this path"
            # A site2site path now runs between two equivalent atoms and ends at
            # the same height it started, so its two endpoints are the same
            # state by symmetry and their energy difference is a free error bar
            # on the barrier above it. A gap that is not small next to E_d means
            # the two ends did not converge to the same electronic state, and
            # the barrier inherits that error whatever the path looks like.
            if v["class"] == "site2site" and abs(v["endpoint_gap_eV"]) > \
                    max(0.05, 0.3 * abs(v["E_d_eV"])):
                flag += "  <- endpoints disagree; not one state"
            print(f"{tag:<12} {v['E_d_eV']:>9.3f} {v['endpoint_gap_eV']:>+10.3f} "
                  f"{v['class']:<10} {v['points']:>4} {v['failed_jobs']:>7}{flag}")
        bad = [t for t, v in results.items() if not v["usable"]]
        if bad:
            print(f"\n{len(bad)} candidate(s) have no barrier on their path and need")
            print("a denser one, or endpoints that are genuinely the same site:")
            for t_ in bad:
                v = results[t_]
                print(f"  {t_:<12} endpoint is {v['endpoint_gap_eV']:+.3f} eV from "
                      f"the start -- the two sites are not equivalent")
    if any(v["rejected_spin_contaminated"] for v in results.values()):
        print("\nrejected for spin contamination (S**2 > "
              f"{S2_MAX}) -- wrong SCF solution:")
        for t_, v in sorted(results.items()):
            for n in v["rejected_spin_contaminated"]:
                print(f"  {n}")
    if missing:
        print("\nincomplete (no barrier reported):")
        for tag, got, total, failed in missing:
            print(f"  {tag:<12} {got} usable path points, {failed}/{total} jobs failed")

    if "HATCN" in results:
        print(f"\nCHECK: HATCN came out at {results['HATCN']['E_d_eV']:.3f} eV. "
              f"psi4 gave 0.286 on the same")
        print("  path. Agreement there is what licenses the rest of the column.")

    if not results:
        # An empty results file would read later as "the campaign ran and found
        # nothing", which is not what an unrun campaign means.
        print("\nNo barriers to write. Run RUN_ALL.bat on the workstation first,")
        print("then point this script at the returned folder.")
        return
    out = os.path.join(os.path.dirname(root), "runs",
                       "diffusion_barriers_gaussian.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(results, open(out, "w"), indent=1)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

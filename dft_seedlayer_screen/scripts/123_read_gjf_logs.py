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
# Gaussian prints "S**2 =" only in the quadratic-convergence summary, so a
# pattern that matches only that form reads nothing at all from the jobs that
# converged under DIIS -- 71 of the first 92 v8 outputs. The annihilation line
# is printed for every UKS job.
S2_RE = re.compile(r"S\*\*2 before annihilation\s+(\d+\.\d+)")
S2_ALT = re.compile(r"<S\*\*2>=\s*(\d+\.\d+)")
NAME_RE = re.compile(r"_t(\d+)_z(\d+)(_ref)?$")
DZ_RE = re.compile(r"dz=([+-]?\d+\.\d+)")


# A doublet has S**2 = 0.75 exactly, and a point that has collapsed onto a
# different electronic state shows it. But there is no one number to test
# against: the fixed 0.7525 calibrated on PBE rejects 21 of the first 92 PBE0
# jobs, including ones whose energies agree with their neighbours to a
# milli-Hartree, because a hybrid carries more contamination everywhere.
#
# What separates a bad point from its folder is that it differs from the rest of
# that folder, so the test is made relative. The absolute number is kept only as
# a caveat on a whole folder: F4TCNQ runs 0.7759-0.8455 across every one of its
# points, which is a real open-shell Ag-to-acceptor charge transfer rather than
# one stray job, and a single-reference doublet describes it poorly.
S2_SPREAD = 0.02      # above the folder median -> a different state
S2_FOLDER = 0.77      # folder median above this -> flag, do not reject


def relaxed(hs):
    """Energy at the relaxed adatom height, from a scan of (dz, E) pairs.

    Taking the lower of two heights does not relax anything: whichever wins is
    an endpoint of the scan, so the minimum is never bracketed. In the v8 run
    that mattered by up to 0.456 eV -- more than any barrier measured -- and
    which height won flipped along a single path, so it put structure into the
    curve that was not in the physics.

    With three or more heights the minimum can be bracketed and fitted. A
    parabola through the lowest point and its two neighbours gives the relaxed
    energy; if the lowest point is still an end of the scan the minimum lies
    outside it, and the raw value is returned with a flag so the caller can say
    so rather than quoting a fit that extrapolated.
    """
    hs = sorted(hs)
    i = min(range(len(hs)), key=lambda k: hs[k][1])
    if len(hs) < 3 or i == 0 or i == len(hs) - 1:
        return hs[i][1], False
    (x0, y0), (x1, y1), (x2, y2) = hs[i - 1], hs[i], hs[i + 1]
    d = (x0 - x1) * (x0 - x2) * (x1 - x2)
    if abs(d) < 1e-12:
        return y1, False
    a = (x2 * (y1 - y0) + x1 * (y0 - y2) + x0 * (y2 - y1)) / d
    b = (x2 * x2 * (y0 - y1) + x1 * x1 * (y2 - y0) + x0 * x0 * (y1 - y2)) / d
    c = (x1 * x2 * (x1 - x2) * y0 + x2 * x0 * (x2 - x0) * y1
         + x0 * x1 * (x0 - x1) * y2) / d
    if a <= 0:                       # not a minimum -- do not extrapolate off it
        return y1, False
    xv = -b / (2 * a)
    if not (x0 <= xv <= x2):
        return y1, False
    return c - b * b / (4 * a), True


def read_log(p):
    """(energy, S**2, dz) of the last converged SCF, or Nones if unusable."""
    try:
        txt = open(p, errors="replace").read()
    except OSError:
        return None, None, None
    if "Normal termination" not in txt:
        return None, None, None
    hits = E_RE.findall(txt)
    if not hits:
        return None, None, None
    s2 = S2_RE.findall(txt) or S2_ALT.findall(txt)
    dz = DZ_RE.findall(txt)
    return (float(hits[-1]), (float(s2[-1]) if s2 else None),
            (float(dz[0]) if dz else None))


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
        pts, failed, total, contaminated, seen = {}, 0, 0, [], []
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".gjf"):
                continue
            total += 1
            stem = fn[:-4]
            mo = NAME_RE.search(stem)
            if not mo:
                continue
            t = int(mo.group(1))
            e = s2 = dz = None
            for ext in (".out", ".log"):      # G16W writes .out, g16 on Linux .log
                cand = os.path.join(d, stem + ext)
                if os.path.exists(cand):
                    e, s2, dz = read_log(cand)
                    break
            if e is None:
                failed += 1
                continue
            seen.append((t, stem, e, s2, dz))

        # The spin test needs the whole folder before it can say what is normal
        # for it, so the points are collected first and filtered here.
        s2s = sorted(x[3] for x in seen if x[3] is not None)
        med = s2s[len(s2s) // 2] if s2s else None
        folder_hot = med is not None and med > S2_FOLDER
        scan = {}
        for t, stem, e, s2, dz in seen:
            if med is not None and s2 is not None and not folder_hot \
                    and s2 > med + S2_SPREAD:
                contaminated.append((stem, s2))
                continue
            # The re-run of t=0 repeats a height already in the scan, so it is
            # folded in as the lower of the two rather than as a new height.
            key = dz if dz is not None else len(scan.get(t, []))
            b = scan.setdefault(t, {})
            b[key] = min(b[key], e) if key in b else e
        unbracketed = []
        for t, b in scan.items():
            pts[t], fitted = relaxed(list(b.items()))
            if not fitted:
                unbracketed.append(t)

        tag, path_class = cls.get(folder, (folder, "?"))
        if 0 not in pts or len(pts) < 3:
            missing.append((tag, len(pts), total, failed))
            continue
        order = sorted(pts)
        ed = (max(pts.values()) - pts[0]) * H2EV
        # A hop between equivalent sites is symmetric about its midpoint, so
        # E(t) and E(1-t) must agree. Half of the largest disagreement is the
        # barrier's own error bar, measured on the same footing as the barrier.
        n = len(order)
        asym = max(abs(pts[order[i]] - pts[order[n - 1 - i]])
                   for i in range(n)) * H2EV / 2 if path_class == "site2site" else 0.0
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
                        "spin_contaminated_folder": folder_hot,
                        "s2_median": med,
                        "path_minimum": int(min(pts, key=pts.get)),
                        "height_unbracketed_points": sorted(unbracketed),
                        "asymmetry_eV": round(asym, 4),
                        "resolved": bool(ed > 3 * asym) if asym > 0 else None,
                        "lowest_point": int(lowest),
                        "endpoint_gap_eV": round((pts[order[-1]] - pts[0]) * H2EV, 4),
                        "usable": not mono}

    if results:
        # Screening does not need the absolute barrier, and the absolute
        # barrier is the part this protocol gets least right: a rigid substrate,
        # a small basis and a single molecule instead of a film all push the
        # same way on every candidate. Differences against a common reference
        # keep what screening uses and drop most of what it cannot trust.
        ref_tag = next((t for t in ("benzene", "pyridine", "HATCN")
                        if t in results), None)
        ref_ed = results[ref_tag]["E_d_eV"] if ref_tag else None
        print(f"{'candidate':<12} {'E_d (eV)':>9} {'+/-':>7} "
              f"{'vs ' + (ref_tag or '-'):>11} {'class':<10} {'pts':>4} {'fail':>5}")
        print("-" * 72)
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
            if v["spin_contaminated_folder"]:
                flag += f"  <- S**2 ~ {v['s2_median']:.3f} throughout"
            # The path minimum landing between the endpoints means the adatom
            # prefers somewhere the input structure did not put it, so the site
            # the barrier was measured from is not a binding site at all.
            if 0 < v["path_minimum"] < v["points"] - 1:
                flag += f"  <- minimum at point {v['path_minimum']}, not an end"
            if v["height_unbracketed_points"]:
                flag += (f"  <- height not bracketed at "
                         f"{v['height_unbracketed_points']}")
            if v["resolved"] is False:
                flag += "  <- barrier below its own error"
            rel = (f"{v['E_d_eV'] - ref_ed:>+11.3f}" if ref_ed is not None
                   else f"{'-':>11}")
            print(f"{tag:<12} {v['E_d_eV']:>9.3f} {v['asymmetry_eV']:>7.3f} "
                  f"{rel} {v['class']:<10} {v['points']:>4} "
                  f"{v['failed_jobs']:>5}{flag}")
        bad = [t for t, v in results.items() if not v["usable"]]
        if bad:
            print(f"\n{len(bad)} candidate(s) have no barrier on their path and need")
            print("a denser one, or endpoints that are genuinely the same site:")
            for t_ in bad:
                v = results[t_]
                print(f"  {t_:<12} endpoint is {v['endpoint_gap_eV']:+.3f} eV from "
                      f"the start -- the two sites are not equivalent")
    if any(v["rejected_spin_contaminated"] for v in results.values()):
        print("\nrejected: S**2 more than "
              f"{S2_SPREAD} above the folder median -- a different SCF state:")
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

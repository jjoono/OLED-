"""Turn the Gaussian logs back into barriers.

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
NAME_RE = re.compile(r"_t(\d+)_z(\d+)$")


def read_log(p):
    """Last converged SCF energy, or None if the job did not terminate normally."""
    try:
        txt = open(p, errors="replace").read()
    except OSError:
        return None
    if "Normal termination" not in txt:
        return None
    hits = E_RE.findall(txt)
    return float(hits[-1]) if hits else None


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
        pts, failed, total = {}, 0, 0
        for fn in sorted(os.listdir(d)):
            if not fn.endswith(".gjf"):
                continue
            total += 1
            stem = fn[:-4]
            mo = NAME_RE.search(stem)
            if not mo:
                continue
            t = int(mo.group(1))
            e = read_log(os.path.join(d, stem + ".log"))
            if e is None:
                failed += 1
                continue
            if t not in pts or e < pts[t]:
                pts[t] = e

        tag, path_class = cls.get(folder, (folder, "?"))
        if 0 not in pts or len(pts) < 3:
            missing.append((tag, len(pts), total, failed))
            continue
        ed = (max(pts.values()) - pts[0]) * H2EV
        results[tag] = {"E_d_eV": round(ed, 4), "class": path_class,
                        "points": len(pts), "failed_jobs": failed,
                        "program": "gaussian"}

    if results:
        print(f"{'candidate':<12} {'E_d (eV)':>9} {'class':<10} {'pts':>4} {'failed':>7}")
        print("-" * 48)
        for tag, v in sorted(results.items(), key=lambda kv: -kv[1]["E_d_eV"]):
            print(f"{tag:<12} {v['E_d_eV']:>9.3f} {v['class']:<10} "
                  f"{v['points']:>4} {v['failed_jobs']:>7}")
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

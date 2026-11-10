"""Combine the E_d results from however many runs produced them.

The campaign can be split across parallel Gaussian instances, across the psi4
driver and the Gaussian one, and across restarts. Each writes its own JSON. This
puts them in one table and says where every number came from, since a barrier
from Gaussian and one from psi4 are the same quantity computed by different
programs and that belongs in the record.

    python scripts/103_merge_Ed.py                    # every *.json in runs/
                                                     # except the rejected xTB one
    python scripts/103_merge_Ed.py part1.json part2.json

Writes runs/Ed_merged.json and prints the table, sorted by binding class then
barrier. A system present in more than one file is reported once, with the
disagreement shown rather than silently resolved.
"""
import glob, json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pathgeom import CANDIDATES, RUNS

OUT = os.path.join(RUNS, "Ed_merged.json")


def load(paths):
    seen = {}
    for p in paths:
        try:
            d = json.load(open(p))
        except (OSError, ValueError):
            continue
        if not isinstance(d, dict):
            continue          # runs/ holds plenty of JSON that is not this
        src = os.path.basename(p)
        for tag, v in d.items():
            if not isinstance(v, dict) or v.get("E_d_eV") is None:
                continue
            seen.setdefault(tag, []).append((src, v))
    return seen


def main():
    if len(sys.argv) > 1:
        paths = [p if os.path.isabs(p) else os.path.join(RUNS, p)
                 for p in sys.argv[1:]]
    else:
        paths = sorted(glob.glob(os.path.join(RUNS, "*.json")))
        # The xTB screen was rejected: it returns 0.000 eV on both systems where
        # a DFT barrier exists, against 0.286 and 0.170. Those numbers must not
        # find their way back in through a wildcard.
        paths = [p for p in paths
                 if os.path.basename(p) != os.path.basename(OUT)
                 and "xtb" not in os.path.basename(p).lower()]
    seen = load(paths)
    if not seen:
        print(f"no results found in: {', '.join(os.path.basename(p) for p in paths)}")
        return

    merged, conflicts = {}, []
    for tag, entries in seen.items():
        vals = {round(v["E_d_eV"], 3) for _, v in entries}
        src, v = entries[0]
        merged[tag] = dict(v, sources=[s for s, _ in entries])
        if len(vals) > 1:
            conflicts.append((tag, [(s, x["E_d_eV"], x.get("program", "psi4"))
                                    for s, x in entries]))

    order = {t: i for i, (t, _, _, _) in enumerate(CANDIDATES)}
    rows = sorted(merged.items(),
                  key=lambda kv: (kv[1].get("class", ""), -kv[1]["E_d_eV"]))
    print(f"{'candidate':<12} {'E_d (eV)':>9} {'class':<10} {'anchor':<7} "
          f"{'program':<9} source")
    print("-" * 72)
    for tag, v in rows:
        print(f"{tag:<12} {v['E_d_eV']:>9.3f} {v.get('class',''):<10} "
              f"{v.get('anchor',''):<7} {v.get('program','psi4'):<9} "
              f"{','.join(v['sources'])}")

    done = set(merged)
    missing = [t for t, _, _, _ in CANDIDATES if t not in done]
    print(f"\n{len(done)}/{len(CANDIDATES)} candidates have a barrier")
    if missing:
        print(f"still missing: {' '.join(missing)}")
    if conflicts:
        print("\nSame system computed more than once, values differ:")
        for tag, e in conflicts:
            print(f"  {tag}: " + "  ".join(f"{v:.3f} ({s}, {p})" for s, v, p in e))
        print("  Kept the first. A Gaussian and a psi4 number for one system are")
        print("  the same quantity from different programs; a gap between them is")
        print("  the method spread and worth keeping in view, not averaging away.")

    json.dump(merged, open(OUT, "w"), indent=1)
    print(f"\nwrote {os.path.relpath(OUT, os.path.dirname(RUNS))}")


if __name__ == "__main__":
    main()

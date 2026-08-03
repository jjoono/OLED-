#!/bin/bash
# Slab E_b for Ag on the HATCN monolayer:
#     E_b = E(HATCN monolayer) + E(isolated Ag) - E(HATCN + Ag)
# Runs the Ag reference and the rigid Ag-N distance scan, then reports the curve.
# Every job is a single point; nothing here needs a geometry optimisation.
set -u
source /root/miniforge3/etc/profile.d/conda.sh
conda activate cp2k
export CP2K_DATA_DIR=/root/miniforge3/envs/cp2k/share/cp2k/data
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}

SLABS="$(cd "$(dirname "$0")/../slabs" && pwd)"

run() {                     # run <dir> <name>
  local d="$1" n="$2"
  if grep -aq "Total FORCE_EVAL" "$SLABS/$d/$n.out" 2>/dev/null; then
    echo "[skip] $n (already done)"; return
  fi
  local s=$SECONDS
  ( cd "$SLABS/$d" && cp2k.psmp -i "$n.inp" -o "$n.out" >/dev/null 2>&1 )
  local e
  e=$(grep -a "Total FORCE_EVAL" "$SLABS/$d/$n.out" 2>/dev/null | tail -1 | awk '{print $NF}')
  echo "[done] $n  $((SECONDS-s)) s  E = ${e:-FAILED}"
}

run ag_atom  ag_atom
run hatcn_ml hatcn_ml
# far -> near; each point starts from the previous one's converged wavefunction,
# which is what actually fixes the charge sloshing at short Ag-N distance.
PREV=""
for r in 3p5 3p0 2p6 2p4 2p2 2p0; do
  n="hatcn_ag_r$r"
  [ -n "$PREV" ] && cp -f "$SLABS/hatcn_ag_scan/$PREV-RESTART.wfn" \
                          "$SLABS/hatcn_ag_scan/$n-RESTART.wfn" 2>/dev/null
  run hatcn_ag_scan "$n"
  PREV="$n"
done

echo
python3 - "$SLABS" <<'PY'
import sys, os, re, glob
S = sys.argv[1]
HA = 27.211386


def energy(p):
    if not os.path.exists(p):
        return None
    # Prefer the T->0 extrapolated energy: with Fermi smearing the raw total is a
    # free energy and carries -TS, which does not cancel in E_b between systems
    # with different densities of states at E_F.
    txt = open(p, errors="ignore").read().splitlines()
    for line in reversed(txt):
        if "extrapolated to T->0" in line:
            return float(line.split()[-1])
    for line in reversed(txt):
        if "Total FORCE_EVAL" in line:
            return float(line.split()[-1])
    return None


e_ml = energy(f"{S}/hatcn_ml/hatcn_ml.out")
e_ag = energy(f"{S}/ag_atom/ag_atom.out")
print(f"E(HATCN monolayer) = {e_ml}")
print(f"E(isolated Ag)     = {e_ag}")
if e_ml is None or e_ag is None:
    print("references missing -- cannot form E_b")
    raise SystemExit

print(f"\n{'Ag-N (A)':>10}{'E_total (Ha)':>18}{'E_b (eV)':>12}")
rows = []
for f in sorted(glob.glob(f"{S}/hatcn_ag_scan/*.out")):
    m = re.search(r"r(\d+)p(\d+)", os.path.basename(f))
    if not m:
        continue
    r = float(f"{m.group(1)}.{m.group(2)}")
    e = energy(f)
    if e is None:
        print(f"{r:>10.1f}{'FAILED':>18}")
        continue
    eb = (e_ml + e_ag - e) * HA
    rows.append((r, eb))
    print(f"{r:>10.1f}{e:>18.6f}{eb:>12.3f}")

if rows:
    r_best, eb_best = max(rows, key=lambda t: t[1])
    print(f"\nslab E_b = {eb_best:.3f} eV at Ag-N = {r_best:.1f} A")
    print("cluster reference (PBE-D3BJ/def2-SVP, counterpoise): 1.03 eV at 2.3 A")
    print("\nThe two are NOT expected to agree to 0.01 eV: different basis")
    print("(DZVP-MOLOPT/GTH vs def2-SVP all-electron), and the slab number carries")
    print("no counterpoise correction, so it is biased high by residual BSSE. What")
    print("the slab tests is whether the periodic environment changes the picture")
    print("-- i.e. whether the strong nitrile anchoring survives in a real monolayer.")
PY

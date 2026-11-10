#!/usr/bin/env bash
# Install everything script 97 needs on a cluster with nothing on it.
#
# No root, no module system, no admin request: miniforge unpacks into your home
# directory and the environment lives inside it. The only requirement is
# outbound HTTPS on the machine you run this on, which login nodes normally
# have even when compute nodes do not -- install here, submit the job after.
#
#   bash setup_cluster.sh                  # installs to $HOME/miniforge3
#   PREFIX=/scratch/$USER/mf bash setup_cluster.sh    # small home quota
#
# Roughly 4 GB and ten minutes. Re-running is safe; it skips what is present.

set -euo pipefail
PREFIX="${PREFIX:-$HOME/miniforge3}"
ENV_NAME="${ENV_NAME:-qc}"

echo "==> target prefix : $PREFIX"
echo "==> environment   : $ENV_NAME"

if [ ! -x "$PREFIX/bin/conda" ]; then
  echo "==> downloading miniforge"
  URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  TMP="$(mktemp -d)"
  trap 'rm -rf "$TMP"' EXIT
  if command -v curl >/dev/null 2>&1; then
    curl -fsSL "$URL" -o "$TMP/mf.sh"
  else
    wget -q "$URL" -O "$TMP/mf.sh"
  fi
  bash "$TMP/mf.sh" -b -p "$PREFIX"
else
  echo "==> miniforge already at $PREFIX"
fi

CONDA="$PREFIX/bin/conda"

if [ ! -d "$PREFIX/envs/$ENV_NAME" ]; then
  echo "==> creating the $ENV_NAME environment (this is the slow part)"
  # psi4 and xtb are the only ones that matter for the barrier campaign; the
  # rest are what the other scripts in this repo import.
  "$CONDA" create -y -n "$ENV_NAME" -c conda-forge \
      python=3.11 psi4 xtb ase rdkit numpy scipy matplotlib openpyxl \
      simple-dftd3 dftd3-python
else
  echo "==> environment $ENV_NAME already exists"
fi

PY="$PREFIX/envs/$ENV_NAME/bin/python"

echo "==> verifying"
"$PY" - <<'PYCHK'
import sys
bad = []
for mod in ("psi4", "numpy", "scipy", "ase"):
    try:
        m = __import__(mod)
        print(f"  [ ok ] {mod:8s} {getattr(m, '__version__', '?')}")
    except Exception as e:
        print(f"  [FAIL] {mod:8s} {e}")
        bad.append(mod)
sys.exit(1 if bad else 0)
PYCHK

echo "==> a real single point, to prove the DFT stack works end to end"
"$PY" - <<'PYCHK'
import psi4
psi4.set_memory("2 GB"); psi4.set_num_threads(2)
psi4.core.be_quiet()
psi4.set_options({"basis": "def2-svp", "scf_type": "df", "reference": "uks"})
e = psi4.energy("pbe-d3bj", molecule=psi4.geometry(
    "0 2\nAg 0.0 0.0 0.0\nsymmetry c1\nno_reorient\nno_com\n"))
print(f"  [ ok ] Ag atom, PBE-D3BJ/def2-SVP: {e:.6f} Eh")
PYCHK

cat <<EOF

==> done

  python : $PY

  Run the campaign directly:
      $PY scripts/97_diffusion_barriers.py

  Or one system at a time, to see it working first:
      $PY scripts/97_diffusion_barriers.py pyridine

  Under a scheduler, edit and submit scripts/run_Ed.slurm instead.
EOF

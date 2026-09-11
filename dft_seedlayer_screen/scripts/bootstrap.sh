#!/usr/bin/env bash
# Verify a fresh local checkout can actually run this project.
#
# Written because the cloud container kept being restored from snapshots and the
# environment had to be re-verified from scratch each time. On a local machine it
# runs once, after cloning, and tells you exactly what is missing and what paths
# need editing -- several scripts carry cloud absolute paths that will not exist
# on your machine.
#
#   bash dft_seedlayer_screen/scripts/bootstrap.sh
#
# Exits non-zero if anything required is missing. Nothing here modifies the repo.

set -uo pipefail
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PROJ="$REPO_ROOT/dft_seedlayer_screen"
FAIL=0
WARN=0

say()  { printf '%s\n' "$*"; }
ok()   { printf '  [ ok ] %s\n' "$*"; }
bad()  { printf '  [FAIL] %s\n' "$*"; FAIL=$((FAIL+1)); }
warn() { printf '  [warn] %s\n' "$*"; WARN=$((WARN+1)); }

say ""
say "=============================================================="
say " Ag seed-layer project -- local environment check"
say " repo: $REPO_ROOT"
say "=============================================================="

# ---------------------------------------------------------------- conda
say ""
say "[1/6] conda"
# `conda` is usually a shell function defined by conda.sh in an interactive rc
# file, so `command -v conda` fails under a non-interactive `bash bootstrap.sh`
# even on a machine where conda works fine. Source it from the usual locations
# first, then check -- otherwise this reports both environments missing when they
# are present, which is exactly what the first run of this script did.
CONDA_SH=""
for cand in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" \
            /root/miniforge3 /opt/conda "${CONDA_PREFIX:-}"; do
    if [ -n "$cand" ] && [ -f "$cand/etc/profile.d/conda.sh" ]; then
        CONDA_SH="$cand/etc/profile.d/conda.sh"
        break
    fi
done
if [ -z "$CONDA_SH" ] && command -v conda >/dev/null 2>&1; then
    CONDA_SH="$(conda info --base 2>/dev/null)/etc/profile.d/conda.sh"
fi
if [ -n "$CONDA_SH" ] && [ -f "$CONDA_SH" ]; then
    # shellcheck disable=SC1090
    source "$CONDA_SH"
    ok "conda: $(conda --version 2>&1)  [$CONDA_SH]"
else
    bad "conda not found. Install miniforge (see LOCAL_HANDOFF.md section 2.2)."
    CONDA_SH=""
fi

# ---------------------------------------------------------------- qc env
say ""
say "[2/6] qc environment (psi4 / xtb / ase / rdkit)"
if [ -n "$CONDA_SH" ] && conda env list 2>/dev/null | grep -qE '^qc\s'; then
    ok "env 'qc' exists"
    # `conda run` does not forward stdin, so a heredoc here runs an EMPTY script
    # and silently reports success -- the first version of this check verified
    # nothing. Pass the code with -c instead.
    conda run -n qc python -c '
import importlib, sys
missing = []
for m in ("psi4", "ase", "rdkit", "numpy", "scipy"):
    try:
        mod = importlib.import_module(m)
        print("  [ ok ] %s %s" % (m, getattr(mod, "__version__", "?")))
    except Exception as e:
        missing.append(m)
        print("  [FAIL] %s: %s" % (m, type(e).__name__))
sys.exit(1 if missing else 0)
' || bad "qc env imports failed"
    if conda run -n qc which xtb >/dev/null 2>&1; then
        ok "xtb: $(conda run -n qc xtb --version 2>&1 | grep -oE 'version [0-9.]+' | head -1)"
    else
        bad "xtb not on PATH in env 'qc'"
    fi
else
    bad "env 'qc' missing. conda env create -f $PROJ/env/qc.yml"
fi

# ---------------------------------------------------------------- cp2k env
say ""
say "[3/6] cp2k environment"
CP2K_BIN=""
if [ -n "$CONDA_SH" ] && conda env list 2>/dev/null | grep -qE '^cp2k\s'; then
    CP2K_PREFIX="$(conda env list | awk '$1=="cp2k"{print $NF}')"
    CP2K_BIN="$CP2K_PREFIX/bin/cp2k.psmp"
    if [ -x "$CP2K_BIN" ]; then
        ok "cp2k.psmp: $CP2K_BIN"
        ok "version: $("$CP2K_BIN" --version 2>/dev/null | grep -oE 'CP2K version [0-9.]+' | head -1)"
    else
        bad "cp2k.psmp not executable at $CP2K_BIN"
    fi
    CP2K_DATA="$CP2K_PREFIX/share/cp2k/data"
    if [ -d "$CP2K_DATA" ]; then
        ok "CP2K_DATA_DIR: $CP2K_DATA"
    else
        bad "CP2K data dir missing: $CP2K_DATA"
    fi
else
    bad "env 'cp2k' missing. conda env create -f $PROJ/env/cp2k.yml"
fi

# ------------------------------------------------------- hardcoded paths
say ""
say "[4/6] cloud absolute paths that need editing"
# Skip __pycache__ (compiled copies, regenerated) and this script itself (its
# only "cloud paths" are the fallback list and the message strings above).
HITS=$(grep -rln '/root/miniforge3\|/home/user/OLED-' "$PROJ/scripts" 2>/dev/null \
       | grep -v '__pycache__' | grep -v 'bootstrap.sh' | sort)
if [ -z "$HITS" ]; then
    ok "no cloud paths found"
else
    warn "these files carry cloud paths and will fail until edited:"
    while IFS= read -r f; do
        printf '         %s\n' "${f#$REPO_ROOT/}"
    done <<< "$HITS"
    say ""
    say "       Most important: scripts/_ckpt.py sets"
    say "         REPO = \"/home/user/OLED-\"   ->  change to $REPO_ROOT"
    say "       The cp2k binary path appears as"
    say "         /root/miniforge3/envs/cp2k/bin/cp2k.psmp"
    [ -n "$CP2K_BIN" ] && say "         ->  change to $CP2K_BIN"
fi

# ---------------------------------------------------------------- data
say ""
say "[5/6] project data"
for d in runs structures slabs scripts env; do
    if [ -d "$PROJ/$d" ]; then
        ok "$d/ ($(find "$PROJ/$d" -type f | wc -l | tr -d ' ') files)"
    else
        bad "$d/ missing"
    fi
done
NJSON=$(find "$PROJ/runs" -name '*.json' 2>/dev/null | wc -l | tr -d ' ')
[ "$NJSON" -gt 0 ] && ok "result JSONs: $NJSON" || bad "no result JSONs in runs/"

# checkpoint of the in-flight calculation
if [ -f "$PROJ/runs/tcpm_slab.json" ]; then
    ok "in-flight TCPM slab checkpoint present:"
    python3 -c '
import json, sys
for k, v in sorted(json.load(open(sys.argv[1])).items()):
    print("         %s: %s" % (k, v))
' "$PROJ/runs/tcpm_slab.json" 2>/dev/null || warn "could not read the checkpoint"
fi

# ---------------------------------------------------------------- resources
say ""
say "[6/6] machine"
CORES=$(nproc 2>/dev/null || echo '?')
MEMGB=$(awk '/MemTotal/{printf "%.0f", $2/1048576}' /proc/meminfo 2>/dev/null || echo '?')
DISKGB=$(df -BG --output=avail "$REPO_ROOT" 2>/dev/null | tail -1 | tr -dc '0-9' || echo '?')
say "  cores: $CORES   RAM: ${MEMGB} GB   free disk: ${DISKGB} GB"
if [ "$MEMGB" != "?" ] && [ "$MEMGB" -lt 20 ] 2>/dev/null; then
    warn "under 20 GB visible. If this is WSL2, raise it in %USERPROFILE%\\.wslconfig"
    warn "(memory=24GB) -- WSL2 defaults to half the host RAM."
else
    [ "$MEMGB" != "?" ] && ok "RAM is enough to raise the CP2K grid budget:"
    say "         export CP2K_MEM_BUDGET_GB=20     # cloud was forced down to 9"
fi
[ "$CORES" != "?" ] && say "         export OMP_NUM_THREADS=$(( CORES > 12 ? 12 : CORES ))"

# ---------------------------------------------------------------- verdict
say ""
say "=============================================================="
if [ "$FAIL" -eq 0 ]; then
    say " READY  ($WARN warning(s))"
    say ""
    say " Next: read LOCAL_HANDOFF.md section 4.3 (retracted claims) before"
    say " anything else, then section 5 item 1 -- diagnosing the isolated-Ag"
    say " reference in the TCPM slab, which currently gives an impossible"
    say " +2.72 eV binding at 4.0 A."
else
    say " NOT READY -- $FAIL problem(s) above"
fi
say "=============================================================="
say ""
exit $(( FAIL > 0 ? 1 : 0 ))

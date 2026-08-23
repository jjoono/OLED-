#!/usr/bin/env bash
#SBATCH --job-name=Ed_seed
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=Ed_%j.log
#
# Edit the four resource lines above to match your partition, then:
#     sbatch scripts/run_Ed.slurm
#
# The campaign is embarrassingly serial -- one system after another, each a
# handful of single points -- so one node with many cores is the right shape.
# psi4 reads the core count from the CPU affinity mask, which SLURM sets from
# cpus-per-task, so nothing needs to be passed in.
#
# Every system is checkpointed to runs/diffusion_barriers_all.json as it
# finishes, and a re-run skips what is already there. If the job hits its wall
# clock, resubmit it and it continues.

set -euo pipefail

PREFIX="${PREFIX:-$HOME/miniforge3}"
ENV_NAME="${ENV_NAME:-qc}"
PY="$PREFIX/envs/$ENV_NAME/bin/python"

# psi4 writes large temporary integral files. Point them at node-local disk if
# there is any -- a shared filesystem makes this several times slower.
export PSI_SCRATCH="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/psi4_$SLURM_JOB_ID"
mkdir -p "$PSI_SCRATCH"
trap 'rm -rf "$PSI_SCRATCH"' EXIT

# Give psi4 most of the allocation but not all of it; it accounts for the
# integral buffers only, and overshooting gets the job killed by the scheduler.
export PSI4_MEM_GB="${PSI4_MEM_GB:-$(( ${SLURM_MEM_PER_NODE:-65536} / 1024 * 8 / 10 ))}"

cd "$(dirname "$0")/.."
echo "python   : $PY"
echo "cores    : $(nproc --all) allocated $(python3 -c 'import os;print(len(os.sched_getaffinity(0)))' 2>/dev/null || echo '?')"
echo "psi4 mem : ${PSI4_MEM_GB} GB"
echo "scratch  : $PSI_SCRATCH"

srun "$PY" scripts/97_diffusion_barriers.py

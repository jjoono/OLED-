#!/bin/bash
# usage: ./run_cp2k.sh <jobdir>/<name>.inp
source /root/miniforge3/etc/profile.d/conda.sh
conda activate cp2k
export CP2K_DATA_DIR=/root/miniforge3/envs/cp2k/share/cp2k/data
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-4}
inp="$1"; cd "$(dirname "$inp")"
cp2k.psmp -i "$(basename "$inp")" -o "$(basename "${inp%.inp}").out"

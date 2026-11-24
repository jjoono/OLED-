#!/bin/bash
cd "$(dirname "$0")/../.."
while pgrep -f "tddft_anion12.sh" >/dev/null; do sleep 60; done
TD_ANION_NROOTS=12 OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=4000 python3 scripts/166_design_tddft.py Hz-O-Hz >> runs/design/tddft_anion12.log 2>&1

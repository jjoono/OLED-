#!/bin/bash
# anion with 12 roots (reaches the visible) for the heptazines and HATCN as the same-level reference
cd "$(dirname "$0")/../.."
while pgrep -f "tddft_hz.sh|166_design_tddft.py Hz-NCN-Hz CN3" >/dev/null; do sleep 60; done
TD_ANION_NROOTS=12 OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=4000 python3 scripts/166_design_tddft.py Hz-NCN-Hz HATCN >> runs/design/tddft_anion12.log 2>&1

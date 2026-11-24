#!/bin/bash
# after the heptazine-dimer PBE0 job: transparency (TDA-PBE0, neutral + anion)
cd "$(dirname "$0")/../.."
while pgrep -f "177_scaffold_pbe0.py Hz-NCN-Hz" >/dev/null; do sleep 60; done
OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=4000 python3 scripts/166_design_tddft.py Hz-NCN-Hz CN3-heptazine Hz-O-Hz >> runs/design/tddft_hz.log 2>&1

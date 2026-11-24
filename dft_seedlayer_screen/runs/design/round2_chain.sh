#!/bin/bash
# after the heptazine TD-DFT job: PBE0 E_b, then TD-DFT, for the round-2 picks
cd "$(dirname "$0")/../.."
# (no wait: started directly once the CPU was free)
SCAFFOLD_OUT=runs/design/scaffold_pbe0_hz3.json OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 \
  python3 scripts/177_scaffold_pbe0.py "Hz-C(CN)2-Hz" Hz-S-Hz >> runs/design/round2.log 2>&1
OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=4000 python3 scripts/166_design_tddft.py "Hz-C(CN)2-Hz" Hz-S-Hz >> runs/design/round2.log 2>&1

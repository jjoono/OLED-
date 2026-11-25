#!/bin/bash
# after the restart: Hz-S-Hz remaining PBE0 sites, then its TD-DFT (neutral, anion 12 roots)
cd "$(dirname "$0")/../.."
SCAFFOLD_OUT=runs/design/scaffold_pbe0_hz3.json OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 \
  python3 scripts/177_scaffold_pbe0.py Hz-S-Hz >> runs/design/round2.log 2>&1
TD_NEUTRAL_ONLY=1 OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/166_design_tddft.py Hz-S-Hz >> runs/design/round2.log 2>&1
TD_ANION_NROOTS=12 OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/166_design_tddft.py Hz-S-Hz >> runs/design/round2.log 2>&1

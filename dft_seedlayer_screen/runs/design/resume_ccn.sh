#!/bin/bash
# after the restart: Hz-C(CN)2-Hz TD-DFT (PBE0 sites already done)
cd "$(dirname "$0")/../.."
TD_NEUTRAL_ONLY=1 OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/166_design_tddft.py "Hz-C(CN)2-Hz" >> runs/design/round2.log 2>&1
TD_ANION_NROOTS=12 OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/166_design_tddft.py "Hz-C(CN)2-Hz" >> runs/design/round2.log 2>&1

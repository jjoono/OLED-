#!/bin/bash
# electrical (low-Tel GFN2 + PBE0 four-point) for the round-2 dimers, then the cross-guess SCF check
cd "$(dirname "$0")/../.."
env OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/180_lambda_lowT.py "Hz-C(CN)2-Hz" Hz-S-Hz >> runs/design/electrical.log 2>&1
while pgrep -f scf_check.sh >/dev/null; do sleep 120; done
env OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/181_ion_scf_check.py T30 "Hz-C(CN)2-Hz" Hz-S-Hz >> runs/design/electrical.log 2>&1

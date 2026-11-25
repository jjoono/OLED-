#!/bin/bash
# cross-guess ion SCF check: HATCN (T30) now, the dimers once 180 has their single points; NPB cation (300K)
cd "$(dirname "$0")/../.."
R="OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500"
env OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/181_ion_scf_check.py T30 HATCN >> runs/design/electrical.log 2>&1
env OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/181_ion_scf_check.py 300K NPB >> runs/design/electrical.log 2>&1
while pgrep -f 180_lambda_lowT >/dev/null; do sleep 120; done
env OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/181_ion_scf_check.py T30 Hz-NCN-Hz Hz-O-Hz >> runs/design/electrical.log 2>&1

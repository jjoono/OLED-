#!/bin/bash
# NPB cation SCF check, rerun after the first attempt hit ENOSPC (orphaned DF scratch files); waits for the dimer checks to free cores
cd "$(dirname "$0")/../.."
while [ "$(pgrep -fc 181_ion_scf_check)" -ge 2 ]; do sleep 120; done
env OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/181_ion_scf_check.py 300K NPB >> runs/design/electrical.log 2>&1

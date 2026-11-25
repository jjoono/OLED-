#!/bin/bash
# remaining ion SCF checks with the scratch-freeing 181: Hz-O-Hz once at most one older 181 job is left,
# NPB (largest DF tensor, ~7 GB) once nothing else runs. Replaces scf_check.sh tail and scf_npb.sh.
cd "$(dirname "$0")/../.."
E="env OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500"
while [ "$(pgrep -fc 'python3 scripts/181_ion_scf_check')" -ge 2 ]; do sleep 120; done
$E python3 scripts/181_ion_scf_check.py T30 Hz-O-Hz >> runs/design/electrical.log 2>&1
while pgrep -f 'python3 scripts/181_ion_scf_check' >/dev/null; do sleep 120; done
$E python3 scripts/181_ion_scf_check.py 300K NPB >> runs/design/electrical.log 2>&1

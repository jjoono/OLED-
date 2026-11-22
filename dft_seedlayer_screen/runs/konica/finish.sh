#!/bin/bash
# open-site stage for the compounds the KeyError skipped, then the summary
cd "$(dirname "$0")/../.."
run() { OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/173_konica_table2.py --open "$@" >> runs/konica/open_fin_$1.log 2>&1; }
run No11 &
run No03 No09 No06 &
wait
python3 scripts/173_konica_table2.py --summary > runs/konica/summary.txt 2>&1
echo DONE >> runs/konica/summary.txt

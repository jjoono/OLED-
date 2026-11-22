#!/bin/bash
# after the three PBE0 workers: re-place crowded sites, then the summary
cd "$(dirname "$0")/../.."
while pgrep -f runs/konica/chain.sh >/dev/null; do sleep 60; done
run() { OMP_NUM_THREADS=4 QC_THREADS=4 python3 scripts/173_konica_table2.py --open "$@" >> runs/konica/open_$1.log 2>&1; }
run No09 No11 No04 No05 No01 &
run No07 No08 No10 &
run No06 No02 No03 &
wait
python3 scripts/173_konica_table2.py --summary > runs/konica/summary.txt 2>&1
echo DONE >> runs/konica/summary.txt

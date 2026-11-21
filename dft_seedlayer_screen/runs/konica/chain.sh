#!/bin/bash
# three workers, each builds, surveys and scores its own compounds
cd "$(dirname "$0")/../.."
run() { OMP_NUM_THREADS=4 QC_THREADS=4 python3 scripts/173_konica_table2.py --build --survey --pbe0 "$@" >> runs/konica/worker_$1.log 2>&1; }
run No09 No11 No04 No05 No01 &
run No07 No08 No10 &
run No06 No02 No03 &
wait
python3 scripts/173_konica_table2.py --summary > runs/konica/summary.txt 2>&1
echo DONE >> runs/konica/summary.txt

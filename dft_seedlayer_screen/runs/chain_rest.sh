#!/bin/bash
# resume after container restart: ion-pair N1s, then FNO CCSD(T)
cd "$(dirname "$0")/.."
export QC_THREADS=4 QC_MEM_MB=9000
CC_FNO=1 CC_MEM_MB=5000 python3 scripts/171_ccsdt_benchmark.py pyridine pyrazine benzene >> runs/ccsdt_benchmark.log 2>&1
CC_FNO=1 CC_MEM_MB=5000 python3 scripts/171_ccsdt_benchmark.py MeCN >> runs/ccsdt_benchmark.log 2>&1
RETRY_UNCONVERGED=1 python3 scripts/172_n1s_dscf.py >> runs/n1s_dscf.log 2>&1   # ion pair (anion retries skipped)

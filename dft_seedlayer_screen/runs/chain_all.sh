#!/bin/bash
# queue: Bphen:Cs Ag2 -> Ag_n PBE0 -> N1s Delta-SCF -> CCSD(T) -> Bphen:Cs spin checks
cd "$(dirname "$0")/.."
export QC_THREADS=4 QC_MEM_MB=9000
SKIP_SPIN_CHECK=1 python3 scripts/169_bphencs_whole.py >> runs/bphencs_whole.log 2>&1
python3 scripts/170_agn_wetting.py --pbe0 > runs/agn_wetting.log 2>&1
python3 scripts/172_n1s_dscf.py > runs/n1s_dscf.log 2>&1
QC_MEM_MB=11000 python3 scripts/171_ccsdt_benchmark.py > runs/ccsdt_benchmark.log 2>&1
python3 scripts/169_bphencs_whole.py >> runs/bphencs_whole.log 2>&1

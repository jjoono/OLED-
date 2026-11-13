#!/bin/bash
# after 169 (Bphen:Cs) and 170 (Ag_n PBE0): N1s Delta-SCF, then CCSD(T) benchmark
cd "$(dirname "$0")/.."
until [ -f runs/agn_wetting.log ] && ! pgrep -f "scripts/170_agn" >/dev/null && ! pgrep -f "scripts/169_bphencs" >/dev/null; do sleep 60; done
while pgrep -f "171_ccsdt_benchmark.py MeCN --test" >/dev/null; do sleep 60; done
QC_THREADS=4 QC_MEM_MB=9000 python3 scripts/172_n1s_dscf.py > runs/n1s_dscf.log 2>&1
QC_THREADS=4 QC_MEM_MB=11000 python3 scripts/171_ccsdt_benchmark.py > runs/ccsdt_benchmark.log 2>&1

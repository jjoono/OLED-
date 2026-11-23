#!/bin/bash
# after the running 177 job: drop the ENOSPC-failed site and score it again
cd "$(dirname "$0")/../.."
while pgrep -f "177_scaffold_pbe0.py triptycene-HAT" >/dev/null; do sleep 30; done
python3 - <<'P'
import json
p='runs/design/scaffold_pbe0.json'; d=json.load(open(p))
pb=d['triptycene-HAT']['pbe0']
for k in [k for k,v in pb.items() if 'error' in v]: del pb[k]
json.dump(d,open(p,'w'),indent=1)
P
OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=4000 python3 scripts/177_scaffold_pbe0.py triptycene-HAT >> runs/design/scaffold_pbe0.log 2>&1

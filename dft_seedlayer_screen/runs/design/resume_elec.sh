#!/bin/bash
# after the container restart: finish Hz-O-Hz (resumes its single points), then HATCN with the kick
cd "$(dirname "$0")/../.."
OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/179_electrical.py Hz-O-Hz >> runs/design/electrical.log 2>&1
python3 - <<'P'
import json; p='runs/design/electrical.json'; d=json.load(open(p))
if 'HATCN_nokick' not in d: d['HATCN_nokick']=d.pop('HATCN')
json.dump(d,open(p,'w'),indent=1)
P
OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/179_electrical.py HATCN >> runs/design/electrical.log 2>&1

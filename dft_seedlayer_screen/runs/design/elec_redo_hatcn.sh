#!/bin/bash
# HATCN again with the symmetry-breaking kick, after the running batch
cd "$(dirname "$0")/../.."
while pgrep -f "179_electrical.py HATCN Hz-NCN-Hz" >/dev/null; do sleep 60; done
python3 - <<'P'
import json; p='runs/design/electrical.json'; d=json.load(open(p)); d['HATCN_nokick']=d.pop('HATCN'); json.dump(d,open(p,'w'),indent=1)
P
OMP_NUM_THREADS=4 QC_THREADS=4 QC_MEM_MB=3500 python3 scripts/179_electrical.py HATCN >> runs/design/electrical.log 2>&1

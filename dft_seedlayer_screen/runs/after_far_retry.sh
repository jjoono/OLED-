#!/bin/bash
# when the running (old-code) retry of Ag_aza_far finishes, mark it retried and stop that process,
# so the queue moves on to CCSD(T); the ion pair runs in the last queue line
cd "$(dirname "$0")/.."
until [ "$(grep -c '^Ag_aza_far' runs/n1s_dscf.log)" -ge 2 ]; do sleep 30; done
python3 -c "
import json;p='runs/n1s_dscf.json';d=json.load(open(p));d['Ag_aza_far']['retried']=True;json.dump(d,open(p,'w'),indent=1)"
kill 23267

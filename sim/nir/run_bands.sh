#!/bin/bash
# 700/750/800/850 nm bands: record emitters + Lumtec-available emitters, reduced case set
cd /home/user/OLED-/sim/nir
export CASES_REDUCED=1
for em in DCNSPTPA DTN BF2 PtTPBP Pt794 PtazaTBP; do t0=$SECONDS; python3 emitter_opt.py $em; echo "$em done in $((SECONDS-t0)) s"; done
echo ALLBANDSDONE

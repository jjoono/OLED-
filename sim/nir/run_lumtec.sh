#!/bin/bash
cd /home/user/OLED-/sim/nir
for em in DTN BF2 PtTPBP; do t0=$SECONDS; python3 emitter_opt.py $em; echo "$em done in $((SECONDS-t0)) s"; done
echo ALLLUMTECDONE

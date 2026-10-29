#!/bin/bash
cd /home/user/OLED-/sim/front_lum
[ -f results_fab/TREF.json ] || python3 ce_fab.py TREF
for c in "F2 4" "F1 4" "F2 2" "F2 6" "F1 2"; do set -- $c
  [ -f results_fab/${1}_N${2}.json ] || python3 ce_fab.py $1 $2
done
echo ALLFABDONE

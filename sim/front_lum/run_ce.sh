#!/bin/bash
cd /home/user/OLED-/sim/front_lum
[ -f results_ce/REF.json ] || python3 ce_opt.py REF
for c in "0 0.01" "2 0.01" "4 0.01" "6 0.01" "4 0" "6 0"; do set -- $c
  [ -f results_ce/B_N${1}_mua${2}.json ] || python3 ce_opt.py B $1 $2
done
echo ALLCEDONE

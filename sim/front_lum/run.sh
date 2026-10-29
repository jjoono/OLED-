#!/bin/bash
cd /home/user/OLED-/sim/front_lum
for c in "REF 0" "A 0" "A 2" "A 4" "A 6" "B 0" "B 2" "B 4" "B 6"; do
  set -- $c; [ -f results/${1}_N${2}.json ] || python3 front_lum.py $1 $2
done
echo ALLDONE

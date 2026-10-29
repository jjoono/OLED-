#!/bin/bash
cd /home/user/OLED-/sim/inverted
for c in INV NORM; do for e in planar mla15 mla177; do
  [ -f results/${c}_${e}.json ] || python3 inv_opt.py $c $e
  [ -f results/${c}_${e}_ITO150.json ] || python3 inv_opt.py $c $e fix
done; done
echo ALLINVDONE

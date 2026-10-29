#!/bin/bash
cd /home/user/OLED-/sim/front_lum
until grep -q ALLFABDONE run_fab.log; do sleep 60; done
for N in 4 6; do [ -f results_fab/F1_N${N}.json ] || python3 ce_fab.py F1 $N; done
echo ALLFAB2DONE

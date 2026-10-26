#!/bin/bash
cd /home/user/OLED-/sim/nir
for te in IOH Ag ITO IZO Au; do for n in 1.8 1.5; do for nd in 0 4; do
  f=results/${te}_n${n}_dbr${nd}_q0.3_mla.json
  [ -f $f ] || python3 nir_opt.py $te $n $nd 0.3 mla
done; done; done
echo ALLDONE

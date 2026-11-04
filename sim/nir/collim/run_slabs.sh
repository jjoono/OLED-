#!/bin/bash
cd /home/user/OLED-/sim/nir/collim; mkdir -p slabs
for S in 0.3 1 3 10; do for g in 0.0 0.6 0.9; do
  f=slabs/S${S}_g${g}.npz; [ -f $f ] || python3 slab_mc.py $S $g 20000 $f
done; done
echo SLABSDONE

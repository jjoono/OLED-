#!/bin/bash
cd /home/user/OLED-/sim/nir
for em in Pt-fprpz2 PyCNTPA-24 PyCNTPA-12 Ir-tmd; do python3 emitter_opt.py $em; done
echo ALLEMDONE

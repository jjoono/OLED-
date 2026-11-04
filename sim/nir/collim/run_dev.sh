#!/bin/bash
cd /home/user/OLED-/sim/nir/collim
for c in IZO_n1.8_dbr0_mla IZO_n1.8_dbr4_mla ITO_n1.5_dbr0_flat; do
  [ -f dev/Pt_$c.npz ] || python3 device_data.py ../results_emitter/Pt-fprpz2_$c.json 0.78 dev/Pt_$c.npz
done
echo DEVDONE

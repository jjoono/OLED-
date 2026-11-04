#!/bin/bash
cd /home/user/OLED-/sim/nir/collim; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
D=dev/Pt_IZO_n1.8_dbr0_mla.npz
python3 opt_collim.py $D 740 80 0.80 4 Leq10   Pt_N4_Leq10 > logs/Pt_N4_Leq10.log 2>&1 &
python3 opt_collim.py $D 740 80 0.80 8 Leq10   Pt_N8_Leq10 > logs/Pt_N8_Leq10.log 2>&1 &
python3 opt_collim.py $D 740 80 0.80 8 E_z3_a1 Pt_N8_E3    > logs/Pt_N8_E3.log 2>&1 &
python3 opt_collim.py $D 740 80 0.80 8 C20     Pt_N8_C20   > logs/Pt_N8_C20.log 2>&1 &
wait
python3 planar_cavity.py IZO 1.5 4 Leq10 740 80 0.80 0.78 Pt_planar_IZO_n1.5_N4_Leq10 > logs/planar.log 2>&1 &
python3 opt_collim.py $D 740 20 0.80 8 Leq10 FW20_N8_Leq10 > logs/FW20.log 2>&1 &
python3 opt_collim.py $D 740 40 0.80 8 Leq10 FW40_N8_Leq10 > logs/FW40.log 2>&1 &
wait; echo BATCH3DONE

import sys, numpy as np
sys.path.insert(0,'/home/user/OLED-/sim/mla'); import raytrace as rt
BT,BR=rt.bsdf(1.0,1.75,N=20000,seed=175,shape='cap'); np.savez('mla_hemi_1.75.npz',BT=BT,BR=BR); print('done')

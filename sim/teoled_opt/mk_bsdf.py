import sys, numpy as np
sys.path.insert(0,'/home/user/OLED-/sim/mla'); import raytrace as rt
BT,BR=rt.bsdf(1.0,1.65,N=20000,seed=165,shape='cap'); np.savez('mla_hemi_1.65.npz',BT=BT,BR=BR); print('done')

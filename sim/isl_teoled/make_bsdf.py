import sys, numpy as np
sys.path.insert(0,'/home/user/OLED-/sim/mla'); import raytrace as rt
from multiprocessing import Pool
def f(n): BT,BR=rt.bsdf(1.0,n,N=20000,seed=int(n*100),shape='cap'); return n,BT,BR
if __name__=='__main__':
    with Pool(3) as p: res=p.map(f,[1.5,1.6,1.7])
    np.savez('mla_hemi_bsdf.npz',**{f'BT_{n}':BT for n,BT,BR in res},**{f'BR_{n}':BR for n,BT,BR in res})
    print('done')

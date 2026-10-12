import numpy as np, sys
from scipy.io import loadmat
from scipy.interpolate import PchipInterpolator
import hftrace as T
R='/home/user/OLED-/'
C=loadmat(R+'opt_4band_result_25by25.mat'); H=loadmat(R+'opt_hemisphere_result.mat'); F=loadmat(R+'freeform_EQEtotal_result.mat')
HX,HY=H['HEMI_X'].ravel(),H['HEMI_Y'].ravel()
def prof(xv):
    px=np.concatenate(([0.0],np.sort(xv[0:5]),[1.0])); py=np.concatenate(([1.0],xv[5:10],[0.0]))
    px,u=np.unique(px,return_index=True); py=py[u]; t=np.linspace(0,1,20001)
    return t, PchipInterpolator(px,py)(t)*xv[12]
names=['0-20','20-40','40-60','60-80','total']
designs={}
for j in range(5):
    xv=np.concatenate((HX,HY,H['hemi_x'][j])); designs['hemi_'+names[j]]=(xv,)
for j in range(4): designs['free_'+names[j]]=(C['band_x'][j],)
designs['free_total']=(F['gBestX'].ravel(),)
out={}
for k,(xv,) in designs.items():
    rho,h=prof(xv)
    BT,BR,lost=T.bsdf(rho,h,n=1.51,N=20000)
    out[k+'_BT']=BT; out[k+'_BR']=BR; out[k+'_x']=xv; out[k+'_lost']=lost
    print(k,'p=%.4f'%((BT.sum(1)*np.cos(np.radians(np.arange(90)+.5))*np.sin(np.radians(np.arange(90)+.5))).sum()/ (np.cos(np.radians(np.arange(90)+.5))*np.sin(np.radians(np.arange(90)+.5))).sum()),'lost max %.4f'%lost.max(),flush=True)
    np.savez('bsdfs.npz',**out)

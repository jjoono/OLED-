import numpy as np, matplotlib
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

matplotlib.use('Agg')
import matplotlib.pyplot as plt
d=np.load(str(OUT_DIR / 'izo_data.npz'))
f=np.load(str(OUT_DIR / 'izo_fit_O2.npz'))
wl=d['o2wl']; P=d['o2psi']; L=d['o2del']; dep=d['o2dep']
Pce=d['o2psi_ce']; Lce=d['o2del_ce']
Pm=f['psi_c']; Lm=f['del_c']
ang=[45,50,55,60,65]
fig,ax=plt.subplots(3,2,figsize=(15,12))
cmap=matplotlib.colormaps['tab10']
for i in range(5):
    c=cmap(i/9)
    ax[0,0].plot(wl,P[:,i],'.',ms=2,color=c)
    ax[0,0].plot(wl,Pm[:,i],'-',lw=0.9,color=c)
    ax[0,1].plot(wl,L[:,i],'.',ms=2,color=c,label=f'{ang[i]}°')
    ax[0,1].plot(wl,Lm[:,i],'-',lw=0.9,color=c)
    ax[1,0].plot(wl,P[:,i],'.',ms=2,color=c)
    ax[1,0].plot(wl,Pce[:,i],'-',lw=0.9,color=c)
    ax[1,1].plot(wl,L[:,i],'.',ms=2,color=c)
    ax[1,1].plot(wl,Lce[:,i],'-',lw=0.9,color=c)
    ax[2,0].plot(wl,dep[:,i],'-',lw=0.8,color=c)
ax[0,0].set_title('O2 data(dots) vs MY final model(lines): Psi')
ax[0,1].set_title('Delta'); ax[0,1].legend(fontsize=7)
ax[1,0].set_title('O2 data vs CompleteEASE model curves: Psi')
ax[1,1].set_title('Delta')
ax[2,0].set_title('measured %depol'); ax[2,0].set_ylim(-5,25)
for a in ax.ravel(): a.grid(alpha=0.2); a.set_xlabel('wl (nm)')
dd_my=(Lm-L+180)%360-180; dd_ce=(Lce-L+180)%360-180
txt=''
for lo,hi in [(245,350),(350,450),(450,600),(600,800),(800,1000),(1000,1200),(1200,1450),(1450,1690)]:
    m=(wl>=lo)&(wl<hi)
    txt+='%4d-%4d: MY psi %5.2f del %6.2f | CE psi %5.2f del %6.2f\n'%(lo,hi,
        np.sqrt(np.nanmean((Pm-P)[m]**2)),np.sqrt(np.nanmean(dd_my[m]**2)),
        np.sqrt(np.nanmean((Pce-P)[m]**2)),np.sqrt(np.nanmean(dd_ce[m]**2)))
ax[2,1].axis('off'); ax[2,1].text(0.02,0.95,txt,family='monospace',fontsize=9,va='top')
plt.tight_layout()
plt.savefig(str(OUT_DIR / 'izo_O2_look.png'),dpi=140)
print(txt)
print('saved izo_O2_look.png')

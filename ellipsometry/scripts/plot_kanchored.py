import numpy as np, matplotlib
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ellipsometry_fit as ef

d=np.load(str(OUT_DIR / 'Perov_kanchored_fit.npz'))
wl=d['wl']; pm=d['psi_m']; dm=d['del_m']; pc=d['psi_c']; dc=d['del_c']
n=d['n']; k=d['k']; W=d['W']; wl_si=d['wl_si']; k_si=d['k_si']
angles=[65,70,75]
dc_p,dm_p=ef._align_delta(dc,dm)

fig,ax=plt.subplots(2,2,figsize=(14,9))
fig.suptitle('Perovskite — UV-Vis k-anchored SE fit (final)\n'
             f'k(510-570nm) pinned to SI Fig.6a (UV-Vis) | d fixed by SEM: bulk={float(d["d_bulk"]):.0f}nm '
             f'+ rough {float(d["d_rough"]):.0f}nm(f={float(d["f_rough"]):.2f}) = 270nm | KK-consistent Gen-Osc',
             fontsize=9,fontweight='bold')
cmap=matplotlib.colormaps['tab10']; cols=[cmap(i/9) for i in range(3)]
axp,axd,axn,axk=ax[0,0],ax[0,1],ax[1,0],ax[1,1]
for i,(a,c) in enumerate(zip(angles,cols)):
    hi=W[:,i]>0.5; lo=~hi
    axp.plot(wl[lo],pm[lo,i],'.',ms=2,color='0.85')
    axp.plot(wl[hi],pm[hi,i],'o',ms=1.6,color=c,alpha=0.35)
    axp.plot(wl,pc[:,i],'-',lw=1.1,color=c,label=f'{a}°')
    axd.plot(wl[lo],dm_p[lo,i],'.',ms=2,color='0.85')
    axd.plot(wl[hi],dm_p[hi,i],'o',ms=1.6,color=c,alpha=0.35)
    axd.plot(wl,dc_p[:,i],'-',lw=1.1,color=c,label=f'{a}°')
axp.set_title('Ψ (grey = depol-excluded)'); axp.set_ylabel('Ψ (deg)')
axd.set_title('Δ'); axd.set_ylabel('Δ (deg)')
for a_ in (axp,axd): a_.set_xlabel('Wavelength (nm)'); a_.legend(fontsize=7); a_.grid(True,alpha=0.2)

axn.plot(wl,n,'b-',lw=1.6,label='n (SE, k-anchored)')
axn.axhline(2.3,color='gray',ls='--',lw=0.9,label='paper assumption n=2.3')
axn.axvspan(510,560,color='green',alpha=0.08,label='emission band')
axn.set_title('Perovskite n(λ)'); axn.set_xlabel('Wavelength (nm)'); axn.set_ylabel('n')
axn.set_xlim(380,1100); axn.legend(fontsize=8); axn.grid(True,alpha=0.2)

axk.plot(wl,k,'r-',lw=1.6,label='k (model, KK-consistent)')
axk.plot(wl_si,k_si,'ko',ms=4,mfc='none',label='UV-Vis (SI Fig.6a, digitized)')
axk.set_title('Perovskite k(λ) — anchored to UV-Vis')
axk.set_xlabel('Wavelength (nm)'); axk.set_ylabel('k')
axk.set_xlim(380,700); axk.set_ylim(bottom=0); axk.legend(fontsize=8); axk.grid(True,alpha=0.2)

plt.tight_layout()
plt.savefig(str(OUT_DIR / 'Perov_kanchored_fit.png'),dpi=150)
print('Figure -> Perov_kanchored_fit.png')

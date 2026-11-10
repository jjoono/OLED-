import numpy as np, matplotlib
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ellipsometry_fit as ef

d=np.load(str(OUT_DIR / 'Perov_reliable_fit.npz'))
wl=d['wl']; pm=d['psi_m']; dm=d['del_m']; pc=d['psi_c']; dc=d['del_c']
n=d['n']; k=d['k']; nlo=d['n_lo']; nhi=d['n_hi']; klo=d['k_lo']; khi=d['k_hi']
W=d['W']; angles=[65,70,75]
dc_p,dm_p=ef._align_delta(dc,dm)

fig,ax=plt.subplots(2,2,figsize=(14,9))
fig.suptitle('Perovskite on damaged GraHIL — reliability-maximized fit\n'
             f'bulk={float(d["d_bulk"]):.0f}nm rough={float(d["d_rough"]):.0f}nm(f={float(d["f_rough"]):.2f}) '
             f'total≈270nm(SEM)  |  KK-corrected Gen-Osc, depol+angle-quality weights, '
             'band = angle-pair bootstrap + alias + SiO2 + blue-weight envelope',
             fontsize=9,fontweight='bold')
cmap=matplotlib.colormaps['tab10']; cols=[cmap(i/9) for i in range(3)]
axp,axd,axn,axk=ax[0,0],ax[0,1],ax[1,0],ax[1,1]
for i,(a,c) in enumerate(zip(angles,cols)):
    # grey-out low-weight (unreliable) points
    hi=W[:,i]>0.5; lo=~hi
    axp.plot(wl[lo],pm[lo,i],'.',ms=2,color='0.85')
    axp.plot(wl[hi],pm[hi,i],'o',ms=1.6,color=c,alpha=0.35)
    axp.plot(wl,pc[:,i],'-',lw=1.1,color=c,label=f'{a}°')
    axd.plot(wl[lo],dm_p[lo,i],'.',ms=2,color='0.85')
    axd.plot(wl[hi],dm_p[hi,i],'o',ms=1.6,color=c,alpha=0.35)
    axd.plot(wl,dc_p[:,i],'-',lw=1.1,color=c,label=f'{a}°')
axp.set_title('Ψ (grey = excluded by depol)'); axp.set_xlabel('Wavelength (nm)'); axp.set_ylabel('Ψ (deg)')
axd.set_title('Δ'); axd.set_xlabel('Wavelength (nm)'); axd.set_ylabel('Δ (deg)')
for a_ in (axp,axd): a_.legend(fontsize=7); a_.grid(True,alpha=0.2)

axn.fill_between(wl,nlo,nhi,color='b',alpha=0.15,label='uncertainty envelope')
axn.plot(wl,n,'b-',lw=1.6,label='n (best estimate)')
axn.axvline(549,color='gray',ls=':',lw=0.8)
axn.set_title('Perovskite n(λ)'); axn.set_xlabel('Wavelength (nm)'); axn.set_ylabel('n')
axn.legend(fontsize=8); axn.grid(True,alpha=0.2)

axk.fill_between(wl,klo,khi,color='r',alpha=0.15,label='uncertainty envelope')
axk.plot(wl,k,'r-',lw=1.6,label='k (best estimate)')
axk.axvline(560,color='green',ls=':',lw=0.8,label='560 nm (k→0 enforced)')
axk.set_title('Perovskite k(λ)'); axk.set_xlabel('Wavelength (nm)'); axk.set_ylabel('k')
axk.set_ylim(bottom=0); axk.legend(fontsize=8); axk.grid(True,alpha=0.2)

plt.tight_layout()
plt.savefig(str(OUT_DIR / 'Perov_reliable_fit.png'),dpi=150)
print('Figure -> Perov_reliable_fit.png')

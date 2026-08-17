import numpy as np, matplotlib, csv
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'
plt.rcParams['axes.unicode_minus']=False
import ellipsometry_fit as ef

# final nk + band
wl_c,nb,kb,nlo,nhi=[],[],[],[],[]
with open(str(OUT_DIR / 'Perov_FINAL_nk.csv')) as f:
    r=csv.reader(f); next(r)
    for row in r:
        wl_c.append(float(row[0])); nb.append(float(row[1])); kb.append(float(row[2]))
        nlo.append(float(row[3])); nhi.append(float(row[4]))
wl_c=np.array(wl_c); nb=np.array(nb); kb=np.array(kb); nlo=np.array(nlo); nhi=np.array(nhi)

d=np.load(str(OUT_DIR / 'Perov_final_window_fit.npz'))
wl=d['wl']; pm=d['psi_m']; dm=d['del_m']; pc=d['psi_c']; dc=d['del_c']; W=d['W']
wl_si=d['wl_si']; k_si=d['k_si']
angles=[65,70,75]
dc_p,dm_p=ef._align_delta(dc,dm)

fig,ax=plt.subplots(2,2,figsize=(14.5,9.5))
fig.suptitle('#2 시료 (Si / SiO$_2$ 2nm / damaged GraHIL / Perovskite) — 최종 n,k 결과\n'
             'Perovskite (FA$_{0.7}$MA$_{0.1}$GA$_{0.2}$)$_{0.87}$Cs$_{0.13}$PbBr$_3$,  두께 260 nm (SEM 창 260–280 자유),  '
             'k: UV-Vis 앵커 (Nature SI Fig.6a),  KK-consistent Gen-Osc',
             fontsize=10.5,fontweight='bold')
cmap=matplotlib.colormaps['tab10']; cols=[cmap(i/9) for i in range(3)]
axp,axd,axn,axk=ax[0,0],ax[0,1],ax[1,0],ax[1,1]

for i,(a,c) in enumerate(zip(angles,cols)):
    hi=W[:,i]>0.5; lo=~hi
    axp.plot(wl[lo],pm[lo,i],'.',ms=2,color='0.88')
    axp.plot(wl[hi],pm[hi,i],'o',ms=1.6,color=c,alpha=0.3)
    axp.plot(wl,pc[:,i],'-',lw=1.1,color=c,label=f'{a}°')
    axd.plot(wl[lo],dm_p[lo,i],'.',ms=2,color='0.88')
    axd.plot(wl[hi],dm_p[hi,i],'o',ms=1.6,color=c,alpha=0.3)
    axd.plot(wl,dc_p[:,i],'-',lw=1.1,color=c,label=f'{a}°')
axp.set_title('Ψ  (점=측정, 선=모델, 회색=depol 제외)',fontsize=10)
axd.set_title('Δ',fontsize=10)
axp.set_ylabel('Ψ (deg)'); axd.set_ylabel('Δ (deg)')
for a_ in (axp,axd):
    a_.set_xlabel('Wavelength (nm)'); a_.legend(fontsize=7); a_.grid(True,alpha=0.2)

m=(wl_c>=380)&(wl_c<=1100)
axn.fill_between(wl_c[m],nlo[m],nhi[m],color='tab:blue',alpha=0.18,label='불확실성 밴드 (두께창+청색k+통계)')
axn.plot(wl_c[m],nb[m],color='tab:blue',lw=1.8,label='n (best estimate)')
axn.axhline(2.3,color='gray',ls='--',lw=0.9,label='논문 가정 n=2.3')
axn.axvspan(510,560,color='green',alpha=0.07)
axn.annotate('발광대역\nn=2.35–2.47 (±3%)',xy=(535,2.44),xytext=(640,2.52),fontsize=9,
             arrowprops=dict(arrowstyle='->',lw=0.9))
axn.annotate('n(600)=2.28 ±2.8%',xy=(600,2.275),xytext=(720,2.33),fontsize=9,
             arrowprops=dict(arrowstyle='->',lw=0.9))
axn.set_title('Perovskite  n(λ)',fontsize=10)
axn.set_xlabel('Wavelength (nm)'); axn.set_ylabel('n')
axn.set_xlim(380,1100); axn.legend(fontsize=8,loc='lower right'); axn.grid(True,alpha=0.2)

mk=(wl_c>=380)&(wl_c<=700)
axk.plot(wl_c[mk],kb[mk],color='tab:red',lw=1.8,label='k (KK-consistent 모델)')
axk.plot(wl_si,k_si,'ko',ms=4.5,mfc='none',mew=1.1,label='UV-Vis 실측 (SI Fig.6a 디지타이즈)')
axk.axvspan(510,560,color='green',alpha=0.07)
axk.annotate('엑시톤 k=0.132 @524nm',xy=(524,0.135),xytext=(560,0.25),fontsize=9,
             arrowprops=dict(arrowstyle='->',lw=0.9))
axk.annotate('k=0 (λ>555nm)',xy=(575,0.004),xytext=(600,0.08),fontsize=9,
             arrowprops=dict(arrowstyle='->',lw=0.9))
axk.set_title('Perovskite  k(λ) — UV-Vis 앵커',fontsize=10)
axk.set_xlabel('Wavelength (nm)'); axk.set_ylabel('k')
axk.set_xlim(380,700); axk.set_ylim(bottom=0); axk.legend(fontsize=8); axk.grid(True,alpha=0.2)

plt.tight_layout(rect=[0,0,1,0.94])
plt.savefig(str(OUT_DIR / 'Perov_FINAL_summary.png'),dpi=150)
print('Figure -> Perov_FINAL_summary.png')

import numpy as np, matplotlib, csv
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False

d=np.load(str(OUT_DIR / 'izo_data.npz'))
ch=np.load(str(OUT_DIR / 'izo_chain_O2.npz'))
two=np.load(str(OUT_DIR / 'izo_two_O2.npz'))
wl=d['o2wl']; P=d['o2psi']; L=d['o2del']
Pm=two['psi_c']; Lm=two['del_c']
ang=[45,50,55,60,65]
cmap=matplotlib.colormaps['tab10']

fig,ax=plt.subplots(2,2,figsize=(14.5,10))
for i in range(5):
    c=cmap(i/9)
    ax[0,0].plot(wl,P[:,i],'.',ms=2,color=c)
    ax[0,0].plot(wl,Pm[:,i],'-',lw=0.9,color=c,label=f'{ang[i]}°' if i%2==0 else None)
    ax[0,1].plot(wl,L[:,i],'.',ms=2,color=c)
    ax[0,1].plot(wl,Lm[:,i],'-',lw=0.9,color=c)
ax[0,0].set_title('O₂_IZO3: 데이터(점) vs best 전역모델(2층, 선) — Ψ')
ax[0,1].set_title('Δ')
ax[0,0].legend(fontsize=8)
for a in ax[0]: a.set_xlabel('wavelength (nm)'); a.grid(alpha=0.2)

wlc=ch['wl']; nc=ch['n']; kc=ch['k']; rc=ch['res']
good=rc<0.08
ax[1,0].plot(wlc,nc,color='0.75',lw=1)
ax[1,0].plot(wlc[good],nc[good],'b.',ms=3,label='n (모델-프리 추출, 신뢰점)')
ax[1,0].plot(two['wl'],two['n_bot'],'g--',lw=1.2,label='n (2층모델 bottom)')
ax[1,0].plot(two['wl'],two['n_top'],'g:',lw=1.2,label='n (2층모델 top)')
ax[1,0].axhline(1,color='gray',lw=0.5)
ax[1,0].annotate('플라즈마 엣지\n(자유전자, 저감쇠)',xy=(1600,1.18),xytext=(1150,1.45),fontsize=9,
                 arrowprops=dict(arrowstyle='->',lw=0.9))
ax[1,0].set_title('IZO n(λ) — O₂ 200°C 1h'); ax[1,0].set_ylim(0.8,2.6)
ax[1,0].legend(fontsize=8); ax[1,0].grid(alpha=0.2); ax[1,0].set_xlabel('wavelength (nm)')

ax[1,1].plot(wlc,kc,color='0.75',lw=1)
ax[1,1].plot(wlc[good],kc[good],'r.',ms=3,label='k (모델-프리 추출, 신뢰점)')
ax[1,1].set_title('IZO k(λ)'); ax[1,1].set_ylim(0,0.35)
ax[1,1].legend(fontsize=8); ax[1,1].grid(alpha=0.2); ax[1,1].set_xlabel('wavelength (nm)')

fig.suptitle('O₂ 1h 200°C IZO3 (#2.xlsx) — 두께 d≈183–191nm, NIR 플라즈마 엣지 보유\n'
             '모델-프리 추출(d=191nm 고정): 가시광 n≈1.9, k≈0 · NIR에서 n 급락(고이동도 캐리어) · '
             '균일 단층 B-spline(CE) 실패 원인',fontsize=11,fontweight='bold')
plt.tight_layout(rect=[0,0,1,0.93])
plt.savefig(str(OUT_DIR / 'IZO_O2_final.png'),dpi=140)
print('fig saved')

with open(str(OUT_DIR / 'IZO_O2_nk.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['wl_nm','n','k','ncs_residual','flag(res<0.08 reliable)'])
    for a,b,c,r in zip(wlc,nc,kc,rc):
        w.writerow(['%.2f'%a,'%.4f'%b,'%.4f'%c,'%.4f'%r,int(r<0.08)])
print('csv saved')

import numpy as np, matplotlib, csv
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family']='Malgun Gothic'; plt.rcParams['axes.unicode_minus']=False
f=np.load(str(OUT_DIR / 'ag_final.npz'))
wl=f['wl']
# Ag25: fixed-25 branch (user prior, literature-consistent)
n25,k25=f['n25b'],f['k25b']
n5,k5=f['n5'],f['k5']
n3,k3=f['n3'],f['k3']

fig,ax=plt.subplots(2,2,figsize=(14,10))
# Palik Ag reference points
palik_wl=[350,400,500,600,700,800,1000]
palik_n=[0.07,0.05,0.05,0.06,0.08,0.09,0.13]
palik_k=[1.42,2.11,2.87,3.75,4.52,5.20,6.83]

ax[0,0].plot(wl,n25,'-',color='tab:blue',lw=1.8,label='Ag25 (연속막, d=25nm 고정)')
ax[0,0].plot(wl,n5,'-',color='tab:green',lw=1.8,label='Ag5 (연속막, d≈3.1nm)')
ax[0,0].plot(wl,n3,'-',color='tab:red',lw=1.8,label='Ag3 (섬층 유효매질, d_eff≈13nm)')
ax[0,0].plot(palik_wl,palik_n,'k^',ms=6,mfc='none',label='벌크 Ag (Palik)')
ax[0,0].set_title('n(λ)'); ax[0,0].set_ylim(0,4); ax[0,0].legend(fontsize=8)
ax[0,1].plot(wl,k25,'-',color='tab:blue',lw=1.8)
ax[0,1].plot(wl,k5,'-',color='tab:green',lw=1.8)
ax[0,1].plot(wl,k3,'-',color='tab:red',lw=1.8)
ax[0,1].plot(palik_wl,palik_k,'k^',ms=6,mfc='none')
ax[0,1].set_title('k(λ)')
e1_25=n25**2-k25**2; e1_5=n5**2-k5**2; e1_3=n3**2-k3**2
ax[1,0].plot(wl,e1_25,color='tab:blue',lw=1.8)
ax[1,0].plot(wl,e1_5,color='tab:green',lw=1.8)
ax[1,0].plot(wl,e1_3,color='tab:red',lw=1.8)
ax[1,0].axhline(0,color='k',lw=0.7)
ax[1,0].set_title('ε₁(λ) — 음수=금속성'); ax[1,0].set_ylim(-60,10)
for a in [ax[0,0],ax[0,1],ax[1,0]]: a.set_xlabel('wavelength (nm)'); a.grid(alpha=0.25)
txt='''HATCN 실측 두께 (bare, stage1):
  HATCN30: 21.0 nm (명목 70%)
  HATCN2/4/6: 11.0/11.8/12.7 nm  ⚠ 명목과 큰 괴리
   → 고정 오프셋 ~10nm + 0.35×명목 (증착 오버슈트 의심)

Ag 두께 (전역피팅, HATCN2/4/6/30 기판 순):
  Ag25: d=25nm 고정(사용자 prior; 자유피팅시 15-22
        — 불투명 금속의 d-ε 퇴화) rms=0.064
  Ag5:  3.06 / 3.08 / 3.13 / 3.83 nm  rms=0.040
        → 명목 5의 ~62-77%, 금속성 연속막!
  Ag3:  12.8 / 12.5 / 12.9 / 15 nm(유효높이)
        rough 2-4nm, rms=0.033
        → percolation 미만 섬층 (ε₁>0, 유전체적)

핵심: percolation 문턱이 명목 3↔5nm 사이.
HATCN seed 위에서 실두께 ~3nm에 연속 금속막 달성
(산화물 위 통상 6-10nm 대비 우수) — Ag-CN 배위
핵생성과 부합.'''
ax[1,1].axis('off'); ax[1,1].text(0.02,0.98,txt,family='Malgun Gothic',fontsize=9.5,va='top')
fig.suptitle('Thin Ag on HATCN — 최종 nk·두께 (5각도 NCS 전역피팅, 4기판 공유 ε)',fontsize=12,fontweight='bold')
plt.tight_layout(rect=[0,0,1,0.95])
plt.savefig(str(OUT_DIR / 'AgHATCN_final.png'),dpi=140)
with open(str(OUT_DIR / 'AgHATCN_nk.csv'),'w',newline='') as fo:
    w=csv.writer(fo)
    w.writerow(['wl_nm','n_Ag25','k_Ag25','n_Ag5','k_Ag5','n_Ag3eff','k_Ag3eff'])
    for i in range(len(wl)):
        w.writerow(['%.2f'%wl[i],'%.4f'%n25[i],'%.4f'%k25[i],'%.4f'%n5[i],'%.4f'%k5[i],'%.4f'%n3[i],'%.4f'%k3[i]])
print('saved AgHATCN_final.png / AgHATCN_nk.csv')

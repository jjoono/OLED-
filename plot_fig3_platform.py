"""
plot_fig3_platform.py  --  Fig 3 후보 그림 생성

(a) 각도별 투과 프로파일: flat(class A) vs 최적 DBR vs ideal filter, 목표밴드 표시
(b) 손실 a 의존성: flat / ideal 의 밴드 비율 (플랫폼 이득이 손실에 어떻게 눌리는지)
(c) 대역폭 의존성: DBR 이 벽을 넘으려면 얼마나 단색이어야 하는가  <-- 핵심 설계조건
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

N_SUB = 1.51
TH_LO, TH_HI = 40., 60.
TH_SUB_LO = np.degrees(np.arcsin(np.sin(np.radians(TH_LO))/N_SUB))
TH_SUB_HI = np.degrees(np.arcsin(np.sin(np.radians(TH_HI))/N_SUB))
WALL = np.sin(np.radians(TH_HI))**2 - np.sin(np.radians(TH_LO))**2

d1 = np.load('angular_recycling_result.npz')
d3 = np.load('angular_recycling_bandwidth.npz')

fig, ax = plt.subplots(1, 3, figsize=(13.5, 3.9))

# ---- (a) 각도 프로파일 ----
th = d1['th']
ax[0].axvspan(TH_SUB_LO, TH_SUB_HI, color='0.85', label='target band')
ax[0].plot(th, d1['T_flat'], lw=2, label='flat interface (class A)')
ax[0].plot(th, d1['T_dbr'], lw=2, label='optimized DBR (500-600 nm)')
Ti = np.zeros_like(th); Ti[(th >= TH_SUB_LO) & (th <= TH_SUB_HI)] = 1
ax[0].plot(th, Ti, 'k--', lw=1.3, label='ideal angular filter')
ax[0].set_xlim(0, 50); ax[0].set_ylim(0, 1.05)
ax[0].set_xlabel(r'$\theta$ in substrate (deg)'); ax[0].set_ylabel('transmittance')
ax[0].set_title('(a) single-pass angular response')
ax[0].legend(fontsize=7.5, loc='upper right'); ax[0].grid(alpha=.3)

# ---- (b) 손실 의존성 ----
a_list = np.array([0.0, 0.01, 0.02, 0.05, 0.10, 0.20])
flat_b = np.array([33.8, 33.3, 32.8, 31.3, 29.1, 25.5])
ideal_b = np.array([100.0, 94.2, 89.1, 76.6, 62.1, 45.0])
ax[1].plot(a_list*100, ideal_b, 'o-', lw=2, label='ideal angular filter')
ax[1].plot(a_list*100, flat_b, 's-', lw=2, label='flat / MLA / scatterer')
ax[1].axhline(WALL*100, color='r', ls='--', lw=1.3, label='class-A wall (33.7%)')
ax[1].set_xlabel('mirror loss per bounce  a (%)')
ax[1].set_ylabel('into target band (% of generated)')
ax[1].set_title('(b) loss dependence'); ax[1].legend(fontsize=7.5); ax[1].grid(alpha=.3)

# ---- (c) 대역폭 의존성 (핵심) ----
dl, bd, sel = d3['dlam'], np.array(d3['band'])*100, np.array(d3['sel'])*100
ax[2].plot(dl, sel, 'o-', lw=2, color='C3', label='DBR selectivity')
ax[2].axhline(WALL*100, color='r', ls='--', lw=1.3, label='class-A wall (33.7%)')
ax[2].axhline(100, color='k', ls=':', lw=1, label='ideal filter')
ax[2].set_xlabel(r'source bandwidth $\Delta\lambda$ (nm)')
ax[2].set_ylabel('band / total escaped (%)')
ax[2].set_title('(c) how monochromatic must the source be?')
ax[2].legend(fontsize=7.5); ax[2].grid(alpha=.3); ax[2].set_ylim(0, 105)

plt.tight_layout()
plt.savefig('fig3_platform.png', dpi=165)
print('saved -> fig3_platform.png')
print(f'bandwidth sweep: dlam={list(dl)}  sel(%)={list(np.round(sel,1))}')

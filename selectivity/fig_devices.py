"""Fig. S11: device inputs at 550 nm -- substrate-launched distribution v0(theta) and angle-resolved reflectance R(theta)."""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from openpyxl import Workbook
import sys; sys.path.insert(0, '..')
import figstyle as FS
FS.apply()
DV = np.load('data/devices.npz'); th = np.arange(180) * 0.5 + 0.25
W = np.diff(np.sin(np.radians(np.linspace(0, 90, 181))) ** 2)
COL = {'Al': '#4d4d4d', 'Ag': '#2a78d6', 'LL15': '#1baf7a'}; LAB = {'Al': 'Al', 'Ag': 'Ag', 'LL15': 'Low loss'}
fig, ax = plt.subplots(1, 2, figsize=(FS.W, 2.35), gridspec_kw=dict(wspace=0.36)); wb = Workbook(); ws = wb.active; ws.title = 'devices_550nm'
ws.append(['theta_sub_deg'] + [f'{d}_v0_radiance' for d in COL] + [f'{d}_R' for d in COL])
cols = []
for d, c in COL.items():
    P = DV[f'{d}_P'][75]; R = DV[f'{d}_R'][75]
    L = P / W; L = L / L[:40].mean()                     # radiance relative to its near-normal value
    ax[0].plot(th, L, color=c, label=LAB[d]); ax[1].plot(th, R, color=c, label=LAB[d]); cols.append((L, R))
    A = 1 - float(R @ W / W.sum())
    ax[1].text(3, 0.385 + 0.065 * list(COL).index(d), f"{LAB[d]}: $A'$ = {A:.3f}", color=c, fontsize=6.8)
for i in range(180):
    ws.append([float(th[i])] + [float(cols[k][0][i]) for k in range(3)] + [float(cols[k][1][i]) for k in range(3)])
ax[0].axvline(np.degrees(np.arcsin(1 / 1.51)), color='#8a8a8a', lw=0.7, ls=':'); ax[1].axvline(np.degrees(np.arcsin(1 / 1.51)), color='#8a8a8a', lw=0.7, ls=':')
ax[0].set_xlabel(r'Substrate angle $\theta$ (°)'); ax[0].set_ylabel('First-pass radiance (normalized)'); FS.label(ax[0], 'a'); ax[0].legend(fontsize=6.8, loc='upper left')
ax[1].set_xlabel(r'Substrate angle $\theta$ (°)'); ax[1].set_ylabel(r'Device reflectance $R(\theta)$'); FS.label(ax[1], 'b'); ax[1].set_ylim(0.35, 1.0)
for a in ax: a.set_xlim(0, 90); a.set_xticks([0, 30, 60, 90])
FS.save(fig, '../figures_oe/figS11_device_inputs')
wb.save('../figures_oe/figS11_device_inputs.xlsx')

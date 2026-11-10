"""Fig. S11: device inputs at 550 nm -- substrate-launched distribution v0(theta) and angle-resolved reflectance R(theta)."""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from openpyxl import Workbook
plt.rcParams.update({'font.size': 7, 'axes.spines.top': False, 'axes.spines.right': False})
DV = np.load('data/devices.npz'); th = np.arange(180) * 0.5 + 0.25
W = np.diff(np.sin(np.radians(np.linspace(0, 90, 181))) ** 2)
COL = {'Al': '#52514e', 'Ag': '#2a78d6', 'LL15': '#1baf7a'}; LAB = {'Al': 'Al', 'Ag': 'Ag', 'LL15': 'low loss'}
fig, ax = plt.subplots(1, 2, figsize=(5.6, 2.2), gridspec_kw=dict(wspace=0.35)); wb = Workbook(); ws = wb.active; ws.title = 'devices_550nm'
ws.append(['theta_sub_deg'] + [f'{d}_v0_radiance' for d in COL] + [f'{d}_R' for d in COL])
cols = []
for d, c in COL.items():
    P = DV[f'{d}_P'][75]; R = DV[f'{d}_R'][75]
    L = P / W; L = L / L[:40].mean()                     # radiance relative to its near-normal value
    ax[0].plot(th, L, color=c, label=LAB[d]); ax[1].plot(th, R, color=c, label=LAB[d]); cols.append((L, R))
    A = 1 - float(R @ W / W.sum())
    ax[1].text(2, 0.40 + 0.06 * list(COL).index(d), f"{LAB[d]}: A' = {A:.3f}", color=c, fontsize=5.5)
for i in range(180):
    ws.append([float(th[i])] + [float(cols[k][0][i]) for k in range(3)] + [float(cols[k][1][i]) for k in range(3)])
ax[0].axvline(np.degrees(np.arcsin(1 / 1.51)), color='#b5b4ad', lw=0.6, ls=':'); ax[1].axvline(np.degrees(np.arcsin(1 / 1.51)), color='#b5b4ad', lw=0.6, ls=':')
ax[0].set_xlabel('substrate angle θ (deg)'); ax[0].set_ylabel('radiance into substrate (rel. to 0–20°)'); ax[0].set_title('(a) first-pass radiance v⁽⁰⁾/w', loc='left'); ax[0].legend(frameon=False, fontsize=5.5)
ax[1].set_xlabel('substrate angle θ (deg)'); ax[1].set_ylabel('device reflectance R(θ)'); ax[1].set_title('(b) angle-resolved reflectance', loc='left'); ax[1].set_ylim(0.35, 1.0)
for a in ax: a.set_xlim(0, 90)
fig.savefig('../figures_oe/figS11_device_inputs.png', dpi=300, bbox_inches='tight'); fig.savefig('../figures_oe/figS11_device_inputs.pdf', bbox_inches='tight')
wb.save('../figures_oe/figS11_device_inputs.xlsx')

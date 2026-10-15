"""Single-pass escape probability p versus substrate index:
hemispherical MLA (LightTools BSDF, n_MLA = n_sub) vs planar interface (Fresnel) vs 1/n^2.
p = int B_T(th) cos(th) sin(th) dth / int cos(th) sin(th) dth  (Yablonovitch 1982, Eq. 7)."""
import numpy as np, scipy.io as sio, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from openpyxl import Workbook

th = np.radians(np.arange(90) + 0.5); W = np.cos(th) * np.sin(th)


thf = np.linspace(0, np.pi / 2, 200001)[1:]; Wf = np.cos(thf) * np.sin(thf)


def p_planar(n):
    """Exact planar-interface average on a fine grid (avoids 1-degree binning steps)."""
    s = n * np.sin(thf); ok = s < 1; ct = np.sqrt(np.clip(1 - s ** 2, 0, 1)); ci = np.cos(thf)
    rs = ((n * ci - ct) / (n * ci + ct)) ** 2; rp = ((n * ct - ci) / (n * ct + ci)) ** 2
    return float((np.where(ok, 1 - 0.5 * (rs + rp), 0.0) * Wf).sum() / Wf.sum())


def fresnel_T(n):
    s = n * np.sin(th); ok = s < 1; ct = np.sqrt(np.clip(1 - s ** 2, 0, 1)); ci = np.cos(th)
    rs = ((n * ci - ct) / (n * ci + ct)) ** 2; rp = ((n * ct - ci) / (n * ct + ci)) ** 2
    return np.where(ok, 1 - 0.5 * (rs + rp), 0.0)


B = sio.loadmat('lt_hemisphere_bsdf.mat')['BSDF_MLA']           # [180 out, 90 in, 15 n]
NS = np.round(np.arange(1.30, 2.001, 0.05), 2)
p_mla = np.array([(B[:90, :, k].sum(0) * W).sum() / W.sum() for k in range(NS.size)])
nn = np.linspace(1.3, 2.0, 141)
p_flat = np.array([p_planar(n) for n in nn])

fig, ax = plt.subplots(figsize=(4.2, 3.4))
ax.plot(NS, p_mla, color='0.25', ls='--', lw=1.6, label='hemispherical MLA (BSDF)')
ax.plot(nn, p_flat, color='#C6512F', lw=1.6, label='planar interface (Fresnel)')
ax.plot(nn, 1 / nn ** 2, color='#3B6EA5', ls=':', lw=1.4, label=r'$1/n_{\rm sub}^2$ (T = 1)')
ax.set_xlim(1.3, 2.0); ax.set_ylim(0, 1)
ax.set_xlabel(r'$n_{\rm sub}$'); ax.set_ylabel(r'$p$')
ax.tick_params(direction='in', top=True, right=True)
ax.legend(frameon=False, loc='upper right')
fig.savefig('p_vs_nsub.png', dpi=300, bbox_inches='tight'); fig.savefig('p_vs_nsub.pdf', bbox_inches='tight')

wb = Workbook(); ws = wb.active; ws.title = 'MLA_BSDF'
ws.append(['n_sub', 'p_MLA', 'p_planar_Fresnel', '1/n^2', 'MLA / planar'])
for n, p in zip(NS, p_mla):
    pf = p_planar(n)
    ws.append([float(n), float(p), float(pf), float(1 / n ** 2), float(p / pf)])
ws2 = wb.create_sheet('planar_continuous'); ws2.append(['n_sub', 'p_planar_Fresnel', '1/n^2'])
for n, p in zip(nn, p_flat): ws2.append([float(n), float(p), float(1 / n ** 2)])
wb.save('p_vs_nsub.xlsx')
for n in (1.3, 1.5, 1.8, 2.0):
    k = int(np.argmin(abs(NS - n))); pf = p_planar(n)
    print(f'n {n}: MLA {p_mla[k]:.3f}  planar {pf:.3f}  1/n2 {1/n**2:.3f}  MLA/planar {p_mla[k]/pf:.3f}')

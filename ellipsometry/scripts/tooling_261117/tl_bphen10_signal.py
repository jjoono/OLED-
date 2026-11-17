"""Model-free look: measured minus bare substrate, against what a continuous
Bphen film would produce.  If the measured film signal has a shape no film
thickness can scale into, no homogeneous model will ever fit it."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_ce_genosc as C, tl_genosc as GO, jnk_model as M
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
W = _os.path.join(ELLIPS_OUT, '')
CE = json.load(open(W + 'ce_genosc.json'))
pB = np.array(CE['Bphen']['with_roughness']['p']); pN = np.array(CE['NPB']['pick']['p'])

def psidel(wl, layers, ds, ang):
    N, Cc, S = M.se_ncs(wl, layers, ds, ang)
    return np.degrees(np.arccos(N) / 2), np.degrees(np.arctan2(S, Cc))

fig, ax = plt.subplots(2, 2, figsize=(12, 7), sharex=True)
for col, sh, Nf_p, lab in ((0, '2_23', pB, 'Bphen'), (1, '1_23', pN, 'NPB')):
    wl, P, D, ang = L.se(sh, 245, 1000, 1)
    air = np.ones(len(wl), complex); ox = C.N_ce('NTVE_JAW', wl); si = C.N_ce('SI_JAW', wl)
    Pb, Db = psidel(wl, [air, ox, si], [C.D_OX], ang)
    a = 4                                              # 65 deg
    ax[0, col].plot(wl, P[:, a] - Pb[:, a], 'k', lw=2, label='measured - bare')
    ax[1, col].plot(wl, ((D[:, a] - Db[:, a] + 180) % 360) - 180, 'k', lw=2, label='measured - bare')
    for d, c in ((5, '#9ecae1'), (10, '#3182bd'), (20, '#08519c')):
        Pm, Dm = psidel(wl, [air, GO.N_of(Nf_p, wl), ox, si], [d, C.D_OX], ang)
        ax[0, col].plot(wl, Pm[:, a] - Pb[:, a], color=c, lw=1, label='%s %d nm film' % (lab, d))
        ax[1, col].plot(wl, ((Dm[:, a] - Db[:, a] + 180) % 360) - 180, color=c, lw=1)
    ax[0, col].set_title('%s (%s nominal 10 nm), 65 deg' % (sh, lab))
    ax[0, col].set_ylabel('dPsi (deg)'); ax[1, col].set_ylabel('dDelta (deg)'); ax[1, col].set_xlabel('wavelength (nm)')
    for r in (0, 1): ax[r, col].grid(alpha=.3)
    ax[0, col].legend(fontsize=7)
fig.tight_layout(); fig.savefig(W + 'bphen10_signal.png', dpi=120)
print('ok')

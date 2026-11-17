"""Figure: what the Bphen 10 nm pieces look like against the candidate pictures."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_ce_genosc as C, tl_glass as GL, tl_glass_global as GG, tl_bphen10 as B, tl_genosc as GO
import jnk_model as M
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
W = _os.path.join(ELLIPS_OUT, '')
CE = json.load(open(W + 'ce_genosc.json')); p = np.array(CE['Bphen']['with_roughness']['p'])
S2 = {r['sheet'] + r['model']: r for r in json.load(open(W + 'bphen10_sub2.json'))}
GJ = json.load(open(W + 'ce_glass.json'))

def psidel(N, Cc, S):
    return np.degrees(np.arccos(N) / 2), np.degrees(np.arctan2(S, Cc)) % 360

fig, ax = plt.subplots(2, 2, figsize=(12, 7.5), sharex='col')
# --- Si 2_23
wl, P, D, ang = L.se('2_23', 245, 1000, 1); a = 4
air = np.ones(len(wl), complex); ox = C.N_ce('NTVE_JAW', wl)**2; si = C.N_ce('SI_JAW', wl)**2
eb, _ = B.eps_bphen(p, wl)
def run(eo, ds): return psidel(*B.ncs_regions(wl, ang, [(1.0, eo, eo, ds)], GG.CONV['normalised']))
v = S2['2_23rough Si (air) + Bphen']['x']
curves = [('Bphen 10 nm film (nominal)', run([air, eb, ox, si], [10.0, C.D_OX]), '#3182bd', '-'),
          ('best Bphen film, 3.6 nm (MSE 42)', run([air, eb, ox, si], [3.63, C.D_OX]), '#9ecae1', '-'),
          ('rough Si, no Bphen (MSE 14)', run([air, B.bruggeman_eps(si, air, v['fsi']), ox, si], [v['t'], C.D_OX]), '#d62728', '-')]
ax[0, 0].plot(wl, P[:, a], 'k', lw=2.2, label='measured 2_23'); ax[1, 0].plot(wl, D[:, a] % 360, 'k', lw=2.2)
for lab, (Pm, Dm), c, ls in curves:
    ax[0, 0].plot(wl, Pm[:, a], color=c, ls=ls, lw=1.2, label=lab); ax[1, 0].plot(wl, Dm[:, a], color=c, ls=ls, lw=1.2)
ax[0, 0].set_title('2_23 on Si, 65 deg'); ax[0, 0].legend(fontsize=7)
# --- glass 2_34
wl, P, D, ang = L.se('2_34', 300, 1000, 1); a = 2
A, Bc = GJ['cauchy']; ng = A + Bc / (wl / 1000)**2 + 0j; Nb = GO.N_of(p, wl)
def g(d):
    return psidel(*GG.norm(*GL.ncs_glass(wl, Nb, d, ng, ang, GJ['f'])))
ax[0, 1].plot(wl, P[:, a], 'k', lw=2.2, label='measured 2_34'); ax[1, 1].plot(wl, D[:, a] % 360, 'k', lw=2.2)
for d, c, lab in ((10.0, '#3182bd', 'Bphen 10 nm (nominal)'), (3.35, '#2ca02c', 'Bphen 3.4 nm (best film)'), (0.01, '#999', 'bare glass')):
    Pm, Dm = g(d); ax[0, 1].plot(wl, Pm[:, a], color=c, lw=1.2, label=lab); ax[1, 1].plot(wl, Dm[:, a], color=c, lw=1.2)
ax[0, 1].set_title('2_34 on glass, 55 deg'); ax[0, 1].legend(fontsize=7)
for r in range(2):
    for c in range(2): ax[r, c].grid(alpha=.3)
    ax[r, 0].set_ylabel(('Psi', 'Delta')[r] + ' (deg)')
for c in range(2): ax[1, c].set_xlabel('wavelength (nm)')
fig.tight_layout(); fig.savefig(W + 'bphen10_models.png', dpi=130)
print('ok')

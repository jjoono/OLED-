"""Kim et al. TE device: EQE ceilings when the OCR/air face is textured too (PLQY = 1).
planar | ext_recip: reciprocal complete mixing at the textured OCR/air face | both: textured face + device-side mixing |
ideal_cone_dev: Fresnel-free etendue-limited escape + device mixing | sta1: every photon in the OCR escapes (hard cap).
Usage (from selectivity/): python3 advsci_top_texture.py <material library .mat>"""
import sys, types, numpy as np
sys.modules.setdefault('run_display', types.ModuleType('run_display'))
import sel, advsci_model as AM
spec = np.exp(-4*np.log(2)*((AM.LAMS-530.0)/22.0)**2); spec /= spec.sum()
dev = AM.device(0.85, 'Ag'); NB = sel.NB; I = np.eye(NB); w = np.asarray(sel.WL, float); w /= w.sum()
th = np.asarray(sel.TH, float); nO = AM.N_OCR
acc = dict(planar=0, ext_recip=0, both=0, ideal_cone_dev=0, sta1=0)
for i, lam in enumerate(AM.LAMS):
    Rf, Tf = sel.tmm_RT([], lam, sel.TH, nO); R = dev['R'][i]; v = dev['P'][i]; es = dev['eta_sub'][i]
    def eta(p, B): return p @ np.linalg.solve(I - B, v)
    p = Tf
    # external reciprocal complete mixing at the textured OCR/air face, specular device
    q = w*(1-p); q /= q.sum(); ext = eta(p, np.diag(R) @ np.outer(q, 1-p))
    # both: textured top mixes (recip) and device mixes (recip): approximate by device-mixing D after top return
    D = np.outer(R*w, R)/(w*R).sum(); both = eta(p, D @ np.outer(q, 1-p))
    # ideal etendue-limited escape (no Fresnel loss in the cone) + device complete mixing
    pc = (nO*np.sin(np.radians(th)) < 1).astype(float); ideal = eta(pc, D @ np.diag(1-pc))
    for k, val in (('planar', eta(p, np.diag(R*(1-p)))), ('ext_recip', ext), ('both', both), ('ideal_cone_dev', ideal), ('sta1', 1.0)):
        acc[k] += spec[i]*es*val
for k, v in acc.items(): print(f'{k:15s} EQE ceiling {v:.3f}')

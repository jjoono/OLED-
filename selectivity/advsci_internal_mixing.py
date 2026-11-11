"""Ray-optical ceiling of internal mixing for the Kim et al. (Adv. Sci. 2026) top-emitting device.
Uses advsci_model.device() (TE microcavity, Ag 20 nm top, NPB 60, OCR n=1.59 exit medium).
Cases per wavelength (PLQY = 1):
  planar     : specular device, planar OCR/air face (P-O2 ceiling)
  iso        : lossless isotropic randomizer between device and OCR face (advsci_model m = 1)
  recip_dev  : reciprocal rank-one mixing by the device itself, D(i->j) = R_i R_j w_j / sum(w R)
Reports eta into OCR, WG/SPP/abs, A' seen from OCR, Cov_w(R, p), and EQE = air output / dipole power.
Usage: python3 advsci_internal_mixing.py <material library .mat>"""
import sys, json, types, numpy as np
sys.modules.setdefault('run_display', types.ModuleType('run_display'))   # only needed by advsci_model.run (not used)
import sel, advsci_model as AM

spec = np.exp(-4 * np.log(2) * ((AM.LAMS - 530.0) / 22.0) ** 2); spec /= spec.sum()
dev = AM.device(0.85, 'Ag')
NB = sel.NB; I = np.eye(NB); w = sel.WL if np.ndim(sel.WL) else None
w = np.asarray(sel.WL, float); w = w / w.sum()
res = {k: 0.0 for k in ('planar', 'iso', 'recip_dev')}; Aps, covs, pbar = [], [], []
for i, lam in enumerate(AM.LAMS):
    Rf, Tf = sel.tmm_RT([], lam, sel.TH, AM.N_OCR); p = Tf; R = dev['R'][i]; v = dev['P'][i]
    def out(B): return p @ np.linalg.solve(I - B, v)
    Bp = np.diag(R * (1 - p))                                   # specular both sides
    Ts, _ = sel.mixer(1.0); Biso = (Ts * Rf[None, :]) @ Ts      # as in advsci_model.run (m = 1)
    xi = np.linalg.solve(I - R[:, None] * Biso, v); iso = Tf @ (Ts @ xi)
    D = np.outer(R * w, R) / (w * R).sum()                      # D[j, i]
    rd = out(D @ np.diag(1 - p))
    es = dev['eta_sub'][i]
    res['planar'] += spec[i] * es * out(Bp); res['iso'] += spec[i] * es * iso; res['recip_dev'] += spec[i] * es * rd
    Aps.append(1 - (w * R).sum()); pb = (w * p).sum(); pbar.append(pb)
    covs.append(((w * (R - (w * R).sum())) * (p - pb)).sum())
wg, spp, ab = (spec @ dev['br']).tolist()
summ = dict(eta_into_OCR=float(spec @ dev['eta_sub']), WG=wg, SPP=spp, abs=ab, A_prime=float(spec @ np.array(Aps)),
            p_OCR_air=float(spec @ np.array(pbar)), cov_R_p=float(spec @ np.array(covs)),
            R_in_cone=float(spec @ np.array([(w * dev['R'][i])[:39].sum() / w[:39].sum() for i in range(len(AM.LAMS))])),
            R_out_cone=float(spec @ np.array([(w * dev['R'][i])[40:].sum() / w[40:].sum() for i in range(len(AM.LAMS))])),
            EQE=res)
print(json.dumps(summ, indent=1)); json.dump(summ, open('data/advsci_internal_mixing.json', 'w'), indent=1)

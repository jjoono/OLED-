"""Model of Kim et al., Adv. Sci. 2026 (top-emitting microcavity + gentle internal corrugation + OCR):
Al 50/Ag 50 | HAT-CN 7/NPB 75/HAT-CN 7/NPB 75/TCTA 15 | EML 30 (dipole at centre) | Bphen 42/LiF 0.6/Al 1/Ag 20 | NPB 60 | OCR (n 1.59) | air.
OCR is treated as the incoherent exit medium (1 um thick); light trapped in it by TIR at the OCR/air face returns
to the stack (reflectance R_LED seen from OCR) and the corrugation is modelled as a partial randomizer m per round trip
(m = 0: planar P-O2; m = 1: full randomization, the ray-optical ceiling of G-O2).  Light in WG/SPP modes of the thin film
(n_eff > 1.59) is not recovered (ray-optical model).  PLQY = 1, charge balance 1: values are optical ceilings.
Output: EQE, on-axis cd/A, LCF = EQE / EQE_fwd, compared with the reported P-O2 / G-O2 values."""
import os, sys, json, numpy as np, scipy.io as sio
os.environ.setdefault('OMP_NUM_THREADS', '1')
import sel, run_display as RD
from sel import cps2

LIB = sys.argv[1]
M = sio.loadmat(LIB, squeeze_me=True, struct_as_record=False)['material']
lib = lambda k: np.asarray(getattr(M, k)).astype(complex)[:301]
NK = {'Al': lib('l_Al_JO'), 'Ag': lib('l_Ag_McPeak'), 'Ag_thin': lib('l_Ag_15_onHATCN'), 'NPB': lib('l_NPB'), 'HATCN': lib('l_HATCN'),
      'TCTA': lib('l_TCTA_o_JO'), 'Bphen': lib('l_Bphen_JS'), 'LiF': lib('l_LiF'), 'EML': sel.D['n_EML_o']}
LAMS = np.arange(470.0, 601.0, 2.0)
N_OCR = 1.59


def device(theta_h, thin='Ag'):
    es, P, R, br = [], [], [], []
    for lam in LAMS:
        j = int(lam) - 400; n = lambda k: NK[k][j]; iso = lambda k, d: (n(k), n(k), d)
        above = [iso('TCTA', 15), iso('NPB', 75), iso('HATCN', 7), iso('NPB', 75), iso('HATCN', 7), iso('Ag', 50), iso('Al', 50)]
        below = [iso('Bphen', 42), iso('LiF', 0.6), iso('Al', 1), iso(thin, 20), iso('NPB', 60)]
        S = cps2.Stack(float(lam), (n('EML').real, n('EML').real), 30.0, 15.0, above=above, below=below, n_sub=N_OCR, n_top=1.5, h=theta_h)
        r = cps2.solve_pol(S, npts=3000)
        th_f = np.linspace(0, 90, 1801)[:-1] + 0.025
        pf = np.clip(np.real(cps2.sub_angular(S, th_f)), 0, None).reshape(sel.NB, 10).sum(1)
        es.append(r['air'] + r['sub']); P.append(pf / pf.sum()); br.append([r['wg'], r['spp'], r['abs']])
        R.append(np.clip(np.real(cps2.stack_reflectance(S, sel.TH)), 0, 1))
    return dict(eta_sub=np.array(es), P=np.array(P), R=np.array(R), br=np.array(br))


def run(dev, m, spec):
    Rf, Tf = np.zeros((LAMS.size, sel.NB)), np.zeros((LAMS.size, sel.NB))
    for i, lam in enumerate(LAMS):
        Rf[i], Tf[i] = sel.tmm_RT([], lam, sel.TH, N_OCR)
    Ts, Rs = sel.mixer(m)
    air = np.zeros((LAMS.size, sel.NB)); I = np.eye(sel.NB)
    for i in range(LAMS.size):
        U = Ts; BR = (Ts * Rf[i][None, :]) @ U
        x = np.linalg.solve(I - dev['R'][i][:, None] * BR, dev['P'][i]); air[i] = Tf[i] * (U @ x)
    # radiant intensity vs air angle, then luminance per current
    lo = np.degrees(np.arcsin(np.clip(N_OCR * np.sin(np.radians(sel.EDG[:-1])), 0, 1)))
    hi = np.degrees(np.arcsin(np.clip(N_OCR * np.sin(np.radians(sel.EDG[1:])), 0, 1)))
    ok = N_OCR * np.sin(np.radians(sel.EDG[:-1])) < 1
    Iang = RD.rebin(air[:, ok], lo[ok], hi[ok]) / RD.SR[None, :]
    S = spec[:, None] * dev['eta_sub'][:, None] * Iang
    yb = np.interp(LAMS, np.arange(400, 701), np.load('data/cie.npz')['y'])
    lum = 683 * (S * (yb * 1239.84 / LAMS)[:, None]).sum(0) / np.cos(np.radians(RD.AC))
    eqe = float((S * RD.SR[None, :]).sum()); k_lamb = 683 * (spec * yb * 1239.84 / LAMS).sum() / np.pi
    eqe_fwd = lum[0] / k_lamb
    return dict(EQE=eqe, cdA=float(lum[0]), EQE_fwd=float(eqe_fwd), LCF=eqe / eqe_fwd)


if __name__ == '__main__':
    spec = np.exp(-4 * np.log(2) * ((LAMS - 530.0) / 22.0) ** 2); spec /= spec.sum()     # tCzphB-Fl, R2 EL 530/22 nm
    out = {}
    for thin in ('Ag', 'Ag_thin'):
        for th in (0.67, 0.85, 0.95):
            dev = device(th, thin); A = float(spec @ (1 - (dev['R'] * sel.WL).sum(1))); es = float(spec @ dev['eta_sub'])
            wg, spp, ab = (spec @ dev['br']).tolist()
            for m in (0.0, 0.1, 0.3, 1.0):
                r = run(dev, m, spec); r.update(A_prime=A, eta_ocr=es, wg=wg, spp=spp, absorb=ab)
                out[f'{thin}|{th}|{m}'] = r
                print(f"thinAg={thin:8s} Theta {th:.2f} m {m:.1f}: EQE {r['EQE']:.3f}  cd/A {r['cdA']:5.0f}  EQE_fwd {r['EQE_fwd']:.3f}  LCF {r['LCF']:.3f} | "
                      f"into OCR {es:.3f} WG {wg:.3f} SPP {spp:.3f} abs {ab:.3f} A'(from OCR) {A:.3f}", flush=True)
    json.dump(out, open('data/advsci_model.json', 'w'), indent=1)

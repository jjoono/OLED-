"""On-axis cd/A and Lambertian-assumed EQE (EQE_L = pi x photons per sr at 0 deg) for the Kim et al. TE device.
Cases: P2 (NPB/air exit), P-O2 planar (OCR exit), G-O2 ray-optical ceilings (isotropic / reciprocal device-side
complete mixing in the OCR). PLQY = 1, charge balance 1. Usage: python3 advsci_onaxis.py <material .mat>"""
import sys, types, json, numpy as np
sys.modules.setdefault('run_display', types.ModuleType('run_display'))
import sel, advsci_model as AM
from sel import cps2

def ybar(l):   # CIE 1931 y-bar, multi-lobe fit (Wyman et al. 2013)
    g = lambda x, m, s1, s2: np.exp(-0.5 * ((x - m) / np.where(x < m, s1, s2)) ** 2)
    return 0.821 * g(l, 568.8, 46.9, 40.5) + 0.286 * g(l, 530.9, 16.3, 31.1)

def device(th_h, thin, n_exit):
    es, P, R = [], [], []
    for lam in AM.LAMS:
        j = int(lam) - 400; n = lambda k: AM.NK[k][j]; iso = lambda k, d: (n(k), n(k), d)
        above = [iso('TCTA', 15), iso('NPB', 75), iso('HATCN', 7), iso('NPB', 75), iso('HATCN', 7), iso('Ag', 50), iso('Al', 50)]
        below = [iso('Bphen', 42), iso('LiF', 0.6), iso('Al', 1), iso(thin, 20), iso('NPB', 60)]
        S = cps2.Stack(float(lam), (n('EML').real, n('EML').real), 30.0, 15.0, above=above, below=below, n_sub=n_exit, n_top=1.5, h=th_h)
        r = cps2.solve_pol(S, npts=3000)
        th_f = np.linspace(0, 90, 1801)[:-1] + 0.025
        pf = np.clip(np.real(cps2.sub_angular(S, th_f)), 0, None).reshape(sel.NB, 10).sum(1)
        es.append(r['air'] + r['sub']); P.append(pf / pf.sum())
        R.append(np.clip(np.real(cps2.stack_reflectance(S, sel.TH)), 0, 1))
    return np.array(es), np.array(P), np.array(R)

E = np.asarray(sel.EDG, float); w = np.asarray(sel.WL, float); w = w / w.sum(); I = np.eye(sel.NB)
spec = np.exp(-4 * np.log(2) * ((AM.LAMS - 530.0) / 22.0) ** 2); spec /= spec.sum()   # photon spectrum (intrinsic)
K = 683 * ybar(AM.LAMS) * 1239.84 / AM.LAMS

def onaxis(air_bins, n):  # air_bins: power per exit-medium bin (fraction of dipole power); returns photons per sr at 0 deg in air
    s = n * np.sin(np.radians(E)); ok = s[1:] < np.sin(np.radians(6.0))     # air angles < ~6 deg
    a_lo = np.degrees(np.arcsin(np.clip(s[:-1], 0, 1))); a_hi = np.degrees(np.arcsin(np.clip(s[1:], 0, 1)))
    dOm = 2 * np.pi * (np.cos(np.radians(a_lo)) - np.cos(np.radians(a_hi)))
    return air_bins[ok].sum() / dOm[ok].sum()

out = {}
for th_h in (0.85, 0.95):
    for thin in ('Ag', 'Ag_thin'):
        for case, n_exit in (('P2', 1.0), ('OCR', AM.N_OCR)):
            es, P, R = device(th_h, thin, n_exit)
            res = {}
            modes = ['planar'] if case == 'P2' else ['planar', 'iso', 'recip']
            for mode in modes:
                f0 = np.zeros(len(AM.LAMS)); eq = np.zeros(len(AM.LAMS))
                for i, lam in enumerate(AM.LAMS):
                    if case == 'P2':
                        air = es[i] * P[i]
                    else:
                        Rf, Tf = sel.tmm_RT([], lam, sel.TH, n_exit); p = Tf
                        if mode == 'planar':
                            x = np.linalg.solve(I - np.diag(R[i] * (1 - p)), P[i]); air = es[i] * p * x
                        elif mode == 'iso':
                            Ts, _ = sel.mixer(1.0); x = np.linalg.solve(I - R[i][:, None] * ((Ts * Rf[None, :]) @ Ts), P[i]); air = es[i] * Tf * (Ts @ x)
                        else:
                            D = np.outer(R[i] * w, R[i]) / (w * R[i]).sum(); x = np.linalg.solve(I - D @ np.diag(1 - p), P[i]); air = es[i] * p * x
                    f0[i] = onaxis(air, n_exit); eq[i] = air.sum()
                EQE = float(spec @ eq); EQEL = float(np.pi * (spec @ f0)); CE = float((spec * K) @ f0)
                res[mode] = dict(EQE=EQE, EQE_L=EQEL, CE=CE, LCF=EQE / EQEL)
            out[f'{case}|{thin}|{th_h}'] = res
            for m, v in res.items():
                print(f"{case:4s} {thin:7s} Theta {th_h}: {m:6s} EQE {v['EQE']*100:5.1f}%  EQE_L {v['EQE_L']*100:5.1f}%  CE {v['CE']:6.1f} cd/A  LCF {v['LCF']:.3f}", flush=True)
json.dump(out, open('data/advsci_onaxis.json', 'w'), indent=1)

"""Measured 30-deg s-pol reflectances of the reference device (2026-10-27, data/refl_261027):
  reflector_airAg: aux reflector seen from air (Ag quality),  center_auxR: glass | organics | TFE | parylene-C | MoOx | Ag (no cathode),
  onpixel: beam centred on the active area (pixel + surroundings).
Envelope comparison (50-nm windows: mean and peak-to-peak) of center_auxR with candidate absorbers; Ag n,k scaled to match the
air-side measurement.  Usage: python3 refl_center_fit.py"""
import sys, numpy as np
sys.path[:0] = ['/home/user/OLED-/sim/diag', '/home/user/OLED-/sim/fig3c']
import fit_ref_angles as F
from scipy.optimize import minimize_scalar
D = F.D; A = np.radians(30.0); NG = 1.52
def load(f):
    a = np.genfromtxt('data/refl_261027/' + f, delimiter=',', skip_header=2, usecols=(0, 1)); a = a[~np.isnan(a[:, 0])][::-1]
    return a[:, 0], a[:, 1] / 100
def tmm_s(ns, ds, lam, n0, kx):
    """s-pol r of layers ns/ds between semi-infinite n0 (incident) and ns[-1] (exit, d ignored)"""
    k0 = 2 * np.pi / lam; q = lambda n: np.sqrt(n ** 2 * k0 ** 2 - kx ** 2 + 0j)
    M = np.eye(2, dtype=complex)
    for n, d in zip(ns[:-1], ds[:-1]):
        qz = q(n); dl = qz * d; Y = qz
        M = M @ np.array([[np.cos(dl), -1j * np.sin(dl) / Y], [-1j * Y * np.sin(dl), np.cos(dl)]])   # n + ik convention
    Y0, Ys = q(n0), q(ns[-1])
    B, C = M @ np.array([1, Ys]); return (Y0 * B - C) / (Y0 * B + C)
F_, G_ = 172.5 / 185, 160.0 / 210
def center_R(lam, dpar, dag=0.0, kag=0.0, kmo=0.0, kh=0.0, kd=0.0, kcs=0.0, nmo=2.1, tfe=30.0):
    j = int(round(lam - 400)); npar = 1.62 + 0.0065 / (lam / 1000) ** 2
    ag = D.AG[j] + dag + 1j * kag; hat = D.HAT[j] + 1j * kh
    ns = [hat, D.TPo[j], hat, D.TPo[j], hat, D.TPo[j], D.TCo[j], D.EMo[j], D.B3o[j], D.B3o[j] + 1j * kd, 1.7 + 1j * kcs, 1.75, npar, nmo + 1j * kmo, ag, 1.0]
    ds = [5, 70 * G_, 5, 70 * G_, 5, 70 * G_, 10, 25, 120 * F_, 60 * F_, 5, tfe, dpar, 5, 100, 0]
    kx = 2 * np.pi / lam * np.sin(A); r = tmm_s(ns, ds, lam, NG, kx); R2 = abs(r) ** 2
    ct = np.cos(np.arcsin(np.sin(A) / NG)); R1 = ((np.cos(A) - NG * ct) / (np.cos(A) + NG * ct)) ** 2
    return R1 + (1 - R1) ** 2 * R2 / (1 - R1 * R2)
def air_ag(lam, dag, kag):
    j = int(round(lam - 400)); ag = D.AG[j] + dag + 1j * kag; r = tmm_s([ag, 1.0], [100, 0], lam, 1.0, 2 * np.pi / lam * np.sin(A)); return abs(r) ** 2
WIN = [(420, 470), (450, 500), (500, 550), (550, 600), (600, 650), (650, 700)]
def env(L, R):
    return np.array([[R[(L >= a) & (L < b)].mean(), np.ptp(R[(L >= a) & (L < b)])] for a, b in WIN])
if __name__ == '__main__':
    La, Ra = load('reflector_airAg_30s.csv'); m = (La >= 420) & (La <= 700)
    r = minimize_scalar(lambda dn: np.sum((np.array([air_ag(l, dn, 0) for l in La[m]]) - Ra[m]) ** 2), bounds=(0, 0.3), method='bounded')
    dag = r.x; print('Ag: n + %.3f fits the air-side R (rms %.2f %%)' % (dag, 100 * np.sqrt(r.fun / m.sum())))
    Lc, Rc = load('center_auxR_30s.csv'); E = env(Lc, Rc)
    print('measured center  mean/ptp:', ' '.join('%.1f/%.1f' % tuple(100 * e) for e in E))
    lf = np.arange(420, 701, 0.5)
    def model_env(**kw):
        best = None
        for dpar in np.arange(2900, 3600, 25):
            R = np.array([center_R(l, dpar, dag=dag, **kw) for l in lf]); e = env(lf, R)
            c = np.sum((e[:, 0] - E[:, 0]) ** 2) + 0.5 * np.sum((e[:, 1] - E[:, 1]) ** 2)
            if best is None or c < best[0]: best = (c, dpar, e)
        return best
    for name, key, vals in (('none', None, [0]), ('MoOx k', 'kmo', [0.1, 0.3, 0.6, 1.0]), ('HATCN k', 'kh', [0.1, 0.2, 0.4]),
                            ('doped ETL k', 'kd', [0.02, 0.04, 0.08]), ('Cs layer k', 'kcs', [0.6, 1.5, 3.0])):
        for v in vals:
            c, dpar, e = model_env(**({key: v} if key else {}))
            print('%-12s %-5s d_par %4d  rms %.1f  ' % (name, v, dpar, 100 * np.sqrt(c / 6)) + ' '.join('%.1f/%.1f' % tuple(100 * x) for x in e), flush=True)

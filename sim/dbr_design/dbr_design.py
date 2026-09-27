"""ZnS/LiF DBR on a glass substrate (n = 1.5), air outside: quarter-wave at 550 nm vs a linearly chirped stack
optimised for the Lambertian (cos.sin) and Ir(ppy)2acac photon-weighted average reflectance, 400-700 nm."""
import sys, numpy as np, csv
sys.path.insert(0, '/home/user/OLED-/sim/design_rule4'); import fig2d as F
from multiprocessing import Pool
NP = 10; NSUB = 1.5
LAM = np.arange(400.0, 701.0, 2.0); sel = np.isin(F.LAM, LAM)
TH = np.radians(np.arange(90) + 0.5)
g = np.cos(TH)*np.sin(TH); g /= g.sum()
S = np.clip(np.interp(LAM, F.LAM, F.GREEN), 0, None) * LAM; S /= S.sum()   # photon-number weight
ZNS, LIF = F.ZNS[sel], F.LIF[sel]
def layers(dz, dl, fac):
    s = np.linspace(1, fac, NP); out = []
    for k in range(NP): out += [(ZNS, dz*s[k]), (LIF, dl*s[k])]
    return out
def RT(lay, lam=LAM):
    return F.RT(lay, np.full(lam.shape, NSUB, complex), np.ones(lam.shape, complex), TH, lam)
def score(a):
    R, T = RT(layers(*a)); return (float(S @ (R @ g)), ) + tuple(a)
if __name__ == '__main__':
    i = int(np.argmin(abs(F.LAM - 550))); qz, ql = 550/(4*F.ZNS[i].real), 550/(4*F.LIF[i].real)
    grid = [(dz, dl, fac) for dz in range(40, 91, 3) for dl in range(60, 141, 5) for fac in np.arange(1.0, 1.81, 0.05)]
    with Pool(4) as p: res = p.map(score, grid, chunksize=40)
    best = max(res)
    b = best[1:]; loc = [(b[0]+x, b[1]+y, b[2]+z) for x in np.arange(-2, 2.1, 0.5) for y in np.arange(-3, 3.1, 1) for z in np.arange(-0.03, 0.031, 0.01)]
    with Pool(4) as p: best = max(p.map(score, loc, chunksize=20))
    dz, dl, fac = best[1:]
    out = {}
    for name, lay in (('quarter-wave 550 nm', layers(qz, ql, 1.0)), ('chirped (optimised)', layers(dz, dl, fac))):
        R, T = RT(lay); A = 1 - R - T
        cone = TH < np.arcsin(1/NSUB); gc = g*cone; gc /= gc.sum()
        out[name] = dict(lay=lay, R=R, T=T, Rav=S @ (R @ g), Tav=S @ (T @ g), Aav=S @ (A @ g), Rcone=S @ (R @ gc), R0=S @ R[:, 0])
        print(f"{name:22s} <R> Lambertian {out[name]['Rav']:.4f}  (T {out[name]['Tav']:.4f}, A {out[name]['Aav']:.4f})  <R> inside escape cone {out[name]['Rcone']:.4f}  R(0 deg, spectrum) {out[name]['R0']:.4f}")
    print('quarter-wave: ZnS %.1f nm, LiF %.1f nm' % (qz, ql))
    print('chirped: ZnS %.1f -> %.1f nm, LiF %.1f -> %.1f nm (factor %.2f)' % (dz, dz*fac, dl, dl*fac, fac))
    np.savez('/home/user/OLED-/sim/dbr_design/dbr_maps.npz', lam=LAM, th=np.degrees(TH),
             R_qw=out['quarter-wave 550 nm']['R'], T_qw=out['quarter-wave 550 nm']['T'], R_ch=out['chirped (optimised)']['R'], T_ch=out['chirped (optimised)']['T'],
             d_qw=[d for _, d in out['quarter-wave 550 nm']['lay']], d_ch=[d for _, d in out['chirped (optimised)']['lay']])

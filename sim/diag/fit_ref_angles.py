"""Fit the reference (flat glass) device to its angle-resolved EL (data/2026-10-22_glass_d32_p5_gon-spec-deg.csv).
Observable: R(lam, theta) = S(lam, theta) / S(lam, 0), theta = 10..70 deg (+-theta averaged; -80 deg kept, +80 deg is an outlier),
5 nm bins over 490-620 nm.  The ratio is free of the intrinsic emitter spectrum, the spectrometer calibration and the
Purcell factor, so it fixes only the cavity: ETL scale f (B3PyMPM 120 + doped 60), HTL scale g (TAPC 3 x 70),
emission-zone position z (nm from the ETL side of the 25 nm EML), Theta and the IZO thickness.
Result (fit_ref_angles.csv, profile over the ETL scale with the dipole at the EML centre): ETL 161 +- 4 nm (x0.87),
Theta 0.71, rms 2 %; the fit needs a weaker bottom (IZO-side) reflection than 50 nm of the library IZO gives
(IZO -> ~15 nm, or n_IZO x0.7 at 50 nm; TAPC then x0.80..0.99).  Robust quantity: dipole-Ag distance 176..181 nm
(nominal 202.5).  With IZO 50 nm fixed the best rms is only 6-7 %.
Air intensity per sr from a flat n_s substrate (incoherent, multiple glass/air <-> OLED bounces, per polarisation):
  I(th_a) ~ P_sub(th_s) T(th_s) / (1 - R_f R_oled) * cos(th_a) / (sin(th_s) n_s^2 cos(th_s)),  sin th_a = n_s sin th_s."""
import numpy as np, io, sys, itertools
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/diag']
import dev3, cps2
D = dev3.D; NS = 1.52
TA = np.array([0., 10, 20, 30, 40, 50, 60, 70, 80]); TS = np.degrees(np.arcsin(np.sin(np.radians(TA)) / NS)); TS[0] = 0.3
LB = np.arange(490., 621., 5.)


def measured(fn='data/2026-10-22_glass_d32_p5_gon-spec-deg.csv'):
    L = open(fn).read().splitlines(); i = [k for k, l in enumerate(L) if l.startswith('wavelength')][0]
    hdr = L[i].split('\t'); dat = np.loadtxt(io.StringIO('\n'.join(L[i + 1:])), delimiter='\t'); lam = dat[:, 0]
    sig = lambda n: dat[:, hdr.index(n)] - dat[:, hdr.index(n + '_bg')]
    binned = lambda s: np.array([s[(lam >= l - 2.5) & (lam < l + 2.5)].mean() for l in LB])
    S = []
    for a in TA:
        if a == 0: S.append(binned(sig('0.0')))
        elif a == 80: S.append(binned(sig('-80.0')))
        else: S.append(0.5 * (binned(sig('%.1f' % a)) + binned(sig('%.1f' % -a))))
    S = np.array(S)                                   # (angle, lam)
    asym = np.array([binned(sig('%.1f' % a)).sum() / binned(sig('%.1f' % -a)).sum() for a in TA[1:]])
    return S, asym


def stack(j, f, g, z, dizo=50.0, e=1.0, kd=0.0, kh=0.0, kcs=0.0):
    """kd: k of the 10 wt% Cs2CO3-doped ETL, kh: k of the three 5 nm HATCN layers (HATCN/TAPC interface CT absorption),
    kcs: k of the 5 nm 50 wt% Cs2CO3 layer at the Ag interface"""
    B3o, B3e, EMo, EMe, TCo, TCe, TPo, TPe, HAT, AG = (D.B3o[j], D.B3e[j], D.EMo[j], D.EMe[j], D.TCo[j], D.TCe[j], D.TPo[j], D.TPe[j], D.HAT[j], D.AG[j])
    return cps2.Stack(float(D.LAM[j]), (EMo, EMe), 25.0 * e, z,
                      above=[(B3o, B3e, 120.0 * f), (B3o + 1j * kd, B3e + 1j * kd, 60.0 * f), (1.7 + 1j * kcs, 1.7 + 1j * kcs, 5.0), (AG, AG, 100.0)],
                      below=[(TCo, TCe, 10.0 * e), (TPo, TPe, 70.0 * g), (HAT + 1j * kh, HAT + 1j * kh, 5.0 * e), (TPo, TPe, 70.0 * g), (HAT + 1j * kh, HAT + 1j * kh, 5.0 * e),
                             (TPo, TPe, 70.0 * g), (HAT + 1j * kh, HAT + 1j * kh, 5.0 * e), (dev3.IZO[j], dev3.IZO[j], dizo)], n_sub=NS, h=1.0)


def air(S):
    """(vertical, horizontal) air intensity at TA, per sr (arbitrary common unit)"""
    k0 = 2 * np.pi / S.lam; th = np.radians(TS); kx = NS * np.sin(th) * k0
    ct = np.sqrt(1 - (NS * np.sin(th)) ** 2); cs = np.cos(th)
    rf = {'s': (NS * cs - ct) / (NS * cs + ct), 'p': (NS * ct - cs) / (NS * ct + cs)}
    geo = ct / (np.sin(th) * NS ** 2 * cs)
    out = []
    for h in (0.0, 1.0):
        S.h = h; tot = 0
        for pol in ('p', 's'):
            P = np.real(cps2.sub_angular(S, TS, pols=(pol,)))
            ro = np.abs(cps2.stack_rt(S.all_no, S.all_ne, S.all_d, kx, k0, pol)[0]) ** 2
            Rf = np.abs(rf[pol]) ** 2
            tot = tot + P * (1 - Rf) / (1 - Rf * ro)
        out.append(tot * geo)
    return np.array(out)                              # (2, angle)


JB = [int(round(l - 400)) for l in LB]
def model(a):
    return a, np.array([air(stack(j, *a)) for j in JB])      # (lam, 2, angle)


def cost(S, f, g, z, t, dizo):
    """rms of log(R_model / R_meas) over 10..70 deg and the lam bins with > 15 % of the normal peak"""
    Rm = S / S[0]; w = (S[0] / S[0].max() > 0.15).astype(float)
    M = model((f, g, z, dizo))[1]; I = (1 - t) * M[:, 0] + t * M[:, 1]; R = (I / I[:, :1]).T
    e = np.log(R[1:8] / Rm[1:8]) * w
    return float(np.sqrt((e ** 2).sum() / (w.sum() * 7)))


Z = 12.5                                  # dipole distance from the ETL/EML interface [nm]; `python3 fit_ref_angles.py <Z>`


def profile(f):
    """ETL scale fixed, TAPC scale / Theta / IZO thickness free, dipole at Z (default: EML centre)"""
    from scipy.optimize import minimize
    S = measured()[0]
    r = minimize(lambda x: cost(S, f, x[0], Z, x[1], x[2]), [0.85, 0.71, 20.0], method='Nelder-Mead',
                 bounds=[(0.4, 1.5), (0.5, 1.0), (0.0, 90.0)], options=dict(maxiter=300, xatol=1e-3, fatol=1e-6))
    return f, r.fun, *r.x


if __name__ == '__main__':
    Z = float(sys.argv[1]) if len(sys.argv) > 1 else 12.5
    S, asym = measured()
    print('+theta / -theta (lam-integrated), 10..80 deg:', np.round(asym, 3))
    with Pool(4) as p:
        rows = p.map(profile, np.round(np.arange(0.80, 0.961, 0.02) if Z == 12.5 else np.arange(0.86, 1.021, 0.02), 3))
    import csv
    with open('fit_ref_angles.csv' if Z == 12.5 else 'fit_ref_angles_z%g.csv' % Z, 'w', newline='') as fo:
        wr = csv.writer(fo); wr.writerow(['f_ETL', 'ETL_nm', 'rms_log', 'g_TAPC', 'Theta', 'IZO_nm'])
        for r in rows:
            wr.writerow([r[0], round(185 * r[0], 1)] + [round(v, 4) for v in r[1:]])
            print('f %.2f ETL %.0f nm  rms %.4f  TAPC x%.3f  Theta %.3f  IZO %.1f nm' % (r[0], 185 * r[0], *r[1:]), flush=True)

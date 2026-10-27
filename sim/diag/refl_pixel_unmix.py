"""Pixel-only reflectance from the on-pixel measurement (data/refl_261027/onpixel_30s.csv, 30 deg s-pol).
The beam also covers the surroundings, which carry the parylene-cavity ripple; the pixel (opaque 100 nm Ag cathode) cannot.
Model:  R_meas = (1 - s1 - s2) R_pix + s1 R_strip + s2 R_center
  R_center: measured (center_auxR_30s.csv, glass | organics | TFE | parylene | MoOx | Ag)
  R_strip : modelled anode strip next to the pixel (glass | IZO | organics | TFE | parylene d | MoOx | Ag), absorber as fitted
            to R_center (refl_center_fit: HATCN k 0.2 or MoOx k 1.0)
  R_pix   : smooth, R_pix_model(no extra absorption) - (c0 + c1 x + c2 x^2), x = (lam - 575) / 125
Fit 450-700 nm.  -> refl_pixel_unmix.csv"""
import numpy as np, csv
from scipy.optimize import least_squares
import refl_center_fit as C
D = C.D; A = C.A; NG = C.NG
F_, G_ = C.F_, C.G_
DAG = 0.099
def pix_R(lam, dizo=50.0):
    j = int(round(lam - 400)); ag = D.AG[j] + DAG; izo = C.F.dev3.IZO[j]
    ns = [izo, D.HAT[j], D.TPo[j], D.HAT[j], D.TPo[j], D.HAT[j], D.TPo[j], D.TCo[j], D.EMo[j], D.B3o[j], D.B3o[j], 1.7, ag, 1.0]
    ds = [dizo, 5, 70 * G_, 5, 70 * G_, 5, 70 * G_, 10, 25, 120 * F_, 60 * F_, 5, 100, 0]
    r = C.tmm_s(ns, ds, lam, NG, 2 * np.pi / lam * np.sin(A)); R2 = abs(r) ** 2
    ct = np.cos(np.arcsin(np.sin(A) / NG)); R1 = ((np.cos(A) - NG * ct) / (np.cos(A) + NG * ct)) ** 2
    return R1 + (1 - R1) ** 2 * R2 / (1 - R1 * R2)
def strip_R(lam, dpar, kh=0.0, kmo=0.0, dizo=50.0):
    j = int(round(lam - 400)); npar = 1.62 + 0.0065 / (lam / 1000) ** 2; ag = D.AG[j] + DAG; hat = D.HAT[j] + 1j * kh
    izo = C.F.dev3.IZO[j]
    ns = [izo, hat, D.TPo[j], hat, D.TPo[j], hat, D.TPo[j], D.TCo[j], D.EMo[j], D.B3o[j], D.B3o[j], 1.7, 1.75, npar, 2.1 + 1j * kmo, ag, 1.0]
    ds = [dizo, 5, 70 * G_, 5, 70 * G_, 5, 70 * G_, 10, 25, 120 * F_, 60 * F_, 5, 30, dpar, 5, 100, 0]
    r = C.tmm_s(ns, ds, lam, NG, 2 * np.pi / lam * np.sin(A)); R2 = abs(r) ** 2
    ct = np.cos(np.arcsin(np.sin(A) / NG)); R1 = ((np.cos(A) - NG * ct) / (np.cos(A) + NG * ct)) ** 2
    return R1 + (1 - R1) ** 2 * R2 / (1 - R1 * R2)
if __name__ == '__main__':
    Lp, Rp = C.load('onpixel_30s.csv'); Lc, Rc = C.load('center_auxR_30s.csv')
    m = (Lp >= 450) & (Lp <= 700); L = Lp[m]; Y = Rp[m]; RC = np.interp(L, Lc, Rc); x = (L - 575) / 125
    P0 = np.array([pix_R(l) for l in L])
    best = None
    for absb in (dict(kh=0.2), dict(kmo=1.0)):
        for dpar in np.arange(2900, 3601, 20):
            RS = np.array([strip_R(l, dpar, **absb) for l in L])
            def res(p):
                s1, s2, c0, c1, c2 = p
                return (1 - s1 - s2) * (P0 - (c0 + c1 * x + c2 * x ** 2)) + s1 * RS + s2 * RC - Y
            r = least_squares(res, [0.2, 0.2, 0.03, 0, 0], bounds=([0, 0, -0.05, -0.2, -0.2], [0.8, 0.8, 0.3, 0.2, 0.2]))
            c = np.sqrt(np.mean(r.fun ** 2))
            if best is None or c < best[0]: best = (c, dpar, absb, r.x, RS)
    c, dpar, absb, p, RS = best; s1, s2, c0, c1, c2 = p
    Rpix = P0 - (c0 + c1 * x + c2 * x ** 2)
    print('best: absorber %s  d_par %d nm  s_strip %.2f  s_center %.2f  rms %.2f %%' % (absb, dpar, s1, s2, 100 * c))
    for l in (450, 500, 550, 600, 650, 700):
        i = np.argmin(abs(L - l)); print('  %d nm: R_pix %.1f %%  (model without extra absorption %.1f %%, measured on-pixel %.1f %%)' % (l, 100 * Rpix[i], 100 * P0[i], 100 * Y[i]))
    I = np.clip(np.asarray(D.spec.l_I_Irppy2acac, float)[:301], 0, None); W = np.interp(L, D.LAM, I * D.LAM); W /= W.sum()
    print('  Ir(ppy)2acac-weighted: R_pix %.1f %%, model %.1f %%' % (100 * W @ Rpix, 100 * W @ P0))
    with open('refl_pixel_unmix.csv', 'w', newline='') as fo:
        w = csv.writer(fo); w.writerow(['lam', 'R_meas_onpixel', 'R_pix_fit', 'R_pix_model_noabs', 'R_strip_model', 'R_center_meas', 'R_fit_total'])
        Rt = (1 - s1 - s2) * Rpix + s1 * RS + s2 * RC
        for row in zip(L, Y, Rpix, P0, RS, RC, Rt): w.writerow(['%.0f' % row[0]] + ['%.4f' % v for v in row[1:]])

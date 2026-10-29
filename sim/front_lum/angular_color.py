"""Angular luminance and colour (CIE 1931 -> u'v') of a recycling-OLED structure and of a reference, 0-70 deg in air.
CIE CMFs: multi-lobe Gaussian fit (Wyman, Sloan & Shirley 2013).  Usage: python3 angular_color.py <json> <case> <N>"""
import sys, json, numpy as np
import ce_fab as F, ce_opt as C, front_lum as FL, cps2
D = C.D
def g(l, m, s1, s2): s = np.where(l < m, s1, s2); return np.exp(-0.5 * ((l - m) / s) ** 2)
def cmf(l):
    x = 1.056 * g(l, 599.8, 37.9, 31.0) + 0.362 * g(l, 442.0, 16.0, 26.7) - 0.065 * g(l, 501.1, 20.4, 26.2)
    y = 0.821 * g(l, 568.8, 46.9, 40.5) + 0.286 * g(l, 530.9, 16.3, 31.1)
    z = 1.217 * g(l, 437.0, 11.8, 36.0) + 0.681 * g(l, 459.0, 26.0, 13.8); return x, y, z
TA = np.array([0, 10, 20, 30, 40, 50, 60, 70.0])
def radiance(case, x, N, j, S=10.0, gg=0.0, mua=0.01):
    """spectral radiance factor in air at TA (per unit emitted power) and the Purcell F"""
    if case == 'TREF':
        P, *_ = F.point((case, j, tuple(x), N, 2000))
        dF, dB, dAg, dC = x; lam = 400.0 + j
        S_ = cps2.Stack(lam, (D.EMo[j], D.EMe[j]), 25.0, 12.5,
                        above=[(D.TPo[j], D.TPe[j], dB), (D.HAT[j], D.HAT[j], 5.0), (C.IZO[j], C.IZO[j], 10.0), (D.AG[j], D.AG[j], 100.0)],
                        below=[(D.B3o[j], D.B3e[j], dF), (D.HAT[j], D.HAT[j], 5.0), (F.agthin(dAg, j),) * 2 + (dAg,), (1.9, 1.9, dC)], n_sub=1.0, h=C.THETA)
        th = np.maximum(TA, 0.25); I = np.real(cps2.sub_angular(S_, th)) * 180 / np.pi / (2 * np.pi * np.sin(np.radians(th))) / P
        return I / np.cos(np.radians(TA)), P
    P, Ld, Pair, (Ps, Tf, Rf, ns) = F.point((case, j, tuple(x), N, 2000))
    fr, back, lamc, _ = F.layers(case, j, x); fr = fr + C.dbr(N, lamc, j); lam = 400.0 + j
    S1 = cps2.Stack(lam, (D.EMo[j], D.EMe[j]), 25.0, 12.5, above=back, below=fr, n_sub=1.0, n_top=ns, h=C.THETA)
    th = np.maximum(TA, 0.25); Id = np.real(cps2.sub_angular(S1, th)) * 180 / np.pi / (2 * np.pi * np.sin(np.radians(th))) / P
    B = C.BZ['S%g_g%g_a%g' % (S, gg, mua)]; U = np.linalg.solve(np.eye(C.NB) - B * Rf[None, :], B @ Ps)
    tm = np.radians(C.TB); Om = 2 * np.pi * (np.cos(tm - np.radians(1)) - np.cos(tm + np.radians(1)))
    Lm = U / (Om * np.cos(tm)); ok = ns * np.sin(tm) < 1; ta = np.degrees(np.arcsin(ns * np.sin(tm[ok])))
    La = np.interp(TA, np.r_[0, ta], np.r_[(Lm * Tf)[0], (Lm * Tf)[ok]] / ns ** 2, right=0)
    return Id / np.cos(np.radians(TA)) + La, P
def run(case, x, N, label):
    js = list(range(0, 301, 5)); lam = np.array(js) + 400.0; p = C.I[js] * lam; p /= p.sum()
    R = []; Fs = []
    for j in js:
        r, P = radiance(case, x, N, j); R.append(r); Fs.append(P)
    R = np.array(R); Fs = np.array(Fs); k = C.Q * p * Fs / (1 - C.Q + C.Q * (p @ Fs))       # photons per electron per lambda
    E = (k * 1239.84 / lam)[:, None] * R                                                   # W/sr/m2 per A/m2, per lambda
    X, Y, Z = [683 * (c[:, None] * E).sum(0) for c in cmf(lam)]
    u = 4 * X / (X + 15 * Y + 3 * Z); v = 9 * Y / (X + 15 * Y + 3 * Z)
    print(label); print('  angle  ' + ' '.join('%6.0f' % a for a in TA))
    print('  L/J    ' + ' '.join('%6.1f' % y for y in Y) + '  (cd/A)')
    print('  rel    ' + ' '.join('%6.2f' % (y / Y[0]) for y in Y))
    print('  du\'v\'  ' + ' '.join('%6.3f' % np.hypot(u[i] - u[0], v[i] - v[0]) for i in range(len(TA))))
    return Y, u, v
if __name__ == '__main__':
    fn, case, N = sys.argv[1], sys.argv[2], int(sys.argv[3])
    x = list(json.load(open(fn))['structure'].values()); run(case, x, N, fn)

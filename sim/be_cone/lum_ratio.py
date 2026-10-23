"""Front (on-axis) signal on the goniometer detector: reflector over the whole back (B) vs Ag only in the 2x2 mm pixel (A).
Finite substrate MC as finite_substrate.py, but the escaping weight is kept per internal incidence bin and converted to
air directions with the cone array's transmitted-direction matrix (ray tracer, 1-deg incidence bins).
Detector: 75 mm^2 disc at 50 / 80 mm on the normal (half-angle 5.6 / 3.5 deg); extended source neglected (it changes the
integrated signal by < 0.1 %, gonio_check.py)."""
import numpy as np, sys
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/be_cone', '/home/user/OLED-/sim/mla']
import raytrace as rt
import finite_substrate as F

R_UM, H_UM, GAP = 1.4, 1.8, 0.2
AB = np.arange(0, 90.01, 0.5)                                   # air-angle bins 0.5 deg


def tcol(i):
    rng = np.random.default_rng([11, i])
    esc, _, Dd = rt._trace_chunk(H_UM / R_UM, 1.77, 40000, 80, rng, float(F.TH[i]), 'cone', (2 * R_UM + GAP) / R_UM, return_dir=True)
    th = np.degrees(np.arccos(np.clip(Dd[esc, 2], -1, 1)))
    return np.histogram(th, AB)[0] / 40000.0                     # fraction of rays incident in bin i leaving into each air bin


def escbins(args):
    seed, N, scen = args
    rng = np.random.default_rng(seed)
    j = rng.choice(len(WL), N, p=WL); Pc = np.cumsum(PS / PS.sum(1, keepdims=True), 1)
    ib = (Pc[j] < rng.random(N)[:, None]).sum(1).clip(0, 89)
    x = (rng.random(N) - 0.5) * 2 * F.PIX; y = (rng.random(N) - 0.5) * 2 * F.PIX; ph = rng.random(N) * 2 * np.pi; w = np.ones(N)
    E = np.zeros(90); alive = np.ones(N, bool)
    for _ in range(600):
        a = np.flatnonzero(alive)
        if a.size == 0: break
        for half in (0, 1):
            th = np.radians(F.TH[ib[a]] + rng.random(a.size) - 0.5); dd = F.T * np.tan(th)
            x[a] += dd * np.cos(ph[a]); y[a] += dd * np.sin(ph[a])
            out = (np.abs(x[a]) > F.HALF) | (np.abs(y[a]) > F.HALF); w[a[out]] = 0
            if half == 0:
                bt = F.BT[ib[a]]; E += np.bincount(ib[a], weights=w[a] * bt, minlength=90)
                w[a] *= 1 - bt; u = rng.random(a.size)
                ib[a] = (F.BRc[:, ib[a]].T < u[:, None]).sum(1).clip(0, 89); ph[a] = rng.random(a.size) * 2 * np.pi
        inpix = (np.abs(x[a]) <= F.PIX) & (np.abs(y[a]) <= F.PIX); Rd = RS[j[a], ib[a]]
        w[a] *= np.where(inpix, Rd, F.RF[ib[a]]) if scen == 'A' else Rd
        alive &= w > 1e-5
    return E / N


if __name__ == '__main__':
    with Pool(4) as p:
        TA = np.array(p.map(tcol, range(90))).T                     # (air bins, incidence bins)
        res = p.map(F.stack, F.LAMJ)
    ES = np.array([r[0] for r in res]); PS = np.array([r[1] for r in res]); RS = np.array([r[2] for r in res])
    I = np.clip(np.asarray(F.D.spec.l_I_Irppy2acac, float)[:301], 0, None); wp = (F.D.LAM * I)[F.LAMJ]; wp /= wp.sum()
    es = float(wp @ ES); WL = wp * ES / (wp @ ES)
    print('check: ray-traced total transmission vs stored B_T (bins 10/40/60): %s vs %s'
          % (np.round(TA.sum(0)[[10, 40, 60]], 3), np.round(F.BT[[10, 40, 60]], 3)))
    out = {}
    for scen in ('A', 'B'):
        with Pool(4) as p: E = np.mean(p.map(escbins, [(100 * k + 9, 200000, scen) for k in range(8)]), 0)
        A = TA @ (E / np.maximum(TA.sum(0), 1e-12) * TA.sum(0) / np.maximum(F.BT, 1e-12))   # escaping weight -> air bins
        sa = 2 * np.pi * (np.cos(np.radians(AB[:-1])) - np.cos(np.radians(AB[1:])))
        out[scen] = (E.sum(), A, A / sa)
    TCa = 0.5 * (AB[1:] + AB[:-1])
    for L in (50.0, 80.0):
        a = np.degrees(np.arctan(np.sqrt(75 / np.pi) / L)); sel = TCa < a
        pdA, pdB = out['A'][1][sel].sum(), out['B'][1][sel].sum()
        print('detector at %2.0f mm (half-angle %.1f deg): EQE-equivalent signal A %.5f  B %.5f  ratio B/A %.3f'
              % (L, a, es * pdA, es * pdB, pdB / pdA))
    print('front-hemisphere totals: A %.4f  B %.4f  ratio B/A %.3f' % (es * out['A'][0], es * out['B'][0], out['B'][0] / out['A'][0]))
    for th in (0, 20, 40, 60, 75):
        k = np.argmin(abs(TCa - (th + 0.25))); kk = slice(max(k - 4, 0), k + 5)
        print('  I_B/I_A near %2d deg: %.3f' % (th, out['B'][2][kk].mean() / out['A'][2][kk].mean()))

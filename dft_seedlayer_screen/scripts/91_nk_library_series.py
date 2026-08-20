"""Re-run the seed series against the lab's own measured optical constants.

nk_JH_total.mat supplies McPeak and Palik silver plus this lab's measured
HATCN, MoO3, ITO, IZO and organics; they are exported to data/nk/*.csv.

Two corrections to the earlier analysis:
  - the seed indices were assumed (HATCN 1.75 + 0.02i, MoOx 2.10 + 0.01i);
    the measured values are HATCN 1.672 + 0.049i and MoO3 2.112 + 0.003i
  - the ideal-silver floor was a Drude model tuned to n(550) = 0.04; it is
    now McPeak's measured data, n = 0.0438, k = 3.6101 at 550 nm
"""
import csv, json, os
import numpy as np
from scipy.optimize import least_squares

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
NK, RAW = os.path.join(BASE, "data", "nk"), os.path.join(BASE, "data", "TR_20260820", "raw")
KEEP, D_SEED, L = 0.843, 5.0, 550.0
SAMPLES = [("1-2","HATCN",4), ("1-3","HATCN",5), ("1-4","HATCN",6), ("2-1","HATCN",7),
           ("2-2","HATCN",8), ("2-3","HATCN",10), ("2-4","HATCN",12),
           ("1-10","MoOx",4), ("1-11","MoOx",5), ("1-12","MoOx",6),
           ("2-9","MoOx",7), ("2-10","MoOx",8), ("2-11","MoOx",10)]
GLASS_TI = np.array([[350,.9226],[400,.9952],[450,1.],[550,1.],[650,.9987],
                     [700,.9973],[800,.9911],[850,.9920]])


def load_nk(name):
    a = np.loadtxt(os.path.join(NK, f"{name}.csv"), delimiter=",", skiprows=1)
    return lambda l: complex(np.interp(l, a[:,0], a[:,1]), np.interp(l, a[:,0], a[:,2]))


AG_IDEAL = load_nk("Ag_McPeak")
SEED_NK = {"HATCN": load_nk("HATCN"), "MoOx": load_nk("MoO3")}
def n_glass(l): return 1.5220 + 3900./l**2
def t_glass(l): return np.interp(l, GLASS_TI[:,0], GLASS_TI[:,1])


def read_csv(p):
    out = {}
    for row in csv.reader(open(p)):
        if len(row) < 2: continue
        try: out[float(row[0])] = float(row[1])/100.0
        except ValueError: pass
    return out


def tr(n, d, lam):
    n = [complex(x) for x in n]; k0 = 2*np.pi/lam
    M = np.eye(2, dtype=complex)
    for j in range(len(n)-1):
        r = (n[j]-n[j+1])/(n[j]+n[j+1]); t = 2*n[j]/(n[j]+n[j+1])
        I = np.array([[1,r],[r,1]], dtype=complex)/t
        if j+1 < len(n)-1:
            dl = k0*n[j+1]*d[j+1]
            I = I @ np.array([[np.exp(-1j*dl),0],[0,np.exp(1j*dl)]])
        M = M @ I
    return (n[-1].real/n[0].real)*abs(1/M[0,0])**2, abs(M[1,0]/M[0,0])**2


def observable(nk, lam, seed, d_ag):
    ng, ti = n_glass(lam), t_glass(lam); ns = SEED_NK[seed](lam)
    Rb = ((ng-1)/(ng+1))**2
    T1, R1 = tr([1.0, nk, ns, ng], [0, d_ag, D_SEED, 0], lam)
    T1b, R1b = tr([ng, ns, nk, 1.0], [0, D_SEED, d_ag, 0], lam)
    den = 1 - Rb*R1b*ti**2
    return T1*ti*(1-Rb)/den, R1 + KEEP*T1*ti**2*Rb*T1b/den


def invert(lam, T, R, seed, d_ag, guess=(0.4, 3.5)):
    f = lambda x: [observable(complex(abs(x[0]), abs(x[1])), lam, seed, d_ag)[0]-T,
                   observable(complex(abs(x[0]), abs(x[1])), lam, seed, d_ag)[1]-R]
    s = least_squares(f, list(guess), xtol=1e-14, ftol=1e-14)
    return complex(abs(s.x[0]), abs(s.x[1])), float(np.max(np.abs(s.fun)))


def best_device_A(nk, lam, d_ag, n_cpl=2.10):
    best = (9e9, 0.0)
    for dc in np.arange(20, 181, 1.0):
        T, R = tr([1.8, nk, n_cpl, 1.0], [0, d_ag, dc, 0], lam)
        if 1-T-R < best[0]: best = (1-T-R, dc)
    return best


def main():
    ag = AG_IDEAL(L)
    print(f"McPeak silver at {L:.0f} nm:  n = {ag.real:.4f}  k = {ag.imag:.4f}"
          f"   eps1 = {(ag**2).real:.2f}  eps2 = {(ag**2).imag:.4f}")
    print(f"measured seeds:  HATCN {SEED_NK['HATCN'](L):.4f}"
          f"   MoO3 {SEED_NK['MoOx'](L):.4f}\n")

    meas = {}
    for sid, seed, d in SAMPLES:
        pt, pr = os.path.join(RAW, f"{sid}T.csv"), os.path.join(RAW, f"{sid}R.csv")
        if not (os.path.exists(pt) and os.path.exists(pr)): continue
        T, R = read_csv(pt).get(L), read_csv(pr).get(L)
        if T is None or R is None: continue
        nk, res = invert(L, T, R, seed, d)
        a, dc = best_device_A(nk, L, float(d))
        meas.setdefault(seed, {})[d] = dict(T=T, R=R, A=1-T-R, nk=nk,
                                            eps=nk**2, A_dev=a, cpl=dc, resid=res)

    e_ag = (ag**2).imag
    print("DEVICE ONE-PASS ABSORPTION at 550 nm, organic 1.8 / Ag / CPL n=2.1 / air")
    print("CPL thickness re-optimised for every entry\n")
    print(f"{'d_Ag':>5}{'McPeak ideal':>14} | {'HATCN':>9}{'x ideal':>9}{'eps2/blk':>10}"
          f" | {'MoOx':>9}{'x ideal':>9}{'eps2/blk':>10}")
    print("-"*78)
    for d in [3,4,5,6,7,8,10,12,15,20]:
        ai, _ = best_device_A(ag, L, float(d))
        row = f"{d:5d}{100*ai:13.3f}% |"
        for seed in ("HATCN", "MoOx"):
            v = meas.get(seed, {}).get(d)
            if v:
                row += (f"{100*v['A_dev']:8.2f}%{v['A_dev']/ai:8.1f}x"
                        f"{v['eps'].imag/e_ag:10.1f}") + (" |" if seed=="HATCN" else "")
            else:
                row += f"{'-':>9}{'-':>9}{'-':>10}" + (" |" if seed=="HATCN" else "")
        print(row)

    print("\nInversion with the measured seed indices, for comparison with the")
    print("earlier run that assumed HATCN 1.75 + 0.02i:\n")
    print(f"{'seed':<7}{'d_Ag':>5} | {'n':>7}{'k':>7}{'eps1':>8}{'eps2':>7}{'CPL':>7}{'A_dev':>8}")
    print("-"*58)
    for seed in ("HATCN", "MoOx"):
        for d in sorted(meas.get(seed, {})):
            v = meas[seed][d]
            print(f"{seed:<7}{d:5d} | {v['nk'].real:7.3f}{v['nk'].imag:7.3f}"
                  f"{v['eps'].real:8.2f}{v['eps'].imag:7.3f}{v['cpl']:6.0f}nm"
                  f"{100*v['A_dev']:7.2f}%")

    json.dump({s: {str(d): dict(T=v['T'], R=v['R'], A=v['A'], n=v['nk'].real,
                                k=v['nk'].imag, eps1=v['eps'].real, eps2=v['eps'].imag,
                                A_dev=v['A_dev'], cpl=v['cpl'])
                   for d, v in dd.items()} for s, dd in meas.items()},
              open(os.path.join(BASE, "runs", "TR_20260820_measured_nk.json"), "w"),
              indent=2, default=float)


if __name__ == "__main__":
    main()

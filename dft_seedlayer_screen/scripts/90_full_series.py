"""Full seed x thickness series, 2026-08-20.  Cary 6000i + UMA, absolute T and R.

Raw exports live in data/TR_20260820/raw as <id>T.csv / <id>R.csv.

    1-2,1-3,1-4      HATCN 5 nm / Ag 4, 5, 6
    2-1,2-2,2-3,2-4  HATCN 5 nm / Ag 7, 8, 10, 12
    1-9..1-12        MoOx  5 nm / Ag 0, 4, 5, 6
    2-9..2-12        MoOx  5 nm / Ag 7, 8, 10, 12

R is corrected for the 84.3 % back-surface collection measured on bare glass,
and the substrate's own absorption is carried as an internal transmittance.
Each wavelength is then inverted for the silver layer's n and k.
"""
import csv, os, sys, json
import numpy as np
from scipy.optimize import least_squares

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RAW = os.path.join(BASE, "data", "TR_20260820", "raw")
KEEP, D_SEED = 0.843, 5.0
SEED_N = {"HATCN": 1.75, "MoOx": 2.10}
RS = {("HATCN",4):52.2, ("HATCN",5):23.3, ("HATCN",6):18.9, ("HATCN",7):12.7,
      ("HATCN",8):9.1, ("HATCN",10):7.5, ("HATCN",12):5.3,
      ("MoOx",4):138.0, ("MoOx",5):49.1, ("MoOx",6):32.0, ("MoOx",7):14.9,
      ("MoOx",8):10.8, ("MoOx",10):9.6, ("MoOx",12):6.6}
SAMPLES = [("1-2","HATCN",4), ("1-3","HATCN",5), ("1-4","HATCN",6),
           ("2-1","HATCN",7), ("2-2","HATCN",8), ("2-3","HATCN",10), ("2-4","HATCN",12),
           ("1-9","MoOx",0),  ("1-10","MoOx",4), ("1-11","MoOx",5), ("1-12","MoOx",6),
           ("2-9","MoOx",7),  ("2-10","MoOx",8), ("2-11","MoOx",10), ("2-12","MoOx",12)]

JC = np.array([[326.3,.17,1.95],[357.5,.13,1.81],[397.4,.05,2.07],[430.5,.04,2.36],
               [450.9,.04,2.66],[495.9,.05,3.09],[520.9,.05,3.34],[548.6,.06,3.59],
               [582.1,.05,3.93],[616.8,.06,4.15],[659.5,.05,4.48],[704.5,.041,4.84],
               [756.,.042,5.24],[821.1,.043,5.65]])
GLASS_TI = np.array([[350,.9226],[360,.9612],[370,.9698],[380,.9880],[400,.9952],
                     [450,1.],[500,1.],[550,1.],[600,1.],[650,.9987],[700,.9973],
                     [750,.9978],[780,.9958],[800,.9911],[850,.9920]])


def jc(l): return np.interp(l, JC[:,0], JC[:,1]) + 1j*np.interp(l, JC[:,0], JC[:,2])
def n_glass(l): return 1.5220 + 3900./l**2
def t_glass(l): return np.interp(l, GLASS_TI[:,0], GLASS_TI[:,1])


def read_csv(path):
    out = {}
    for row in csv.reader(open(path)):
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
    ng, ti = n_glass(lam), t_glass(lam)
    Rb = ((ng-1)/(ng+1))**2; ns = SEED_N[seed]
    T1, R1 = tr([1.0, nk, ns, ng], [0, d_ag, D_SEED, 0], lam)
    T1b, R1b = tr([ng, ns, nk, 1.0], [0, D_SEED, d_ag, 0], lam)
    den = 1 - Rb*R1b*ti**2
    return T1*ti*(1-Rb)/den, R1 + KEEP*T1*ti**2*Rb*T1b/den


def invert(lam, T, R, seed, d_ag, guess):
    f = lambda x: [observable(complex(abs(x[0]), abs(x[1])), lam, seed, d_ag)[0]-T,
                   observable(complex(abs(x[0]), abs(x[1])), lam, seed, d_ag)[1]-R]
    s = least_squares(f, guess, xtol=1e-14, ftol=1e-14)
    return complex(abs(s.x[0]), abs(s.x[1])), float(np.max(np.abs(s.fun)))


def device_A(nk, lam, d_ag, n_cpl=2.10, d_cpl=None):
    if d_cpl is None:
        best = 9e9
        for dc in np.arange(20, 161, 1.0):
            T, R = tr([1.8, nk, n_cpl, 1.0], [0, d_ag, dc, 0], lam)
            best = min(best, 1-T-R)
        return best
    T, R = tr([1.8, nk, n_cpl, 1.0], [0, d_ag, d_cpl, 0], lam)
    return 1-T-R


def main():
    L = 550.0
    print(f"{'id':<6}{'seed':<7}{'d_Ag':>5}{'Rs':>7} | {'T':>7}{'R':>7}{'A':>7} | "
          f"{'n':>7}{'k':>7}{'eps1':>8}{'eps2':>7} | {'e2/blk':>8}{'rho/blk':>9}"
          f"{'A_dev':>8}")
    print("-"*104)
    out = {}
    for sid, seed, d_ag in SAMPLES:
        pt = os.path.join(RAW, f"{sid}T.csv"); pr = os.path.join(RAW, f"{sid}R.csv")
        if not (os.path.exists(pt) and os.path.exists(pr)):
            miss = [n for n, p in (("T", pt), ("R", pr)) if not os.path.exists(p)]
            print(f"{sid:<6}{seed:<7}{d_ag:5d}{'':>7} | missing {','.join(miss)}")
            continue
        Td, Rd = read_csv(pt), read_csv(pr)
        if L not in Td or L not in Rd: continue
        T, R = Td[L], Rd[L]
        A = 1-T-R
        rs = RS.get((seed, d_ag))
        if d_ag == 0:
            print(f"{sid:<6}{seed:<7}{'bare':>5}{'':>7} | {100*T:6.2f}%{100*R:6.2f}%"
                  f"{100*A:6.2f}% |"); continue
        nk, res = invert(L, T, R, seed, d_ag, [0.4, 3.5])
        e, eb = nk**2, jc(L)**2
        fdc = rs*d_ag*0.1/1.59 if rs else float("nan")
        ad = device_A(nk, L, d_ag)
        out[sid] = dict(seed=seed, d_ag=d_ag, T=T, R=R, A=A, n=nk.real, k=nk.imag,
                        eps1=e.real, eps2=e.imag, e2_bulk=e.imag/eb.imag,
                        rho_bulk=fdc, A_dev=ad, resid=res)
        print(f"{sid:<6}{seed:<7}{d_ag:5d}{rs:7.1f} | {100*T:6.2f}%{100*R:6.2f}%"
              f"{100*A:6.2f}% | {nk.real:7.3f}{nk.imag:7.3f}{e.real:8.2f}{e.imag:7.3f}"
              f" | {e.imag/eb.imag:8.2f}{fdc:9.2f}{100*ad:7.2f}%")
    print(f"\n  bulk Ag at 550 nm: n 0.060  k 3.590  eps1 {(jc(L)**2).real:.2f}"
          f"  eps2 {(jc(L)**2).imag:.3f}")
    json.dump(out, open(os.path.join(BASE, "runs", "TR_20260820_series.json"), "w"),
              indent=2, default=float)


if __name__ == "__main__":
    main()

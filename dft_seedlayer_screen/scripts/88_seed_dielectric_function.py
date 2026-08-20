"""The seed comparison, settled by the sign of eps1.

Absolute T and R on glass/seed 5 nm/Ag 5 nm, Cary 6000i + UMA, 2026-08-20.
Data columns: wavelength_nm, T_percent, R_percent.  R is corrected for the
84.3 % back-surface collection measured on bare glass, and the substrate's own
Fe3+ absorption is carried as an internal transmittance.

Inverting each wavelength for the Ag layer's n and k gives a qualitative, not
quantitative, split:

    HATCN   eps1 negative everywhere, within ~10 % of bulk silver
    MoOx    eps1 POSITIVE everywhere

A continuous metal film cannot have eps1 > 0 in the visible.  The MoOx film is
therefore not a film -- it is islands, and its fitted n,k are the parameters of
a fictitious uniform slab, not material constants.  They must not be fed into a
device TMM.  The HATCN film is bulk silver with the damping raised about
8.6-fold, which is what the size effect predicts and what its sheet resistance
independently says.
"""
import os
import numpy as np
from scipy.optimize import least_squares

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA = os.path.join(BASE, "data", "TR_20260820")
KEEP = 0.843                       # back-surface collection, from bare glass
D_AG = D_SEED = 5.0
SEED_N = {"HATCN": 1.75, "MoOx": 2.10}

JC = np.array([[326.3,.17,1.95],[357.5,.13,1.81],[397.4,.05,2.07],[430.5,.04,2.36],
               [450.9,.04,2.66],[495.9,.05,3.09],[520.9,.05,3.34],[548.6,.06,3.59],
               [582.1,.05,3.93],[616.8,.06,4.15],[659.5,.05,4.48],[704.5,.041,4.84],
               [756.,.042,5.24],[821.1,.043,5.65]])
# substrate internal transmittance, from this session's bare-glass run
GLASS_TI = np.array([[350,.9226],[360,.9612],[370,.9698],[380,.9880],[400,.9952],
                     [450,1.],[500,1.],[550,1.],[600,1.],[650,.9987],[700,.9973],
                     [750,.9978],[780,.9958],[800,.9911],[850,.9920]])


def jc(l): return np.interp(l, JC[:,0], JC[:,1]) + 1j*np.interp(l, JC[:,0], JC[:,2])
def n_glass(l): return 1.5220 + 3900./l**2
def t_glass(l): return np.interp(l, GLASS_TI[:,0], GLASS_TI[:,1])


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


def observable(nk, lam, n_seed, d_ag=D_AG):
    """T and R exactly as the UMA reports them for this stack."""
    ng, ti = n_glass(lam), t_glass(lam)
    Rb = ((ng-1)/(ng+1))**2
    T1, R1 = tr([1.0, nk, n_seed, ng], [0, d_ag, D_SEED, 0], lam)
    T1b, R1b = tr([ng, n_seed, nk, 1.0], [0, D_SEED, d_ag, 0], lam)
    den = 1 - Rb*R1b*ti**2
    return T1*ti*(1-Rb)/den, R1 + KEEP*T1*ti**2*Rb*T1b/den


def invert_file(path, seed, guess):
    d = np.loadtxt(path)
    out, g = [], list(guess)
    for i in range(len(d)-1, -1, -1):                 # long lambda first
        lam, T, R = d[i,0], d[i,1]/100, d[i,2]/100
        f = lambda x: [observable(complex(abs(x[0]), abs(x[1])), lam, SEED_N[seed])[0]-T,
                       observable(complex(abs(x[0]), abs(x[1])), lam, SEED_N[seed])[1]-R]
        s = least_squares(f, g, xtol=1e-14, ftol=1e-14)
        g = list(abs(s.x))
        nk = complex(abs(s.x[0]), abs(s.x[1]))
        out.append(dict(lam=lam, T=T, R=R, A=1-T-R, nk=nk, eps=nk**2,
                        resid=float(np.max(np.abs(s.fun)))))
    return {r["lam"]: r for r in out}


def device_A(nk, lam, n_cpl=2.10, d_cpl=60.0, n_org=1.8, d_ag=D_AG):
    T, R = tr([n_org, nk, n_cpl, 1.0], [0, d_ag, d_cpl, 0], lam)
    return 1-T-R


if __name__ == "__main__":
    H = invert_file(os.path.join(DATA, "HATCN5_Ag5_TR.txt"), "HATCN", (0.5, 3.5))
    M = invert_file(os.path.join(DATA, "MoOx5_Ag5_TR.txt"), "MoOx", (2.5, 2.5))

    print(f"{'lam':>5} | {'--------- HATCN 5 / Ag 5 ---------':^34} | "
          f"{'--- MoOx 5 / Ag 5 ---':^22} | {'bulk':>7}")
    print(f"{'':>5} | {'A':>6}{'n':>7}{'k':>7}{'eps1':>8}{'A_dev':>7} | "
          f"{'A':>6}{'n':>6}{'eps1':>9} | {'eps1':>7}")
    print("-"*82)
    for l in [850,800,750,700,650,600,550,500,470,450,430,410,400,380,360]:
        if l not in H or l not in M: continue
        h, m = H[l], M[l]
        print(f"{l:5.0f} | {100*h['A']:5.1f}%{h['nk'].real:7.3f}{h['nk'].imag:7.3f}"
              f"{h['eps'].real:8.2f}{100*device_A(h['nk'], l):6.2f}% | "
              f"{100*m['A']:5.1f}%{m['nk'].real:6.2f}{m['eps'].real:9.2f} | "
              f"{(jc(l)**2).real:7.2f}")

    he = np.array([H[l]['eps'].real for l in H])
    me = np.array([M[l]['eps'].real for l in M])
    print(f"\n  HATCN eps1: {he.min():7.2f} .. {he.max():7.2f}   all negative = {he.max()<0}")
    print(f"  MoOx  eps1: {me.min():7.2f} .. {me.max():7.2f}   all positive = {me.min()>0}")

    dev = np.array([device_A(H[l]['nk'], l) for l in sorted(H) if 430 <= l <= 700])
    print(f"\n  HATCN device one-pass A, 430-700 nm, CPL 60 nm n=2.1:"
          f"  {100*dev.mean():.2f} +/- {100*dev.std():.2f} %")
    print("  MoOx: no valid slab model, so no device number can be quoted.")

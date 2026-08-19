"""Why the Ag electrode's absorption is NOT set by (n,k,d) alone, and how to get
n,k without ellipsometry.

Part A  environment sensitivity: A_Ag vs capping-layer index and thickness.
Part B  T+R -> n,k point-by-point inversion, with noise propagated to sigma(n),
        sigma(k), sigma(A).  Answers "can absolute T and R replace ellipsometry
        for a 5 nm film?"
Part C  multi-thickness global fit: one specularity p for the whole series.

Absorption in a metal is A = (w/2) * Int eps2 |E|^2 dV.  It is the field INSIDE
the metal that matters, and that field is fixed by the boundary conditions on
both sides -- hence the cap dependence.  Beer-Lambert (A ~ 4*pi*k*d/lam) assumes
the internal field equals the incident field, which is false when R is large.
"""
import numpy as np
from scipy.optimize import least_squares
import json, os

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

JC = np.array([[397.4,0.05,2.07],[413.3,0.05,2.21],[430.5,0.04,2.36],
               [450.9,0.04,2.66],[471.4,0.05,2.83],[495.9,0.05,3.09],
               [520.9,0.05,3.34],[548.6,0.06,3.59],[582.1,0.05,3.93],
               [616.8,0.06,4.15],[659.5,0.05,4.48],[704.5,0.041,4.84]])

WP, G0, HVF = 9.17, 0.021, 0.915      # eV, eV, eV*nm  (hbar*v_F for Ag)


def ag_bulk(lam):
    return np.interp(lam, JC[:, 0], JC[:, 1]) + 1j*np.interp(lam, JC[:, 0], JC[:, 2])


def ag_thin(lam, d, p=0.0):
    """Bulk J&C with the Drude term re-fitted for surface scattering."""
    hw = 1239.84/lam
    e_ib = ag_bulk(lam)**2 + WP**2/(hw**2 + 1j*G0*hw)
    g = G0 + (1.0 - p)*HVF/d
    return np.sqrt(e_ib - WP**2/(hw**2 + 1j*g*hw))


def tr(n, d, lam):
    """Coherent normal-incidence T,R.  n[0], n[-1] semi-infinite."""
    n = [complex(x) for x in n]
    k0 = 2*np.pi/lam
    M = np.eye(2, dtype=complex)
    for j in range(len(n)-1):
        r = (n[j]-n[j+1])/(n[j]+n[j+1])
        t = 2*n[j]/(n[j]+n[j+1])
        I = np.array([[1, r], [r, 1]], dtype=complex)/t
        if j+1 < len(n)-1:
            dl = k0*n[j+1]*d[j+1]
            I = I @ np.array([[np.exp(-1j*dl), 0], [0, np.exp(1j*dl)]])
        M = M @ I
    R = abs(M[1, 0]/M[0, 0])**2
    T = (n[-1].real/n[0].real)*abs(1/M[0, 0])**2
    return T, R


# ---------------------------------------------------------------- Part A
def part_A(lam=550.0, d_ag=5.0):
    nAg = ag_thin(lam, d_ag)
    print("PART A  -- the cap changes A_Ag even though n,k,d are identical")
    print(f"          Ag {d_ag:.0f} nm, lam {lam:.0f} nm, n_Ag = "
          f"{nAg.real:.3f} + {nAg.imag:.3f}i  (same in every row)\n")

    print("  semi-infinite exit medium:")
    print(f"    {'n_exit':>7} {'A_Ag':>8}")
    for ne in [1.00, 1.35, 1.50, 1.80, 2.00, 2.20, 2.50, 3.00]:
        T, R = tr([1.8, nAg, ne], [0, d_ag, 0], lam)
        print(f"    {ne:7.2f} {100*(1-T-R):7.2f} %")

    print("\n  finite capping layer, then air  (organic 1.8 / Ag / CPL / air):")
    hdr = "    " + "".join(f"{t:>9.0f}" for t in [0, 20, 40, 60, 80, 100, 120])
    print(f"    {'n_CPL':>7}" + " " + "d_CPL (nm) ->")
    print(f"    {'':>7}" + hdr[4:])
    for nc in [1.7, 1.9, 2.1, 2.3]:
        row = f"    {nc:7.2f}"
        for tc in [0, 20, 40, 60, 80, 100, 120]:
            if tc == 0:
                T, R = tr([1.8, nAg, 1.0], [0, d_ag, 0], lam)
            else:
                T, R = tr([1.8, nAg, nc, 1.0], [0, d_ag, tc, 0], lam)
            row += f"{100*(1-T-R):8.2f} "
        print(row)


# ---------------------------------------------------------------- Part B
def invert_TR(T_obs, R_obs, d_ag, lam, n_seed=1.75, d_seed=5.0, n_sub=1.52,
              guess=(0.15, 3.0)):
    """Solve for (n,k) of the Ag layer from one (T,R) pair at one wavelength."""
    def resid(x):
        nk = complex(abs(x[0]), abs(x[1]))
        T, R = tr([1.0, nk, n_seed, n_sub], [0, d_ag, d_seed, 0], lam)
        return [T - T_obs, R - R_obs]
    s = least_squares(resid, guess, xtol=1e-14, ftol=1e-14)
    return abs(s.x[0]), abs(s.x[1]), np.max(np.abs(s.fun))


def part_B(lam=550.0, d_ag=5.0, n_mc=400, sigma_TR=0.005, seed=0):
    rng = np.random.default_rng(seed)
    nk_true = ag_thin(lam, d_ag)
    T0, R0 = tr([1.0, nk_true, 1.75, 1.52], [0, d_ag, 5.0, 0], lam)
    A0 = 1 - T0 - R0
    print("\n\nPART B  -- can absolute T+R replace ellipsometry for a 5 nm film?")
    print(f"          truth: n={nk_true.real:.4f} k={nk_true.imag:.4f}  "
          f"T={T0:.4f} R={R0:.4f} A={100*A0:.2f} %")
    print(f"          Monte-Carlo: {n_mc} draws, sigma(T)=sigma(R)="
          f"{100*sigma_TR:.1f} %p absolute\n")

    ns, ks, As = [], [], []
    for _ in range(n_mc):
        Tn = T0 + rng.normal(0, sigma_TR)
        Rn = R0 + rng.normal(0, sigma_TR)
        n, k, res = invert_TR(Tn, Rn, d_ag, lam)
        if res > 1e-8:
            continue
        ns.append(n); ks.append(k); As.append(1 - Tn - Rn)
    ns, ks, As = map(np.array, (ns, ks, As))
    print(f"    n : {ns.mean():.4f} +/- {ns.std():.4f}   "
          f"({100*ns.std()/ns.mean():.1f} % rel)")
    print(f"    k : {ks.mean():.4f} +/- {ks.std():.4f}   "
          f"({100*ks.std()/ks.mean():.1f} % rel)")
    print(f"    A : {100*As.mean():.2f} +/- {100*As.std():.2f} %p  "
          "<- A needs no inversion at all, it is 1-T-R")

    # how badly does a wrong thickness hurt?
    print("\n    sensitivity to the assumed thickness (QCM error):")
    print(f"      {'d assumed':>10} {'n':>8} {'k':>8} {'k*d':>8}")
    for dd in [4.0, 4.5, 5.0, 5.5, 6.0]:
        n, k, _ = invert_TR(T0, R0, dd, lam)
        print(f"      {dd:10.1f} {n:8.4f} {k:8.4f} {k*dd:8.3f}")
    print("      -> k*d drifts +19 % over the same span in which k drifts -20 %,"
          "\n         so k*d is NOT a useful invariant here: 4*pi*k*d/lam = 0.41 at"
          "\n         5 nm, too large for the optically-thin limit.  d must come"
          "\n         from QCM/XRR to +/-0.3 nm or k carries a ~5 % systematic.")


# ---------------------------------------------------------------- Part C
def part_C(lam=550.0, sigma_TR=0.005, seed=1):
    """One free parameter (specularity p) fitted to a whole thickness series."""
    rng = np.random.default_rng(seed)
    d_list = np.array([4., 5., 6., 7., 8., 10., 12.])
    p_true = 0.25
    T_obs, R_obs = [], []
    for d in d_list:
        T, R = tr([1.0, ag_thin(lam, d, p_true), 1.75, 1.52], [0, d, 5.0, 0], lam)
        T_obs.append(T + rng.normal(0, sigma_TR))
        R_obs.append(R + rng.normal(0, sigma_TR))
    T_obs, R_obs = np.array(T_obs), np.array(R_obs)

    def resid(x):
        p = np.clip(x[0], 0.0, 0.999)
        out = []
        for d, To, Ro in zip(d_list, T_obs, R_obs):
            T, R = tr([1.0, ag_thin(lam, d, p), 1.75, 1.52], [0, d, 5.0, 0], lam)
            out += [T - To, R - Ro]
        return out

    s = least_squares(resid, [0.5], bounds=([0.], [0.999]))
    J = s.jac
    cov = np.linalg.inv(J.T @ J) * (2*np.sum(np.array(s.fun)**2)/(2*len(d_list)-1))
    print("\n\nPART C  -- global fit of the thickness series, p as the ONE unknown")
    print(f"          d = {d_list.tolist()} nm, {2*len(d_list)} data points")
    print(f"          p_true = {p_true:.3f}    p_fit = {s.x[0]:.3f} "
          f"+/- {np.sqrt(cov[0,0]):.3f}")
    print("          -> a single microstructural number, over-determined 14:1,"
          "\n             and directly comparable with the Rs series.")
    return float(s.x[0])


if __name__ == "__main__":
    part_A()
    part_B()
    p = part_C()
    os.makedirs(os.path.join(BASE, "runs"), exist_ok=True)
    with open(os.path.join(BASE, "runs", "TR_inversion_demo.json"), "w") as f:
        json.dump({"p_fit_demo": p}, f, indent=2)

"""First absolute T/R measurement inverted to n, k -- 2026-08-20.

glass / HATCN 5 nm / Ag 5 nm, film side toward the beam, Cary 6000i + UMA.
Absolute T at 0/180 deg, absolute R at 6/12 deg, SBW 8 nm.

    measured   T = 72.6 %   R = 16.0 %   at 550 nm

The bare-glass run on the same setup showed the reflection accessory clips
about an eighth of the substrate's back-surface beam (84.3 % collected), so
the measured R is corrected for that before inverting.

RESULT.  The optimistic branch is dead.  eps2 comes out 8.5x bulk against
7.3x from the sheet resistance, so essentially ALL of silver's eps2 grows with
the size effect -- the part that looks like interband absorption in Johnson &
Christy is mostly microstructural, not a fixed material constant.  A/A_bulk =
rho/rho_bulk therefore holds to about 20 %, which makes a four-point probe a
quantitative absorption meter for this film system.
"""
import numpy as np
from scipy.optimize import least_squares

L, NG, KEEP = 550.0, 1.5349, 0.843
RB = ((NG-1)/(NG+1))**2
T_OBS, R_OBS, D_AG, D_SEED, N_SEED = 0.726, 0.160, 5.0, 5.0, 1.75
RS_MEAS, RHO_BULK, G0, WP = 23.3, 1.59, 0.021, 9.17
JC_550 = 0.06 + 3.59j


def tr(n, d, lam=L):
    n = [complex(x) for x in n]; k0 = 2*np.pi/lam
    M = np.eye(2, dtype=complex)
    for j in range(len(n)-1):
        r = (n[j]-n[j+1])/(n[j]+n[j+1]); t = 2*n[j]/(n[j]+n[j+1])
        I = np.array([[1, r], [r, 1]], dtype=complex)/t
        if j+1 < len(n)-1:
            dl = k0*n[j+1]*d[j+1]
            I = I @ np.array([[np.exp(-1j*dl), 0], [0, np.exp(1j*dl)]])
        M = M @ I
    return (n[-1].real/n[0].real)*abs(1/M[0,0])**2, abs(M[1,0]/M[0,0])**2


def observable(nk, d=D_AG):
    """T and R as the instrument reports them, back-surface clipping included."""
    T1, R1 = tr([1.0, nk, N_SEED, NG], [0, d, D_SEED, 0])
    T1b, R1b = tr([NG, N_SEED, nk, 1.0], [0, D_SEED, d, 0])
    T = T1*(1-RB)/(1-RB*R1b)
    back = T1*RB*T1b/(1-RB*R1b)
    return T, R1 + KEEP*back, R1 + back


def invert(T_obs=T_OBS, R_obs=R_OBS, d=D_AG):
    def resid(x):
        T, Rm, _ = observable(complex(abs(x[0]), abs(x[1])), d)
        return [T-T_obs, Rm-R_obs]
    s = least_squares(resid, [0.3, 3.5], xtol=1e-14, ftol=1e-14)
    return complex(abs(s.x[0]), abs(s.x[1])), float(np.max(np.abs(s.fun)))


def device_A(nk, d=D_AG, n_cpl=2.10, d_cpl=60.0, n_org=1.8):
    T, R = tr([n_org, nk, n_cpl, 1.0], [0, d, d_cpl, 0])
    return 1-T-R


if __name__ == "__main__":
    nk, res = invert()
    e, e_bulk = nk**2, JC_550**2
    F_opt = e.imag/e_bulk.imag
    F_dc = RS_MEAS*D_AG*0.1/RHO_BULK

    print(f"measured  T {100*T_OBS:.1f} %   R {100*R_OBS:.1f} %   "
          f"A {100*(1-T_OBS-R_OBS):.1f} %      residual {100*res:.4f} %p\n")
    print(f"  n = {nk.real:.3f}   k = {nk.imag:.3f}"
          f"      eps1 = {e.real:.2f}   eps2 = {e.imag:.3f}")
    print(f"  eps2 / eps2_bulk (optical) = {F_opt:.2f}")
    print(f"  rho  / rho_bulk  (DC)      = {F_dc:.2f}")
    print(f"  ratio optical/DC           = {F_opt/F_dc:.2f}\n")

    print("  the two rival models, for comparison:")
    hw = 1239.84/L
    ed = (-WP**2/(hw**2 + 1j*G0*hw)).imag
    print(f"    only the Drude part scales -> eps2 = "
          f"{ed*F_dc + (e_bulk.imag-ed):.2f}   (measured {e.imag:.2f})")
    print(f"    all of eps2 scales         -> eps2 = "
          f"{e_bulk.imag*F_dc:.2f}   (measured {e.imag:.2f})  <-- this one\n")

    print("  device one-pass absorption with the measured n,k:")
    for lab, nc, dc in [("no cap", 1.0, 0.0), ("CPL 60 nm n=2.1", 2.10, 60.0),
                        ("CPL 65 nm n=2.1", 2.10, 65.0), ("CPL 60 nm n=2.3", 2.30, 60.0)]:
        if dc == 0:
            T, R = tr([1.8, nk, 1.0], [0, D_AG, 0])
            a = 1-T-R
        else:
            a = device_A(nk, n_cpl=nc, d_cpl=dc)
        print(f"    {lab:<20}{100*a:7.2f} %")

    print("\n  sensitivity to the assumed Ag thickness:")
    for d in [4.5, 4.75, 5.0, 5.25, 5.5]:
        nk_d, _ = invert(d=d)
        print(f"    d = {d:.2f} nm ->  n {nk_d.real:.3f}  k {nk_d.imag:.3f}"
              f"   eps2/bulk {(nk_d**2).imag/e_bulk.imag:5.2f}"
              f"   device A {100*device_A(nk_d, d=d):.2f} %")

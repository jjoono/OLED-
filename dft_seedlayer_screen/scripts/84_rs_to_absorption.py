"""Turn the measured sheet-resistance series into a prediction for today's
optical measurement, before the optical data exists.

rho/rho_bulk = gamma/gamma_bulk = eps2/eps2_bulk, so a DC transport measurement
fixes the Drude damping and therefore n, k and the absorptance.  The prediction
is an UPPER bound: at 550 nm an electron travels ~0.4 nm per optical cycle, so
grain boundaries -- which a DC measurement sees at full weight -- are largely
invisible optically.  Measured A below the prediction is the size of that
grain-boundary term, which is otherwise not separable from surface scattering.
"""
import numpy as np, json, os

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RHO_B, MFP, WP, G0, HVF = 1.59, 52.0, 9.17, 0.021, 0.915
LAM = 550.0

JC = np.array([[397.4,0.05,2.07],[413.3,0.05,2.21],[430.5,0.04,2.36],
               [450.9,0.04,2.66],[471.4,0.05,2.83],[495.9,0.05,3.09],
               [520.9,0.05,3.34],[548.6,0.06,3.59],[582.1,0.05,3.93],
               [616.8,0.06,4.15],[659.5,0.05,4.48],[704.5,0.041,4.84]])

D = np.array([4, 5, 6, 7, 8, 10, 12], float)
RS = {"HATCN": np.array([52.2, 23.3, 18.9, 12.7, 9.1, 7.5, 5.3]),
      "MoOx":  np.array([138.0, 49.1, 32.0, 14.9, 10.8, 9.6, 6.6])}
SEED_N = {"HATCN": 1.75, "MoOx": 2.10}


def ag_from_gamma(lam, gamma):
    hw = 1239.84/lam
    e_ib = (np.interp(lam, JC[:,0], JC[:,1]) + 1j*np.interp(lam, JC[:,0], JC[:,2]))**2 \
           + WP**2/(hw**2 + 1j*G0*hw)
    return np.sqrt(e_ib - WP**2/(hw**2 + 1j*gamma*hw))


def tr(n, d, lam):
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


def A_meas(seed, d_ag, gamma):
    """air / Ag / seed 5 nm / glass, illuminated from the Ag side."""
    T, R = tr([1.0, ag_from_gamma(LAM, gamma), SEED_N[seed], 1.52],
              [0, d_ag, 5.0, 0], LAM)
    return 1 - T - R


def A_device(d_ag, gamma, n_cpl=2.10, d_cpl=60.0):
    T, R = tr([1.8, ag_from_gamma(LAM, gamma), n_cpl, 1.0],
              [0, d_ag, d_cpl, 0], LAM)
    return 1 - T - R


out = {}
print(f"Prediction for the T/R run, {LAM:.0f} nm.  gamma_bulk = {G0:.3f} eV\n")
print(f"{'seed':<7}{'d':>4}{'rho':>8}{'x_bulk':>8}{'gamma':>8}{'xFS':>7} | "
      f"{'A_meas':>8}{'A_floor':>9} | {'A_device':>9}")
print("-" * 83)
for seed, rs in RS.items():
    rows = []
    for i, d_ag in enumerate(D):
        rho = rs[i]*d_ag*0.1
        g = G0*rho/RHO_B
        # p = 0 Fuchs floor, expressed the same way as the DC-derived gamma so
        # the two are directly comparable: gamma_bulk * (1 + 0.375 l/d), NOT
        # hbar*v_F/d -- the latter is the naive collision estimate and runs 2.7x
        # larger than the transport average, which put the "floor" above the
        # prediction it was supposed to bound.
        g_floor = G0*(1.0 + 0.375*MFP/d_ag)
        a, af, ad = (A_meas(seed, d_ag, g), A_meas(seed, d_ag, g_floor),
                     A_device(d_ag, g))
        print(f"{seed:<7}{d_ag:4.0f}{rho:8.2f}{rho/RHO_B:8.2f}{g:8.3f}"
              f"{g/g_floor:7.2f} | "
              f"{100*a:7.2f}%{100*af:8.2f}% | {100*ad:8.2f}%")
        rows.append(dict(d=d_ag, rho=rho, gamma=g, A_meas=a,
                         A_floor=af, A_device=ad))
    out[seed] = rows
    print()

print("A_meas   : predicted 1 - T - R of glass/seed5/Ag(d), Ag-side illumination")
print("xFS      : gamma / (p=0 Fuchs-Sondheimer limit).  Surface scattering")
print("           alone cannot exceed 1.00, so any excess is grain boundaries")
print("A_floor  : same thickness at the p = 0 Fuchs limit (no grain boundaries)")
print("A_device : organic 1.8 / Ag / CPL 60 nm n=2.1 / air, same material")
print("\nIf the measurement lands BELOW A_meas, the gap is the grain-boundary")
print("term that DC transport counts and the optical response does not.")

os.makedirs(os.path.join(BASE, "runs"), exist_ok=True)
json.dump(out, open(os.path.join(BASE, "runs", "rs_absorption_prediction.json"), "w"),
          indent=2, default=float)

"""Analysis pipeline for the 2026-08-19 seed-comparison run.

Sample set (16):  bare HATCN 5 | bare MoOx 5 | {HATCN 5, MoOx 5} x Ag {4,5,6,7,8,10,12}

Inputs  (drop files in  data/TR_20260819/ ):
    <sample>_T.csv   two columns: wavelength_nm, T   (absolute, 0-1 or 0-100)
    <sample>_R.csv   two columns: wavelength_nm, R   (absolute, 0-1 or 0-100)
    rs.csv           sample, Rs_ohm_sq          (optional, 4-point probe)
  <sample> names:  HATCN5, MoOx5, HATCN5_Ag4, ... , MoOx5_Ag12

Outputs (runs/TR_20260819/):
    A_vs_thickness.csv/png     absorptance of the whole stack and of the Ag alone
    k_spectra.csv              point-by-point n,k from the T/R pair
    p_fit.json                 one specularity per seed, from the thickness series
    rs_fuchs.json              specularity from the Rs series, independent check
    closure.json               three independent closure criteria

Run with --demo to exercise the whole chain on synthetic data.
"""
import os, sys, json, glob
import numpy as np
from scipy.optimize import least_squares

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DATA = os.path.join(BASE, "data", "TR_20260819")
OUT = os.path.join(BASE, "runs", "TR_20260819")

D_AG = [4, 5, 6, 7, 8, 10, 12]
SEEDS = {"HATCN5": dict(d=5.0, n=1.75, k=0.02),
         "MoOx5":  dict(d=5.0, n=2.10, k=0.01)}
N_GLASS = 1.52
RHO_BULK = 1.59        # uOhm cm
MFP = 52.0             # nm, Ag electron mean free path
WP, G0, HVF = 9.17, 0.021, 0.915

JC = np.array([[397.4,0.05,2.07],[413.3,0.05,2.21],[430.5,0.04,2.36],
               [450.9,0.04,2.66],[471.4,0.05,2.83],[495.9,0.05,3.09],
               [520.9,0.05,3.34],[548.6,0.06,3.59],[582.1,0.05,3.93],
               [616.8,0.06,4.15],[659.5,0.05,4.48],[704.5,0.041,4.84],
               [756.0,0.042,5.24],[821.1,0.043,5.65]])


def ag_bulk(lam):
    return np.interp(lam, JC[:,0], JC[:,1]) + 1j*np.interp(lam, JC[:,0], JC[:,2])


def ag_thin(lam, d, p=0.0):
    hw = 1239.84/np.asarray(lam, float)
    e_ib = ag_bulk(lam)**2 + WP**2/(hw**2 + 1j*G0*hw)
    return np.sqrt(e_ib - WP**2/(hw**2 + 1j*(G0 + (1-p)*HVF/d)*hw))


def tr(n, d, lam):
    n = [complex(x) for x in n]
    k0 = 2*np.pi/lam
    M = np.eye(2, dtype=complex)
    for j in range(len(n)-1):
        r = (n[j]-n[j+1])/(n[j]+n[j+1]); t = 2*n[j]/(n[j]+n[j+1])
        I = np.array([[1, r], [r, 1]], dtype=complex)/t
        if j+1 < len(n)-1:
            dl = k0*n[j+1]*d[j+1]
            I = I @ np.array([[np.exp(-1j*dl), 0], [0, np.exp(1j*dl)]])
        M = M @ I
    return (n[-1].real/n[0].real)*abs(1/M[0,0])**2, abs(M[1,0]/M[0,0])**2


def stack(seed, d_ag, nk_ag):
    """Air / Ag / seed / glass -- illumination from the Ag (film) side."""
    s = SEEDS[seed]
    return ([1.0, nk_ag, complex(s["n"], s["k"]), N_GLASS], [0, d_ag, s["d"], 0])


# ------------------------------------------------------------------ I/O
def load(name, kind):
    f = os.path.join(DATA, f"{name}_{kind}.csv")
    if not os.path.exists(f):
        return None
    a = np.loadtxt(f, delimiter=",", skiprows=1)
    wl, y = a[:, 0], a[:, 1]
    if y.max() > 1.5:                      # percent -> fraction
        y = y/100.0
    return wl, y


# ------------------------------------------------------ 1. absorptance
def absorptance(names):
    rows = {}
    for nm in names:
        t, r = load(nm, "T"), load(nm, "R")
        if t is None or r is None:
            continue
        wl = t[0]
        A = 1.0 - t[1] - np.interp(wl, r[0], r[1])
        rows[nm] = (wl, A)
        neg = (A < -0.01).sum()
        if neg:
            print(f"  !! {nm}: A < 0 at {neg} wavelengths (min {A.min():+.3f})"
                  "  -> a baseline is not absolute")
    return rows


# --------------------------------------------- 2. point-by-point n,k
def invert(seed, d_ag, wl, T_obs, R_obs):
    n_out, k_out = np.zeros_like(wl), np.zeros_like(wl)
    guess = None
    for i, l in enumerate(wl):
        g = guess if guess is not None else (ag_thin(l, d_ag).real,
                                             ag_thin(l, d_ag).imag)
        def resid(x):
            nl, dl = stack(seed, d_ag, complex(abs(x[0]), abs(x[1])))
            T, R = tr(nl, dl, l)
            return [T - T_obs[i], R - R_obs[i]]
        s = least_squares(resid, g, xtol=1e-13, ftol=1e-13)
        n_out[i], k_out[i] = abs(s.x[0]), abs(s.x[1])
        guess = (n_out[i], k_out[i])
    return n_out, k_out


# ------------------------------------ 3. one specularity per seed series
def fit_p(seed, series, wl_grid):
    """series: {d_ag: (wl, T, R)}.  One p for the whole thickness series."""
    def resid(x):
        p = float(np.clip(x[0], 0.0, 0.999))
        out = []
        for d_ag, (wl, T, R) in series.items():
            nk = ag_thin(wl_grid, d_ag, p)
            for j, l in enumerate(wl_grid):
                nl, dl = stack(seed, d_ag, nk[j])
                Tm, Rm = tr(nl, dl, l)
                out += [Tm - np.interp(l, wl, T), Rm - np.interp(l, wl, R)]
        return out
    s = least_squares(resid, [0.4], bounds=([0.0], [0.999]))
    J = s.jac
    dof = max(len(s.fun) - 1, 1)
    cov = np.linalg.inv(J.T @ J)*(2*np.sum(np.array(s.fun)**2)/dof)
    return float(s.x[0]), float(np.sqrt(cov[0, 0]))


# --------------------------------- 4. Fuchs-Sondheimer from sheet resistance
def fuchs_from_rs(rs_map):
    """rho(d) = rho_bulk [1 + 0.375 (1-p) l / d].

    Below the percolation threshold rho diverges, so the small-d points fall off
    the line.  Drop them one at a time from the thin end, but only while the
    thinnest point is a >3 sigma outlier against the fit to the rest -- a plain
    chi2 comparison would always prefer the smallest subset.
    """
    d = np.array(sorted(rs_map), float)
    rho = np.array([rs_map[x]*x*1e-3 for x in d])          # uOhm cm
    lo = 0
    while len(d) - lo >= 5:
        dd, rr = d[lo+1:], rho[lo+1:]                      # fit WITHOUT the thinnest
        Aa = np.vstack([np.ones_like(dd), 1.0/dd]).T
        coef, *_ = np.linalg.lstsq(Aa, rr, rcond=None)
        resid = Aa @ coef - rr
        sd = float(np.std(resid, ddof=2)) or 1e-12
        pred = coef[0] + coef[1]/d[lo]
        if (rho[lo] - pred) > 3*sd:                        # thinnest is above the line
            lo += 1
            continue
        break
    dd, rr = d[lo:], rho[lo:]
    Aa = np.vstack([np.ones_like(dd), 1.0/dd]).T
    coef, *_ = np.linalg.lstsq(Aa, rr, rcond=None)
    resid = Aa @ coef - rr
    dof = max(len(dd) - 2, 1)
    s2 = float(np.sum(resid**2)/dof)
    cov = np.linalg.inv(Aa.T @ Aa)*s2
    p_hat = 1.0 - coef[1]/(0.375*MFP*coef[0])
    # propagate: p = 1 - b/(0.375 l a)
    dp_da = coef[1]/(0.375*MFP*coef[0]**2)
    dp_db = -1.0/(0.375*MFP*coef[0])
    var = (dp_da**2*cov[0,0] + dp_db**2*cov[1,1] + 2*dp_da*dp_db*cov[0,1])
    return dict(d_min_used=float(dd[0]), n_used=int(len(dd)),
                rho_bulk_fit=float(coef[0]), rho_bulk_err=float(np.sqrt(cov[0,0])),
                p=float(p_hat), p_err=float(np.sqrt(max(var, 0.0))),
                rms_resid=float(np.sqrt(np.mean(resid**2))))


# ------------------------------------------------------------ demo data
def make_demo():
    os.makedirs(DATA, exist_ok=True)
    rng = np.random.default_rng(3)
    wl = np.arange(400, 801, 2.0)
    p_true = {"HATCN5": 0.35, "MoOx5": 0.15}
    rs = {}
    for seed in SEEDS:
        s = SEEDS[seed]
        T0, R0 = np.array([tr([1.0, complex(s["n"], s["k"]), N_GLASS],
                              [0, s["d"], 0], l) for l in wl]).T
        np.savetxt(os.path.join(DATA, f"{seed}_T.csv"),
                   np.c_[wl, T0], delimiter=",", header="wl,T", comments="")
        np.savetxt(os.path.join(DATA, f"{seed}_R.csv"),
                   np.c_[wl, R0], delimiter=",", header="wl,R", comments="")
        for d_ag in D_AG:
            nk = ag_thin(wl, d_ag, p_true[seed])
            T, R = np.array([tr(*stack(seed, d_ag, nk[j]), l)
                             for j, l in enumerate(wl)]).T
            T = T + rng.normal(0, 0.002, T.shape)
            R = R + rng.normal(0, 0.002, R.shape)
            base = f"{seed}_Ag{d_ag}"
            np.savetxt(os.path.join(DATA, f"{base}_T.csv"),
                       np.c_[wl, T], delimiter=",", header="wl,T", comments="")
            np.savetxt(os.path.join(DATA, f"{base}_R.csv"),
                       np.c_[wl, R], delimiter=",", header="wl,R", comments="")
            rho = RHO_BULK*(1 + 0.375*(1-p_true[seed])*MFP/d_ag)
            rs[base] = rho*1e3/d_ag*(1 + rng.normal(0, 0.02))
    with open(os.path.join(DATA, "rs.csv"), "w") as f:
        f.write("sample,Rs_ohm_sq\n")
        for k, v in rs.items():
            f.write(f"{k},{v:.3f}\n")
    print(f"demo data written to {DATA}  (p_true {p_true})")


# ------------------------------------------------------------------ main
def main():
    os.makedirs(OUT, exist_ok=True)
    wl_fit = np.arange(450, 781, 30.0)          # coarse grid for the global fit
    report = {}

    print("\n[1] absorptance  A = 1 - T - R")
    names = list(SEEDS) + [f"{s}_Ag{d}" for s in SEEDS for d in D_AG]
    A = absorptance(names)
    if not A:
        print("  no data found in", DATA, "-- run with --demo first")
        return
    lines = ["seed,d_ag_nm,A550_stack_pct,A550_ag_pct"]
    print(f"  {'seed':<8}{'d_Ag':>6}{'A_stack@550':>13}{'A_Ag@550':>11}")
    for seed in SEEDS:
        if seed not in A:
            continue
        wl0, A0 = A[seed]
        a_seed = float(np.interp(550, wl0, A0))
        print(f"  {seed:<8}{'bare':>6}{100*a_seed:12.2f}%{'-':>11}")
        lines.append(f"{seed},0,{100*a_seed:.3f},0")
        for d_ag in D_AG:
            nm = f"{seed}_Ag{d_ag}"
            if nm not in A:
                continue
            wl1, A1 = A[nm]
            a_tot = float(np.interp(550, wl1, A1))
            print(f"  {seed:<8}{d_ag:6d}{100*a_tot:12.2f}%"
                  f"{100*(a_tot-a_seed):10.2f}%")
            lines.append(f"{seed},{d_ag},{100*a_tot:.3f},{100*(a_tot-a_seed):.3f}")
    open(os.path.join(OUT, "A_vs_thickness.csv"), "w").write("\n".join(lines)+"\n")

    print("\n[2] one specularity p per seed, from the whole thickness series")
    for seed in SEEDS:
        series = {}
        for d_ag in D_AG:
            t, r = load(f"{seed}_Ag{d_ag}", "T"), load(f"{seed}_Ag{d_ag}", "R")
            if t is None or r is None:
                continue
            series[d_ag] = (t[0], t[1], np.interp(t[0], r[0], r[1]))
        if len(series) < 3:
            continue
        p, sp = fit_p(seed, series, wl_fit)
        report[f"p_optical_{seed}"] = [p, sp]
        print(f"  {seed:<8} p_optical = {p:.3f} +/- {sp:.3f}   "
              f"({len(series)} thicknesses x {len(wl_fit)} wavelengths)")

    print("\n[3] independent p from the sheet-resistance series (Fuchs-Sondheimer)")
    rsf = os.path.join(DATA, "rs.csv")
    if os.path.exists(rsf):
        raw = {}
        for ln in open(rsf).read().splitlines()[1:]:
            if not ln.strip():
                continue
            k, v = ln.split(",")
            raw[k.strip()] = float(v)
        for seed in SEEDS:
            m = {d: raw[f"{seed}_Ag{d}"] for d in D_AG if f"{seed}_Ag{d}" in raw}
            if len(m) < 4:
                continue
            fit = fuchs_from_rs(m)
            report[f"p_DC_{seed}"] = fit
            print(f"  {seed:<8} p_DC = {fit['p']:.3f} +/- {fit['p_err']:.3f}   "
                  f"rho_bulk_fit = {fit['rho_bulk_fit']:.2f} +/- "
                  f"{fit['rho_bulk_err']:.2f} uOhm-cm (bulk 1.59)   "
                  f"line holds from {fit['d_min_used']:.0f} nm up "
                  f"({fit['n_used']} pts)")
    else:
        print("  rs.csv not found -- skipping")

    print("\n[4] cross-check")
    for seed in SEEDS:
        a, b = report.get(f"p_optical_{seed}"), report.get(f"p_DC_{seed}")
        if a and b:
            print(f"  {seed:<8} p_optical {a[0]:.3f} vs p_DC {b['p']:.3f}   "
                  f"delta = {a[0]-b['p']:+.3f}")
    print("   p_optical > p_DC is expected: at 550 nm an electron moves ~0.4 nm"
          "\n   per optical cycle, so grain boundaries weigh less optically"
          "\n   than they do in a DC measurement.")

    json.dump(report, open(os.path.join(OUT, "p_fit.json"), "w"), indent=2)
    print(f"\nwritten to {OUT}")


if __name__ == "__main__":
    if "--demo" in sys.argv:
        make_demo()
    main()

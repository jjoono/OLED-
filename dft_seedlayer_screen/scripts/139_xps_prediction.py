"""What should the XPS see, if the charge-state mechanism is right?

    python scripts/139_xps_prediction.py

The simulation's one surviving claim is about charge: silver gives an electron
to MoO3 and keeps giving as clusters grow, while on HATCN the charge comes back
once a second atom arrives. The Ag-side XPS signature of that (Auger parameter)
is confounded by the final-state size effect of small clusters, so the clean
test is on the substrate side: every electron silver gives to MoO3 turns a Mo6+
into Mo5+, and Mo 3d resolves the two by ~1.2-1.5 eV.

This turns the computed charges into the number an XPS fit returns -- the Mo5+
fraction of the Mo 3d signal, *added* on top of whatever the pristine film
already has -- as a function of Ag coverage, so the experiment can be judged
against a prediction written down before it was run.

Model, and every assumption in it:

  * electrons per Ag atom: bracketed by the dispersed limit (q of a lone adatom)
    and the clustered limit (q per atom of Ag4), each taken as the range across
    charge schemes (Mulliken, Lowdin, IAO, NPA, Hirshfeld) where script 138 has
    finished, else package G's Mulliken values
  * each electron makes one Mo5+ in the topmost Mo plane, spilling to the next
    plane only when that one is saturated (upper bound on surface localisation,
    which is what DFT shows for MoO3(010) polarons)
  * alpha-MoO3 (010) stacking: Mo planes of 6.83 /nm2, two per bilayer 2.82 A
    apart, bilayers every 6.93 A; an amorphous evaporated film is less dense,
    which makes the signal *larger* for the same transfer
  * signal weighted by exp(-z / (lambda cos theta)); lambda from TPP-2M
  * Ag overlayer attenuates all Mo planes equally to first order, so it drops
    out of the fraction

and the same construction for HATCN: each electron makes one HATCN anion, whose
N 1s shifts to lower binding energy; flat-lying molecules at 0.93 /nm2 per
0.34 nm layer.

What it cannot say: absolute charges are scheme-dependent (that is why the
range is carried), Mo5+ polarons can delocalise over more planes (lowers the
predicted fraction at normal emission), and a substoichiometric film starts
with 5-15% Mo5+, so the observable is the change from the pristine baseline
measured on the same sample before silver.
"""
import os, json

import numpy as np

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
RUNS = os.path.join(ROOT, "runs")
AG_ML = 13.84            # Ag atoms per nm2 in one Ag(111) monolayer
HV = 1486.6              # Al K-alpha


def tpp2m(E, rho, M, Nv, Eg):
    """Inelastic mean free path in Angstrom (Tanuma, Powell, Penn TPP-2M)."""
    Ep = 28.8 * np.sqrt(Nv * rho / M)
    U = Nv * rho / M
    beta = -0.10 + 0.944 / np.sqrt(Ep ** 2 + Eg ** 2) + 0.069 * rho ** 0.1
    gamma = 0.191 * rho ** -0.5
    C = 1.97 - 0.91 * U
    D = 53.4 - 20.8 * U
    return E / (Ep ** 2 * (beta * np.log(gamma * E) - C / E + D / E ** 2))


def charges():
    """(dispersed, clustered) electrons per Ag, each as (lo, hi) across schemes."""
    loc = os.path.join(RUNS, "charge_robustness_local.json")
    cg = json.load(open(os.path.join(RUNS, "cluster_growth.json")))
    out, src = {}, {}
    for t in ("Mo3O9", "HATCN"):
        d1, d4 = [], []
        if os.path.exists(loc):
            L = json.load(open(loc))
            for s in ("mulliken", "lowdin", "iao", "npa", "hirshfeld"):
                a = L.get(f"{t}_Ag1_svp", {}).get(s)
                b = L.get(f"{t}_Ag4_svp", {}).get(s)
                if a is not None:
                    d1.append(a)
                if b is not None:
                    d4.append(b / 4)
        src[t] = "schemes" if d1 and d4 else "Mulliken (package G)"
        if not d1:
            d1 = [cg[t]["q1"]]
        if not d4:
            d4 = [cg[t]["steps"]["4"]["q_total"] / 4]
        out[t] = ((min(d1), max(d1)), (min(d4), max(d4)))
    return out, src


def planes_moo3(n=12):
    z, k = [], 0
    while len(z) < n:
        z += [6.93 * k, 6.93 * k + 2.82]
        k += 1
    return np.array(z[:n]), np.full(n, 6.83)


def planes_hatcn(n=12):
    return 3.4 * np.arange(n), np.full(n, 0.93)


def reduced_fraction(ne, z, dens, lam, theta_deg):
    """Fraction of the core-level signal from reduced centres, filling the
    topmost plane first."""
    w = np.exp(-z / (lam * np.cos(np.radians(theta_deg))))
    left = ne
    red = np.zeros_like(dens)
    for i in range(len(dens)):
        red[i] = min(left, dens[i])
        left -= red[i]
    return float((w * red).sum() / (w * dens).sum())


def electrons_per_ag(frac, ag_ml, substrate, theta_deg=0.0):
    """Invert the model: a fitted reduced fraction at a measured Ag coverage
    (Ag ML from the Ag 3d / substrate intensity ratio) -> electrons per Ag.

    This, not the reduced fraction itself, is what can be compared between
    substrates: HATCN offers 0.93 acceptor molecules per nm2 against 6.83 Mo
    per plane, so the same transfer per silver atom reduces a far larger share
    of the HATCN signal.
    """
    lam = (tpp2m(HV - 232.5, 4.69, 143.94, 24, 3.0) if substrate == "MoO3"
           else tpp2m(HV - 399.0, 1.60, 384.2, 132, 3.3))
    z, dens = planes_moo3() if substrate == "MoO3" else planes_hatcn()
    lo, hi = 0.0, 3.0 * ag_ml * AG_ML
    for _ in range(60):
        mid = (lo + hi) / 2
        if reduced_fraction(mid, z, dens, lam, theta_deg) < frac:
            lo = mid
        else:
            hi = mid
    return lo / (ag_ml * AG_ML)


def main():
    lam_mo = tpp2m(HV - 232.5, 4.69, 143.94, 24, 3.0)      # Mo 3d5/2 in MoO3
    lam_n = tpp2m(HV - 399.0, 1.60, 384.2, 132, 3.3)        # N 1s in HATCN
    q, src = charges()
    print(f"IMFP (TPP-2M): Mo 3d in MoO3 {lam_mo:.1f} A, N 1s in HATCN {lam_n:.1f} A")
    for t in q:
        (a, b), (c, d) = q[t]
        print(f"  {t}: e per Ag  dispersed {a:+.2f}..{b:+.2f}  "
              f"clustered {c:+.2f}..{d:+.2f}   ({src[t]})")
    cov = [0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0]
    res = {"imfp_A": {"Mo3d_MoO3": round(lam_mo, 2), "N1s_HATCN": round(lam_n, 2)},
           "charges": {t: {"dispersed": q[t][0], "clustered": q[t][1],
                           "source": src[t]} for t in q},
           "prediction": {}}
    for t, (z, dens), lam, label in (
            ("Mo3O9", planes_moo3(), lam_mo, "added Mo5+ fraction of Mo 3d"),
            ("HATCN", planes_hatcn(), lam_n, "HATCN anion fraction of N 1s")):
        print(f"\n{label}  (range = scheme spread x cluster-state bound)")
        print(f"{'Ag (ML)':>8}{'Ag (nm)':>9}" + "".join(f"{'th=' + str(a) + 'deg':>20}" for a in (0, 60)))
        rows = []
        for c in cov:
            n_ag = c * AG_ML
            lo = min(q[t][0][0], q[t][1][0]) * n_ag
            hi = max(q[t][0][1], q[t][1][1]) * n_ag
            cells = []
            for ang in (0, 60):
                f_lo = reduced_fraction(max(lo, 0), z, dens, lam, ang)
                f_hi = reduced_fraction(max(hi, 0), z, dens, lam, ang)
                cells.append((f_lo, f_hi))
            rows.append({"ML": c, "nm": round(c * 0.236, 3),
                         "normal": [round(x, 4) for x in cells[0]],
                         "grazing60": [round(x, 4) for x in cells[1]]})
            print(f"{c:>8.2f}{c * 0.236:>9.3f}" + "".join(
                f"{100 * a:>9.1f}-{100 * b:<5.1f}%    " for a, b in cells))
        res["prediction"][t] = rows
    json.dump(res, open(os.path.join(RUNS, "xps_prediction.json"), "w"), indent=1)
    print("\nREAD THIS BEFORE COMPARING SUBSTRATES: at equal Ag coverage HATCN shows")
    print("the LARGER reduced fraction even though each Ag gives it less charge,")
    print("because HATCN has 0.93 acceptor molecules/nm2 against 6.83 Mo per plane.")
    print("Compare electrons per Ag atom, from electrons_per_ag(), never raw fractions.")
    print("Discriminating prediction (clustered regime, >~0.3 ML):")
    for t in q:
        print(f"  {t:<6} {q[t][1][0]:.2f}-{q[t][1][1]:.2f} e per Ag")
    print("round trip check: 5% Mo5+ at 0.3 ML normal emission ->",
          f"{electrons_per_ag(0.05, 0.3, 'MoO3'):.3f} e/Ag")
    print("\nwrote runs/xps_prediction.json")


if __name__ == "__main__":
    main()

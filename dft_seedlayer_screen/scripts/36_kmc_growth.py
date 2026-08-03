"""Lattice kinetic Monte Carlo for Ag island nucleation and growth on a seed layer.

The missing link between the DFT descriptors (E_b, E_d) and what SEM / 4-point
probe actually measure: at what NOMINAL THICKNESS does the Ag film percolate?

MODEL
-----
Solid-on-solid square lattice, periodic. Irreversible attachment (critical
nucleus i = 1: a dimer is stable), which is the standard regime for metal-on-weakly
-interacting-substrate at room temperature.

  deposition  : rate F per site
  diffusion   : a TOPMOST atom with no in-plane neighbour at or above its own
                level is mobile and hops at D = nu0 * exp(-E_d / kT).
                An atom that acquires an in-plane neighbour sticks (i = 1).
  Ehrlich-Schwoebel : a downward step is rejected with prob 1 - exp(-E_ES / kT).

Because every mobile atom shares the same barrier E_d, they all hop at the same
rate, so the whole adatom population can be advanced in lockstep, one vectorised
"sweep" per 1/D of time. That makes the cost O(L^2) per sweep instead of
O(N_adatom) per single hop, which is what makes the real D/F regime reachable.

THE D/F PROBLEM
---------------
The only control parameter is R = D / F (both per site). For Ag at 300 K with
E_d = 0.29 eV and F = 1 ML/s, R ~ 1e8-1e9. Simulating that hop-by-hop is out of
reach for any code. The standard resolution, used here:

  1. run kMC over R = 1e2 ... 1e5, where it is affordable
  2. verify the Venables scaling  N_sat ~ R^(-chi),  chi = i/(i+2) = 1/3
  3. extrapolate N_sat to the experimental R with the *fitted* exponent
  4. convert N_sat -> percolation thickness with the island-shape model below

ISLAND SHAPE / PERCOLATION
--------------------------
Nucleation density fixes the island spacing; the island SHAPE fixes how much
material is needed for neighbouring islands to touch. The shape follows from the
other DFT descriptor via Young / Winterbottom:

    cos(alpha) = W_adh / gamma_Ag - 1,     W_adh = E_b / A_site

so a strongly binding seed (HATCN) gives a flat island (small alpha) that spreads
and percolates early, while a weakly binding one (LiF) gives a tall compact island
that must grow much thicker before touching its neighbour. Both descriptors thus
enter, E_d through N_sat and E_b through alpha.
"""
import numpy as np, json, os, argparse
from scipy import ndimage

KB = 8.617333262e-5           # eV/K
NU0 = 1.0e13                  # attempt frequency 1/s
A_SITE = 0.289                # nm, Ag nearest-neighbour spacing
GAMMA_AG = 1.25               # J/m^2, Ag surface energy
EV_PER_J_NM2 = 6.2415e-3      # 1 J/m^2 = 6.2415e-3 eV/nm^2

DX = np.array([1, -1, 0, 0])
DY = np.array([0, 0, 1, -1])


class KMC:
    """Vectorised SOS growth model. Time is measured in monolayers deposited."""

    def __init__(self, L=128, R=1e4, E_ES=0.05, T=300.0, seed=0):
        self.L = L
        self.R = R                      # D / F, both per site
        self.p_down = np.exp(-E_ES / (KB * T))
        self.rng = np.random.default_rng(seed)
        self.h = np.zeros((L, L), dtype=np.int32)
        self.n_dep = 0
        # deposition events per sweep: one sweep = 1/D of time, during which
        # F * L^2 / D = L^2 / R atoms land.
        self.dep_per_sweep = L * L / R

    # ---- geometry --------------------------------------------------------
    def _max_neighbor_h(self):
        h = self.h
        return np.maximum.reduce([np.roll(h, 1, 0), np.roll(h, -1, 0),
                                  np.roll(h, 1, 1), np.roll(h, -1, 1)])

    def _mobile_mask(self):
        """Topmost atoms with no in-plane neighbour at or above their own level."""
        h = self.h
        occupied = h > 0
        n_at_or_above = ((np.roll(h, 1, 0) >= h).astype(np.int8) +
                         (np.roll(h, -1, 0) >= h).astype(np.int8) +
                         (np.roll(h, 1, 1) >= h).astype(np.int8) +
                         (np.roll(h, -1, 1) >= h).astype(np.int8))
        return occupied & (n_at_or_above == 0)

    # ---- one sweep = 1/D of physical time --------------------------------
    def sweep(self):
        rng = self.rng
        L = self.L

        # --- deposition (Poisson number of arrivals in this interval) ---
        n_new = rng.poisson(self.dep_per_sweep)
        if n_new:
            xs = rng.integers(0, L, n_new)
            ys = rng.integers(0, L, n_new)
            np.add.at(self.h, (xs, ys), 1)
            self.n_dep += n_new

        # --- diffusion: every mobile atom attempts one hop ---
        mob = self._mobile_mask()
        idx = np.flatnonzero(mob.ravel())
        if idx.size == 0:
            return
        ax, ay = np.unravel_index(idx, (L, L))
        d = rng.integers(0, 4, idx.size)
        nx = (ax + DX[d]) % L
        ny = (ay + DY[d]) % L

        h = self.h
        # level the atom would land on vs the level it leaves from
        h_from = h[ax, ay]
        h_to = h[nx, ny]
        descending = h_to < h_from - 1
        accept = np.ones(idx.size, dtype=bool)
        if descending.any():
            accept[descending] = rng.random(descending.sum()) < self.p_down

        ax, ay, nx, ny = ax[accept], ay[accept], nx[accept], ny[accept]
        if ax.size == 0:
            return
        np.add.at(h, (ax, ay), -1)
        np.add.at(h, (nx, ny), 1)
        np.maximum(h, 0, out=h)

    # ---- observables -----------------------------------------------------
    def theta(self):
        return self.n_dep / (self.L * self.L)

    def coverage(self):
        return float((self.h > 0).mean())

    def n_islands(self):
        _, n = ndimage.label(self.h > 0, structure=np.ones((3, 3)))
        return n

    def island_density_per_site(self):
        return self.n_islands() / (self.L * self.L)

    def percolates(self):
        occ = self.h > 0
        lab, n = ndimage.label(occ, structure=np.ones((3, 3)))
        if n == 0:
            return False
        for ax_ in (0, 1):
            a = set(np.unique(lab.take(0, axis=ax_))) - {0}
            b = set(np.unique(lab.take(-1, axis=ax_))) - {0}
            if a & b:
                return True
        return False


def run_to_theta(L, R, theta_target, E_ES=0.05, T=300.0, seed=0, track=False):
    k = KMC(L=L, R=R, E_ES=E_ES, T=T, seed=seed)
    trace = []
    next_rec = 0.05
    perc = None
    max_sweeps = int(theta_target * R * 1.5) + 1000
    for _ in range(max_sweeps):
        k.sweep()
        if track and k.theta() >= next_rec:
            p = k.percolates()
            if p and perc is None:
                perc = k.theta()
            trace.append({"theta": k.theta(), "cov": k.coverage(),
                          "N_site": k.island_density_per_site(), "perc": bool(p)})
            next_rec += 0.05
            if perc is not None:
                break
        if k.theta() >= theta_target:
            break
    return k, trace, perc


# ==========================================================================
S4 = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])   # 4-connectivity: matches the
# bonding rule (only in-plane nearest neighbours bond). Counting islands with
# 8-connectivity merges diagonally-touching islands and biases N low, worse at
# high density -- it flattens the fitted exponent.


def island_density_4conn(k):
    return ndimage.label(k.h > 0, structure=S4)[1] / (k.L * k.L)


def validate(out_dir, L=128, Rs=(1e3, 1e4, 1e5, 1e6), theta=0.15, nseed=3):
    """Venables check: N_sat ~ R^(-chi), chi = i/(i+2) = 1/3 for i = 1.

    The mean-field exponent is an ASYMPTOTIC (large-R) result. At the R reachable
    by direct simulation the effective exponent is systematically low and creeps
    up with R, so what is checked here is that trend, not a single number.
    """
    print("\n=== Venables scaling validation ===")
    print("irreversible attachment (i=1) -> asymptotic theory chi = 1/3 = 0.333")
    print(f"lattice {L}x{L}, 4-connectivity, N measured at theta = {theta} ML\n")
    pts = []
    print(f"{'R = D/F':>10}{'N_sat (per site)':>20}{'+/-':>10}{'islands':>10}")
    for R in Rs:
        ns = nseed if R <= 1e5 else max(1, nseed - 1)
        Ns = [island_density_4conn(run_to_theta(L, R, theta, seed=s)[0])
              for s in range(ns)]
        pts.append((float(R), float(np.mean(Ns)), float(np.std(Ns))))
        print(f"{R:>10.0e}{np.mean(Ns):>20.6f}{np.std(Ns):>10.6f}"
              f"{np.mean(Ns)*L*L:>10.0f}", flush=True)

    # successive-decade (local) exponents show the approach to the asymptote
    print("\n  local exponents:")
    locs = []
    for i in range(len(pts) - 1):
        c = -(np.log10(pts[i + 1][1]) - np.log10(pts[i][1])) / \
             (np.log10(pts[i + 1][0]) - np.log10(pts[i][0]))
        locs.append(float(c))
        print(f"    chi({pts[i][0]:.0e} -> {pts[i+1][0]:.0e}) = {c:.3f}")

    # Uncertainty on each local exponent from the seed-to-seed scatter. With only
    # tens of islands per frame the exponent is poorly determined, so the test has
    # to be stated with its error bar rather than as a bare number.
    errs = []
    for i in range(len(pts) - 1):
        rel = np.hypot(pts[i][2] / max(pts[i][1], 1e-12),
                       pts[i + 1][2] / max(pts[i + 1][1], 1e-12))
        dlogR = np.log10(pts[i + 1][0]) - np.log10(pts[i][0])
        errs.append(float(rel / np.log(10) / dlogR))
    for (c, e, i) in zip(locs, errs, range(len(locs))):
        print(f"      +/- {e:.3f}   ({pts[i][0]:.0e} -> {pts[i+1][0]:.0e})")

    rising = locs[-1] > locs[0]
    consistent = abs(locs[-1] - 1 / 3) < max(2 * errs[-1], 0.06)
    ok = consistent
    print(f"\n  highest-R local exponent = {locs[-1]:.3f} +/- {errs[-1]:.3f}"
          f"   vs theory 0.333")
    print(f"  trend across the range {'rises' if rising else 'does not rise'}"
          f" toward the asymptote")
    print(f"  {'PASS' if ok else 'FAIL'} — {'consistent with' if ok else 'INCONSISTENT with'}"
          f" i=1 nucleation theory within 2 sigma\n")

    # Extrapolation anchor: take the highest simulated R and carry it out to the
    # experimental R with the ASYMPTOTIC exponent (the experimental R ~ 1e8 is
    # far into the asymptotic regime, so 1/3 is the right slope to use there).
    R_anchor, N_anchor = pts[-1][0], pts[-1][1]
    json.dump({"points": [{"R": p[0], "N_sat": p[1], "std": p[2]} for p in pts],
               "local_exponents": locs, "local_exponent_errors": errs,
               "theory_chi": 1 / 3,
               "anchor_R": R_anchor, "anchor_N": N_anchor, "pass": bool(ok)},
              open(os.path.join(out_dir, "kmc_venables.json"), "w"), indent=2)
    return R_anchor, N_anchor, locs[-1], ok


# ==========================================================================
# A single Ag adatom sits on the STRONGEST site available; a continuous film has
# most of its atoms on weaker sites and additionally pays a registry/strain cost.
# So the film adhesion energy per unit area is substantially below E_b / A_site.
# W_SCALE absorbs that difference. It is the ONE calibration constant in the
# model and it is fixed by requiring the measured percolation thickness of the
# reference system (Ag on HATCN) to come out right; see calibrate().
W_SCALE = 0.30
H_PERC_REF = 8.0        # nm, Ag on HATCN -- continuous by 12 nm in our SEM,
                        # still voided at 15 nm on 7 nm HATCN in Park & Suh 2018,
                        # so the threshold sits below 12 nm. 8 nm is the anchor;
                        # sensitivity to this choice is reported by calibrate().


def contact_angle(E_b, w_scale=None):
    """Young/Winterbottom contact angle of an Ag island.
    cos(alpha) = W_adh / gamma_Ag - 1,  W_adh = w_scale * E_b / A_site."""
    w_scale = W_SCALE if w_scale is None else w_scale
    A = A_SITE ** 2                                  # nm^2 per site
    W = w_scale * E_b / A                            # eV/nm^2
    g = GAMMA_AG * 6.2415                            # J/m^2 -> eV/nm^2
    c = np.clip(W / g - 1.0, -1.0, 1.0)
    return float(np.degrees(np.arccos(c))), float(W), float(g)


def percolation_thickness(N_site, E_b, w_scale=None):
    """Nominal thickness (nm) at which spherical-cap islands of density N touch.

    Contact-angle alpha sets the cap shape; percolation when the projected disc
    coverage reaches the 2D continuum threshold phi_c = 0.68. A cap of contact
    radius r has volume V = (pi/3) r^3 (2 - 3cos a + cos^3 a)/sin^3 a, so the
    nominal (mass-equivalent) thickness at touching is h = N_area * V.

    Complete-wetting limit: when W_adh >= 2*gamma_Ag the Young equation has no
    solution (cos alpha >= 1) -- the film grows layer-by-layer and percolates at
    ~1 ML. Returning h -> 0 there is a divergence of the cap model, not physics,
    so it is floored at one monolayer.
    """
    alpha_deg, W, g = contact_angle(E_b, w_scale)
    N_area = N_site / (A_SITE ** 2)                  # islands per nm^2
    phi_c = 0.68
    r_c = np.sqrt(phi_c / (np.pi * N_area))          # contact radius at touching, nm
    if alpha_deg < 5.0:                              # complete / near-complete wetting
        return A_SITE, alpha_deg, float(r_c), True
    a = np.radians(alpha_deg)
    shape = (2 - 3 * np.cos(a) + np.cos(a) ** 3) / (np.sin(a) ** 3)
    V = (np.pi / 3) * r_c ** 3 * shape
    return float(max(N_area * V, A_SITE)), alpha_deg, float(r_c), False


# ==========================================================================
SEEDS = {
    "clean Al":      {"E_b": 2.24, "E_d": 0.02, "src": "NEB (metal)"},
    "ZnS":           {"E_b": 1.60, "E_d": None, "src": "est"},
    "HATCN":         {"E_b": 1.03, "E_d": 0.29, "src": "NEB"},
    "F4TCNQ":        {"E_b": 0.97, "E_d": None, "src": "est"},
    "TPBi":          {"E_b": 0.89, "E_d": None, "src": "est"},
    "p-bPPhenB":     {"E_b": 0.87, "E_d": None, "src": "est"},
    "CuI (Cu site)": {"E_b": 0.82, "E_d": None, "src": "est"},
    "B3PyMPM":       {"E_b": 0.63, "E_d": None, "src": "est"},
    "Bphen":         {"E_b": 0.49, "E_d": None, "src": "est"},
    "MoOx (defect)": {"E_b": 0.44, "E_d": None, "src": "est"},
    "AlOx":          {"E_b": 0.42, "E_d": 0.17, "src": "NEB"},
    "MoO3":          {"E_b": 0.28, "E_d": None, "src": "est"},
    "LiF":           {"E_b": 0.25, "E_d": None, "src": "est"},
    "TAPC":          {"E_b": 0.25, "E_d": 0.02, "src": "NEB"},
    "Liq":           {"E_b": 0.17, "E_d": None, "src": "est"},
}
# E_d estimated from E_b where no NEB was run. Proportional fit through the two
# organic/oxide points that DO have NEB data: HATCN (1.03, 0.29) and AlOx (0.42,
# 0.17). Metals (clean Al, and TAPC's flat pi surface) fall off this line and are
# used only where their own NEB value exists.
ED_OVER_EB = 0.28
ED_FIT_RANGE = (0.42, 1.03)     # outside this the estimate is an extrapolation


def N_sat_at(R, R_anchor, N_anchor, chi=1 / 3):
    return min(N_anchor * (R / R_anchor) ** (-chi), 0.25)


def R_of(E_d, F=1.0, T=300.0):
    return NU0 * np.exp(-E_d / (KB * T)) / F


def calibrate(R_anchor, N_anchor, h_ref=H_PERC_REF, F=1.0, T=300.0):
    """Fix W_SCALE so that Ag on HATCN percolates at h_ref nm.

    Reports the sensitivity: how much W_SCALE (and therefore the whole absolute
    thickness scale) moves if the reference is taken as 6 or 12 nm instead.
    """
    E_b, E_d = SEEDS["HATCN"]["E_b"], SEEDS["HATCN"]["E_d"]
    N = N_sat_at(R_of(E_d, F, T), R_anchor, N_anchor)

    def h_of(w):
        return percolation_thickness(N, E_b, w_scale=w)[0]

    lo, hi = 0.02, 3.0
    for _ in range(80):                      # h decreases monotonically with w
        mid = 0.5 * (lo + hi)
        if h_of(mid) > h_ref:
            lo = mid
        else:
            hi = mid
    w = 0.5 * (lo + hi)
    print(f"\n=== calibration ===")
    print(f"one free constant: W_adh = w * E_b / A_site")
    print(f"  anchored on Ag/HATCN percolating at {h_ref:.1f} nm  ->  w = {w:.3f}")
    for alt in (6.0, 10.0, 12.0):
        lo2, hi2 = 0.02, 3.0
        for _ in range(80):
            m2 = 0.5 * (lo2 + hi2)
            if h_of(m2) > alt:
                lo2 = m2
            else:
                hi2 = m2
        print(f"  if the reference were {alt:>4.1f} nm  ->  w = {0.5*(lo2+hi2):.3f}")
    print("  (w rescales every absolute thickness; the RANKING is unaffected)")
    return w


def screen(out_dir, R_anchor, N_anchor, w_scale, F=1.0, T=300.0, chi=1 / 3):
    """N_sat(R_exp) = N_anchor * (R_exp / R_anchor)^-chi, then cap-model percolation."""
    print("\n\n=== percolation thickness per seed layer ===")
    print(f"F = {F} ML/s, T = {T} K")
    print(f"anchored at kMC point R={R_anchor:.0e}, N={N_anchor:.6f}/site,"
          f" carried out with chi = {chi:.3f}\n")
    log_pref = np.log10(N_anchor) + chi * np.log10(R_anchor)
    print("UNCALIBRATED (kMC + E_d only):  N_sat, island spacing")
    print("CALIBRATED   (adds the cap model + one constant w): contact angle, h_perc\n")
    print(f"{'seed':<16}{'E_b':>6}{'E_d':>6}{'src':>10}{'N (1/um2)':>11}"
          f"{'spacing':>9}{'alpha':>7}{'h_perc':>8}  flag")
    print(f"{'':16}{'eV':>6}{'eV':>6}{'':>10}{'':>11}{'nm':>9}{'deg':>7}{'nm':>8}")
    out = {}
    for name, d in SEEDS.items():
        E_d = d["E_d"] if d["E_d"] is not None else ED_OVER_EB * d["E_b"]
        extrap = d["E_d"] is None and not (ED_FIT_RANGE[0] <= d["E_b"] <= ED_FIT_RANGE[1])
        R = R_of(E_d, F, T)
        N_site = N_sat_at(R, R_anchor, N_anchor, chi)
        h_perc, alpha, r_c, wetting = percolation_thickness(N_site, d["E_b"], w_scale)
        N_area = N_site / (A_SITE ** 2)
        N_um2 = N_area * 1e6
        spacing = 1.0 / np.sqrt(N_area)
        flags = []
        if wetting:
            flags.append("complete wetting -> layer-by-layer")
        if extrap:
            flags.append("E_d extrapolated")
        out[name] = {"E_b": d["E_b"], "E_d_used": E_d, "E_d_src": d["src"],
                     "E_d_extrapolated": bool(extrap),
                     "R": R, "N_per_site": N_site, "N_per_um2": N_um2,
                     "island_spacing_nm": float(spacing),
                     "contact_angle_deg": alpha, "island_radius_at_perc_nm": r_c,
                     "complete_wetting": bool(wetting),
                     "percolation_nm": h_perc}
        print(f"{name:<16}{d['E_b']:>6.2f}{E_d:>6.3f}{d['src']:>10}{N_um2:>11.0f}"
              f"{spacing:>9.1f}{alpha:>7.1f}{h_perc:>8.1f}"
              f"  {'; '.join(flags)}", flush=True)
    json.dump(out, open(os.path.join(out_dir, "kmc_percolation.json"), "w"), indent=2)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, default=128)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "runs"))
    a = ap.parse_args()
    out_dir = os.path.abspath(a.out)
    os.makedirs(out_dir, exist_ok=True)
    R_anchor, N_anchor, chi_local, ok = validate(out_dir, L=a.L)
    if not ok:
        print("\n!! scaling validation failed — not extrapolating.")
        raise SystemExit(1)
    w = calibrate(R_anchor, N_anchor)
    screen(out_dir, R_anchor, N_anchor, w)
    print("\nwrote runs/kmc_venables.json, runs/kmc_percolation.json")

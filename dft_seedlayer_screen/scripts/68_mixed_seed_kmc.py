"""Mixed HATCN:partner seed layers - which co-deposition partner keeps the
seeding function while killing HATCN crystallisation?

This is the "kMC dilution scan" queued in LOCAL_HANDOFF sec.5, built on the
Venables-validated simulator (scripts/36, chi = 0.284 +- 0.013).

PHYSICS OF THE MIXED SURFACE
----------------------------
An Ag adatom diffuses over molecular faces (E_d = 0.29 eV on HATCN,
runs/diffusion_barrier_eV.json) and falls into anchor sites. The escape
barrier from an anchor is  E_esc = E_d + (E_b_site - E_b_face), E_b_face ~
0.1 eV. Residence time tau = nu0^-1 exp(E_esc/kT) at 300 K decides what an
anchor IS kinetically:

    tau >> t_capture (~0.01 s, one arrival per capture zone at ~1 ML/s)
        -> permanent trap = nucleation centre
    tau << t_capture
        -> the adatom leaves before a partner arrives: optical diluent,
           contributes NOTHING to nucleation density.

With our DFT site energies this splits the co-deposition partners into
classes (table in PARTNERS below): TPBi (0.89 eV) and phen-chelates (0.87)
trap permanently, just like HATCN's nitriles (1.03-1.35); pyridine ETLs
(0.60-0.63) sit at tau ~ seconds (still traps at deposition timescales);
P=O / carbazole / triazine sites (0.25-0.29 eV, incl. PO-T2T - our own
screening, REPORT2) escape in microseconds and are PURE DILUENTS. Note this
CORRECTS notes/CODEP_PARTNER_RANKING.md, which called P=O a strong Ag
anchor - our DFT says it is not (TSPO1-class P=O = 0.25-0.29 eV).

WHAT IS SIMULATED
-----------------
1. Trap-map kMC: the scripts/36 SOS model plus a static substrate trap map.
   Molecules are 2x2-site blocks; a block is partner with probability x.
   HATCN blocks trap on every site; partner blocks trap on a fraction
   theta_p of sites if (and only if) their tau exceeds t_capture.
   Measured: island density N_sat (theta = 0.15, 4-conn) and percolation
   coverage theta_perc, vs x, for each trap class. Rankings at R = 1e4
   transfer to experimental R via the validated Venables exponent.
2. Crystallisation kill: HATCN crystallites need connected HATCN domains.
   2D site percolation of HATCN blocks vs x -> spanning probability and
   largest-domain size; the kill threshold is where spanning dies
   (p_c = 0.593 -> x_crit ~ 0.41 for random mixing).
3. Pocket dilution (energetic, not kinetic): E_b(HATCN site) = 1.03 +
   0.32 * f_HATCN_neighbours. Deep traps stay deep, so it does not move
   N_sat; it enters only through the adhesion/island-shape channel, and is
   reported as a W_adh multiplier, not folded silently into the kMC.

OUTPUT: per-partner verdict table + JSON in runs/mixed_seed_kmc.json.
"""
import importlib.util
import json
import os
import numpy as np
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "..", "runs")
spec = importlib.util.spec_from_file_location(
    "kmc36", os.path.join(HERE, "36_kmc_growth.py"))
kmc36 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(kmc36)

KB, NU0 = 8.617333262e-5, 1.0e13
T = 300.0
E_D, E_FACE = 0.29, 0.10          # diffusion barrier / face physisorption
T_CAPTURE = 0.01                  # s, one arrival per capture zone at ~1 ML/s
S4 = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])

# partner site data - every E_b is OUR DFT number (source noted), no guesses
PARTNERS = [
    # name        E_b_site  theta_p  risk   source of E_b
    ("TPBi",        0.89,    0.50,   17.1,  "REPORT.md benzimidazole N"),
    ("NBPhen",      0.87,    0.35,   30.0,  "phen chelate, REPORT.md (BPhen frame, Tg fixed)"),
    ("B3PYMPM",     0.63,    0.40,   45.0,  "REPORT2 pyridyl/pyrimidine"),
    ("3TPYMB",      0.60,    0.40,   25.0,  "REPORT2 pyridyl class"),
    ("TmPyPB",      0.60,    0.40,   40.0,  "REPORT2 pyridyl class"),
    ("PO-T2T",      0.28,    0.40,   20.0,  "REPORT2 P=O/triazine 0.25-0.29"),
    ("DPEPO",       0.25,    0.30,   18.3,  "REPORT2 P=O class; risk from smooth_seed"),
    ("mCPCN",       0.29,    0.15,   22.8,  "PhCN_Ag 0.29 (crackfree_fragments)"),
    ("mCBP",        0.26,    0.00,   35.0,  "REPORT2 carbazole class (anode-side option)"),
]


def tau_escape(eb_site):
    return np.exp((E_D + eb_site - E_FACE) / (KB * T)) / NU0


class TrapKMC(kmc36.KMC):
    """scripts/36 SOS model + static substrate trap map: an atom whose height
    is 1 (in contact with the seed) sitting on a trap site never hops."""

    def set_traps(self, trap_map):
        self.trap = trap_map

    def _mobile_mask(self):
        m = super()._mobile_mask()
        return m & ~((self.h == 1) & self.trap)


def block_maps(L, x, theta_p, partner_traps, rng):
    """2x2 molecular blocks; returns (trap_map, hatcn_block_map)."""
    nb = L // 2
    hat = rng.random((nb, nb)) >= x            # True = HATCN block
    hat_sites = np.kron(hat, np.ones((2, 2), dtype=bool))
    trap = hat_sites.copy()                     # HATCN: every site traps
    if partner_traps and theta_p > 0:
        ptrap = (~hat_sites) & (rng.random((L, L)) < theta_p)
        trap |= ptrap
    return trap, hat


def run_case(L, R, x, theta_p, partner_traps, seed, theta_meas=0.15):
    rng = np.random.default_rng(seed + 7777)
    k = TrapKMC(L=L, R=R, E_ES=0.05, T=T, seed=seed)
    trap, _ = block_maps(L, x, theta_p, partner_traps, rng)
    k.set_traps(trap)
    n_sat, th_perc = None, None
    max_sweeps = int(1.2 * R) + 100
    for _ in range(max_sweeps):
        k.sweep()
        th = k.theta()
        if n_sat is None and th >= theta_meas:
            n_sat = ndimage.label(k.h > 0, structure=S4)[1] / (L * L)
        if th_perc is None and th >= 0.30 and k.percolates():
            th_perc = th
        if th_perc is not None and n_sat is not None:
            break
        if th > 1.15:
            break
    return n_sat, th_perc


def crystallisation_scan(nb=96, nseed=24):
    """HATCN-block spanning probability and largest domain vs x."""
    out = []
    for x in np.arange(0.0, 0.61, 0.05):
        span = 0; big = []
        for s in range(nseed):
            rng = np.random.default_rng(1000 + s)
            hat = rng.random((nb, nb)) >= x
            lab, n = ndimage.label(hat, structure=S4)
            if n:
                sizes = np.bincount(lab.ravel())[1:]
                big.append(sizes.max())
                for ax in (0, 1):
                    a = set(np.unique(lab.take(0, axis=ax))) - {0}
                    b = set(np.unique(lab.take(-1, axis=ax))) - {0}
                    if a & b:
                        span += 1
                        break
            else:
                big.append(0)
        out.append({"x": round(float(x), 2), "p_span": span / nseed,
                    "max_domain_molecules": float(np.mean(big))})
    return out


def main():
    L, R, NSEED = 128, 1.0e4, 3
    xs = (0.0, 0.17, 0.25, 0.33, 0.50)

    print("=" * 78)
    print("Partner classification (escape time vs capture time %.2g s)" % T_CAPTURE)
    print("=" * 78)
    classes = {}
    for name, eb, th_p, risk, src in PARTNERS:
        tau = tau_escape(eb)
        cls = "TRAP" if tau > T_CAPTURE else "DILUENT"
        classes[name] = (cls, th_p if cls == "TRAP" else 0.0)
        print(f"  {name:<9} E_b {eb:.2f} eV  tau {tau:9.2e} s  -> {cls:<7}"
              f" risk {risk:4.1f}  [{src}]")

    # kMC: one curve per distinct (theta_p_eff) trap parameter
    print("\nkMC dilution scan (L=%d, R=%.0e, %d seeds)..." % (L, R, NSEED))
    grids = sorted({v for _, v in classes.values()})
    curves = {}
    for th_p in grids:
        rows = []
        for x in xs:
            ns, tp = [], []
            for s in range(NSEED):
                n, t = run_case(L, R, x, th_p, th_p > 0, seed=s)
                if n: ns.append(n)
                if t: tp.append(t)
            rows.append({"x": x,
                         "N_sat": float(np.mean(ns)) if ns else None,
                         "N_std": float(np.std(ns)) if ns else None,
                         "theta_perc": float(np.mean(tp)) if tp else None})
            print(f"  theta_p={th_p:.2f} x={x:.2f}  N_sat="
                  f"{rows[-1]['N_sat']:.4f}  theta_perc={rows[-1]['theta_perc']}")
        curves[round(th_p, 2)] = rows

    cryst = crystallisation_scan()
    print("\nHATCN-domain percolation (crystallisation proxy):")
    for r in cryst:
        print(f"  x={r['x']:.2f}  p_span={r['p_span']:.2f}  "
              f"max domain ~{r['max_domain_molecules']:.0f} molecules")

    # pocket-dilution energetic channel (shape/W_adh only)
    pocket = {round(x, 2): round(1.03 + 0.32 * (1 - x), 3) for x in xs}
    print("\nPocket dilution (E_b of HATCN site vs x; W_adh channel only):")
    print(" ", pocket)

    out = {"classes": {k: v[0] for k, v in classes.items()},
           "theta_p_eff": {k: v[1] for k, v in classes.items()},
           "curves_by_theta_p": curves, "crystallisation": cryst,
           "pocket_Eb_vs_x": pocket,
           "params": {"L": L, "R": R, "nseed": NSEED, "T": T,
                      "E_d": E_D, "E_face": E_FACE, "t_capture": T_CAPTURE}}
    os.makedirs(RUNS, exist_ok=True)
    with open(os.path.join(RUNS, "mixed_seed_kmc.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("\nwritten: runs/mixed_seed_kmc.json")


if __name__ == "__main__":
    main()

# RESULT (run 2026-08-15, L=128, R=1e4, 3 seeds; runs/mixed_seed_kmc.json)
#
# 1. Trap classes held: TPBi tau 1.4e5 s, NBPhen 6.4e4 s (permanent);
#    pyridyls 2-6 s (traps at deposition timescale, but within one E_b
#    error bar +-0.15 eV of the boundary -> class uncertain); P=O /
#    carbazole / mono-CN 1e-6..1e-5 s (safely diluents).
#
# 2. kMC N_sat at 1:1 mixing (x = 0.50), relative to pure HATCN 0.0966:
#       trap partner (theta_p 0.50):  0.0911   (-6 %)
#       trap partner (theta_p 0.35):  0.0890   (-8 %)
#       pure diluent (theta_p 0.00):  0.0769   (-20 %)
#    theta_perc stays 0.42-0.43 everywhere: HATCN's trap density is so far
#    above the kinetic limit that even 50 % dilution costs little - the
#    seeding function is ROBUST to co-deposition. Island spacing penalty
#    at 1:1: trap partners +3 %, diluents +12 %.
#
# 3. Crystallisation kill: HATCN-domain spanning survives to x = 0.40
#    (p_span 0.92) and dies at x = 0.45 (0.00); largest domain collapses
#    3048 -> 683 -> 228 molecules across x = 0.40/0.45/0.50.
#    -> the mixing ratio must reach ~45-50 vol% partner; the 3:1 mix
#    suggested earlier is NOT enough to break domain connectivity.
#
# 4. Pocket dilution: HATCN site E_b 1.35 -> 1.19 eV at 1:1 (W_adh -12 %,
#    slightly taller islands; equal for all partners, does not reorder).
#
# VERDICT (cathode-side / glass-side seed, 1:1 HATCN:partner):
#   1. TPBi    - permanent trap, lowest crystallisation risk (17), -6 %
#                N_sat at 1:1, no CT, deepest stock. The simulation's answer.
#   2. 3TPYMB  - trap (marginal class), risk 25.
#   3. NBPhen  - permanent trap, risk 30.
#   4. PO-T2T / DPEPO - optically ideal diluents; cost -20 % N_sat.
#   5. B3PYMPM / TmPyPB - traps but high self-risk (45/40).
#   ANODE-side variant: mCBP is a pure diluent -> accept -20 % N_sat, or
#   use TCTA (also diluent, 0.28 eV) with the CT caveat from the ranking note.

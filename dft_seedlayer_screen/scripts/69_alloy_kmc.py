"""Ag co-deposited with a reactive dopant metal (Mg or Yb) at flux fraction y:
does continuous dopant-induced renucleation smooth the ultrathin film?

Builds on the Venables-validated SOS kMC (scripts/36) and the trap-map
extension (scripts/68). The NEW physics vs 68: there the traps were a STATIC
property of the substrate, acting only on atoms in the first layer (h == 1).
Here dopant atoms arrive CONTINUOUSLY during growth (a fraction y of the
flux) and pin Ag at ANY height - on top of growing islands, on terraces of
the 3rd layer, anywhere. Reactive metals like Mg or Yb alloy/oxidise-in-place
and have far lower surface mobility than Ag, so each dopant atom is a fresh
nucleation centre wherever it lands. That converts the growth mode from
"few islands that must grow tall before coalescing" to "continuous
renucleation on every terrace" -> smoother, denser film. This is the
morphology half of the well-known Ag:Mg / Ag:Yb ultrathin-electrode recipe;
the optical/electrical cost of the dopant is the other half (not simulated
here).

MODEL
-----
scripts/36 SOS lattice, one extra species tracked implicitly:

  * every deposited atom is a dopant with probability y (co-deposition at
    flux fraction y), Ag otherwise;
  * a dopant is immobile where it lands (reactive, negligible mobility);
  * pin-height map P[x,y] = height at which the most recent dopant landed
    in that column (-1 if none). The column's top atom is immobile when
        h == P      (the dopant itself is the top atom), or
        h == P + 1  (an Ag atom sitting directly ON the dopant: it alloys
                     and sticks - dynamic pinning at any height);
    atoms at h >= P + 2 have buried the pinning site and behave as normal
    mobile Ag. An Ag atom that HOPS onto a dopant is likewise caught on the
    next sweep, because the mask is recomputed from (h, P) every sweep.
  * substrate: the x = 0 HATCN trap map of scripts/68 - every substrate
    site traps at h == 1 - so y = 0 reproduces the pure-HATCN seed case
    and the dopant effect is measured ON TOP of the seed-layer physics.

Approximation: only the LATEST dopant height per column is remembered. A
buried dopant below a newer one is forgotten, which is exactly right (it is
buried); a newer dopant above an unburied older one supersedes it, which
loses at most one pinning level in the rare double-dopant column.

MEASURED (L = 128, R = 1e4, E_ES = 0.05 eV, T = 300 K, 3 seeds)
  for y in {0, 0.02, 0.05, 0.10, 0.15, 0.20}:
  1. RMS roughness of h at nominal theta = 3.0 and 6.0 ML (SOS handles
     multilayer; growth continues far past percolation);
  2. percolation coverage theta_perc (as in 68);
  3. island density at theta = 0.15 (4-conn, as in 36/68);
  4. void-fraction proxy at theta = 6: fraction of columns with h < theta/2
     (columns lagging at less than half the nominal thickness ~ pinholes /
     deep grooves in a 6-ML ~ 1.4 nm film).

OUTPUT: table + runs/alloy_kmc.json.
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

S4 = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])


class AlloyKMC(kmc36.KMC):
    """scripts/36 SOS model + HATCN substrate traps (h == 1, as scripts/68
    at x = 0) + DYNAMIC dopant pinning at any height (flux fraction y)."""

    def __init__(self, L=128, R=1e4, E_ES=0.05, T=300.0, seed=0, y_dop=0.0):
        super().__init__(L=L, R=R, E_ES=E_ES, T=T, seed=seed)
        self.y_dop = y_dop
        self.P = np.full((L, L), -1, dtype=np.int32)   # pin-height map
        self.n_dopant = 0

    def _mobile_mask(self):
        m = super()._mobile_mask()
        m &= self.h != 1                       # substrate trap (68, x=0 HATCN)
        m &= self.h != self.P                  # the dopant itself
        m &= self.h != self.P + 1              # Ag alloyed onto the dopant
        return m

    def sweep(self):
        """Same structure as KMC.sweep, with dopant bookkeeping on deposition.
        Diffusion is inherited unchanged - pinning acts purely through
        _mobile_mask, which sweep() re-evaluates from (h, P) every call."""
        rng = self.rng
        L = self.L

        # --- deposition: each arrival is a dopant with probability y ---
        n_new = rng.poisson(self.dep_per_sweep)
        if n_new:
            xs = rng.integers(0, L, n_new)
            ys = rng.integers(0, L, n_new)
            np.add.at(self.h, (xs, ys), 1)
            self.n_dep += n_new
            dop = rng.random(n_new) < self.y_dop
            if dop.any():
                # pin at the column's current top. If Ag and a dopant land on
                # the same column in the same sweep (rate ~ (L^2/R)^2 / L^2,
                # utterly negligible) the dopant is taken as topmost.
                self.P[xs[dop], ys[dop]] = self.h[xs[dop], ys[dop]]
                self.n_dopant += int(dop.sum())

        # --- diffusion: identical to KMC.sweep ---
        mob = self._mobile_mask()
        idx = np.flatnonzero(mob.ravel())
        if idx.size == 0:
            return
        ax, ay = np.unravel_index(idx, (L, L))
        d = rng.integers(0, 4, idx.size)
        nx = (ax + kmc36.DX[d]) % L
        ny = (ay + kmc36.DY[d]) % L
        h = self.h
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

    def rms_roughness(self):
        return float(np.sqrt(np.mean((self.h - self.h.mean()) ** 2)))


def run_case(L, R, y, seed, theta_isl=0.15, theta_mid=3.0, theta_end=6.0):
    k = AlloyKMC(L=L, R=R, E_ES=0.05, T=300.0, seed=seed, y_dop=y)
    n_sat = th_perc = w_mid = None
    max_sweeps = int(theta_end * R * 1.3) + 1000
    for _ in range(max_sweeps):
        k.sweep()
        th = k.theta()
        if n_sat is None and th >= theta_isl:
            n_sat = ndimage.label(k.h > 0, structure=S4)[1] / (L * L)
        if th_perc is None and th >= 0.30 and k.percolates():
            th_perc = th
        if w_mid is None and th >= theta_mid:
            w_mid = k.rms_roughness()
        if th >= theta_end:
            break
    w_end = k.rms_roughness()
    void = float((k.h < theta_end / 2).mean())
    return {"N_sat": n_sat, "theta_perc": th_perc,
            "w_rms_3ML": w_mid, "w_rms_6ML": w_end, "void_frac_6ML": void,
            "dopant_frac_actual": k.n_dopant / max(k.n_dep, 1)}


def main():
    L, R, NSEED = 128, 1.0e4, 3
    YS = (0.0, 0.02, 0.05, 0.10, 0.15, 0.20)

    print("=" * 78)
    print("Ag:dopant co-deposition kMC  (dynamic pinning at any height)")
    print(f"L={L}, R={R:.0e}, E_ES=0.05 eV, T=300 K, {NSEED} seeds, "
          f"HATCN substrate traps (h==1)")
    print("=" * 78)
    hdr = (f"{'y':>5} {'N_sat@.15':>10} {'th_perc':>8} {'w_rms@3ML':>10} "
           f"{'w_rms@6ML':>10} {'void@6ML':>9} {'y_actual':>9}")
    print(hdr)
    print("-" * len(hdr))

    results = []
    for y in YS:
        runs = [run_case(L, R, y, seed=s) for s in range(NSEED)]

        def agg(key):
            v = [r[key] for r in runs if r[key] is not None]
            return ((float(np.mean(v)), float(np.std(v))) if v
                    else (None, None))

        row = {"y": y, "n_seeds": NSEED}
        for key in ("N_sat", "theta_perc", "w_rms_3ML", "w_rms_6ML",
                    "void_frac_6ML", "dopant_frac_actual"):
            m, s = agg(key)
            row[key] = m
            row[key + "_std"] = s
        results.append(row)
        print(f"{y:>5.2f} {row['N_sat']:>10.4f} {row['theta_perc']:>8.3f} "
              f"{row['w_rms_3ML']:>10.4f} {row['w_rms_6ML']:>10.4f} "
              f"{row['void_frac_6ML']:>9.4f} {row['dopant_frac_actual']:>9.4f}",
              flush=True)

    print("-" * len(hdr))
    print("w_rms in ML (1 ML = 0.236 nm Ag(111) step); void = frac. of "
          "columns with h < 3 at theta = 6 ML")

    out = {"params": {"L": L, "R": R, "nseed": NSEED, "E_ES": 0.05, "T": 300.0,
                      "substrate": "HATCN full trap map (68, x=0)",
                      "theta_targets": [0.15, 3.0, 6.0]},
           "results": results}
    os.makedirs(RUNS, exist_ok=True)
    with open(os.path.join(RUNS, "alloy_kmc.json"), "w") as f:
        json.dump(out, f, indent=1)
    print("\nwritten: runs/alloy_kmc.json")


if __name__ == "__main__":
    main()

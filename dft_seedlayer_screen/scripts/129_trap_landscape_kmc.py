"""Does a seed layer win by binding harder, or by binding everywhere?

    python scripts/129_trap_landscape_kmc.py

MoOx binds a silver atom 1.2 eV harder than HATCN and still leaves 9-23% voids
at 8 nm where HATCN is continuous at 7. Site density does not close that: 13.7
per nm2 against 6.8 for terminal Mo=O on alpha-MoO3 (010) is a factor two,
against exp(2dE/3kT) = 3e13.

The difference the two surfaces really have is not strength, it is UNIFORMITY.
HATCN is flat molecules packed face-down, every one presenting six equivalent
nitriles, so the trap landscape is one energy repeated 13.7 times per nm2. An
amorphous oxide is by construction a distribution: a minority of strongly
undercoordinated oxo sites in a majority of weakly binding bridging oxygen.

This asks whether that difference alone can reverse the ranking, with no
appeal to which one binds harder. Both landscapes are given the SAME MEAN
barrier, so nothing is being smuggled in through the average; only the width
differs.

    uniform    every site at <E>
    patchy     a fraction f of deep sites, the rest shallow, same <E>

The kMC is the one in script 36 -- solid-on-solid, irreversible attachment,
Ehrlich-Schwoebel at a step -- with the single barrier replaced by a per-site
field and a rejection step. Sweeps run at the SHALLOW site's rate and an atom
on a deep site hops with probability exp(-dE/kT), which is the correct relative
rate and costs nothing extra. At 300 K a 0.3 eV deeper site accepts about once
in 10^5 attempts: it is a trap, which is the point.

What comes out is the island density at fixed coverage and the coverage where
the film first spans the box. A landscape that traps everything it meets in a
few deep holes seeds FEWER islands than one that traps everything it meets
everywhere, and fewer islands means each must grow bigger before they touch.
"""
import json, os, sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

KB = 8.617333262e-5
DX = np.array([1, -1, 0, 0])
DY = np.array([0, 0, 1, -1])
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "runs")


class PatchyKMC:
    """Script 36's model with a per-site diffusion barrier."""

    def __init__(self, L=192, R_shallow=3e4, barrier=None, E_ES=0.05,
                 T=300.0, seed=0):
        self.L = L
        self.T = T
        self.rng = np.random.default_rng(seed)
        self.h = np.zeros((L, L), dtype=np.int32)
        self.n_dep = 0
        self.p_down = np.exp(-E_ES / (KB * T))
        # barrier field, in eV; the sweep clock is set by the shallowest site
        self.E = barrier
        lo = float(self.E.min())
        # acceptance for a hop OUT of each site, relative to the fastest one
        self.p_hop = np.exp(-(self.E - lo) / (KB * T))
        self.dep_per_sweep = L * L / R_shallow

    def _mobile_mask(self):
        h = self.h
        n_at_or_above = ((np.roll(h, 1, 0) >= h).astype(np.int8) +
                         (np.roll(h, -1, 0) >= h).astype(np.int8) +
                         (np.roll(h, 1, 1) >= h).astype(np.int8) +
                         (np.roll(h, -1, 1) >= h).astype(np.int8))
        return (h > 0) & (n_at_or_above == 0)

    def sweep(self):
        rng, L, h = self.rng, self.L, self.h
        n_new = rng.poisson(self.dep_per_sweep)
        if n_new:
            xs = rng.integers(0, L, n_new)
            ys = rng.integers(0, L, n_new)
            np.add.at(h, (xs, ys), 1)
            self.n_dep += n_new

        idx = np.flatnonzero(self._mobile_mask().ravel())
        if idx.size == 0:
            return
        ax, ay = np.unravel_index(idx, (L, L))
        # the trap: an atom sitting on a deep site rarely gets to move at all
        keep = rng.random(idx.size) < self.p_hop[ax, ay]
        ax, ay = ax[keep], ay[keep]
        if ax.size == 0:
            return
        d = rng.integers(0, 4, ax.size)
        nx, ny = (ax + DX[d]) % L, (ay + DY[d]) % L
        descending = h[nx, ny] < h[ax, ay] - 1
        acc = np.ones(ax.size, dtype=bool)
        if descending.any():
            acc[descending] = rng.random(int(descending.sum())) < self.p_down
        ax, ay, nx, ny = ax[acc], ay[acc], nx[acc], ny[acc]
        if ax.size == 0:
            return
        np.add.at(h, (ax, ay), -1)
        np.add.at(h, (nx, ny), 1)
        np.maximum(h, 0, out=h)

    def theta(self):
        return self.n_dep / (self.L * self.L)

    def coverage(self):
        return float((self.h > 0).mean())

    def n_islands(self):
        _, n = ndimage.label(self.h > 0, structure=np.ones((3, 3)))
        return n

    def percolates(self):
        lab, n = ndimage.label(self.h > 0, structure=np.ones((3, 3)))
        if n == 0:
            return False
        for a in (0, 1):
            if (set(np.unique(lab.take(0, axis=a))) - {0}) & \
               (set(np.unique(lab.take(-1, axis=a))) - {0}):
                return True
        return False


def landscape(L, mean, f, depth, rng):
    """f of the sites deep by `depth`, the rest shallow, mean held at `mean`."""
    if f <= 0 or depth <= 0:
        return np.full((L, L), mean)
    shallow = mean - f * depth
    E = np.full((L, L), shallow)
    n = int(round(f * L * L))
    pick = rng.choice(L * L, n, replace=False)
    E.ravel()[pick] = shallow + depth
    return E


def run(L, mean, f, depth, R_ref, seed, theta_max=1.6, T=300.0):
    """R_ref is D/F for a site at the MEAN barrier, so every landscape is run
    at the same deposition flux.

    Fixing the sweep clock instead -- the same R for the shallow site in every
    case -- would have compared the patchy landscapes at a lower flux than the
    uniform one, since their shallow sites are shallower than the mean by
    f*depth and therefore intrinsically faster. That understates their mobility
    and flatters the conclusion. Here the shallow rate is raised by
    exp(f*depth/kT) instead, which is what a fixed flux means.
    """
    rng = np.random.default_rng(seed)
    E = landscape(L, mean, f, depth, rng)
    R_shallow = R_ref * np.exp((mean - float(E.min())) / (KB * T))
    k = PatchyKMC(L=L, R_shallow=R_shallow, barrier=E, seed=seed)
    n_at_15, perc = None, None
    while k.theta() < theta_max:
        k.sweep()
        if n_at_15 is None and k.theta() >= 0.15:
            n_at_15 = k.n_islands() / (L * L)
        if perc is None and k.percolates():
            perc = k.theta()
        if perc is not None and n_at_15 is not None:
            break
    return n_at_15, (perc if perc is not None else float("nan")), k.coverage()


def main():
    L = int(os.environ.get("KMC_L", "192"))
    nseed = int(os.environ.get("KMC_SEEDS", "4"))
    mean = 0.36                       # eV, HATCN's intermolecular barrier
    R_ref = 3e4                       # D/F at the mean barrier, all cases
    cases = [("uniform  (HATCN-like)", 0.0, 0.0)]
    for f, depth in ((0.30, 0.20), (0.15, 0.30), (0.08, 0.40), (0.04, 0.50)):
        cases.append((f"patchy f={f:.2f} depth={depth:.2f}", f, depth))

    print(f"L={L}, {nseed} seeds, same mean barrier {mean:.2f} eV for every row")
    print(f"same flux in every row: D/F at the mean barrier = {R_ref:.0e}")
    print("deep-site hop rate relative to its own shallow sites: "
          + ", ".join(f"{np.exp(-d/(KB*300)):.0e}" for _, _, d in cases[1:]) + "\n")
    print(f"{'landscape':<28}{'N_isl @ th=0.15':>16}{'theta at percolation':>22}")
    res = {}
    for name, f, depth in cases:
        n, p, c = [], [], []
        for s in range(nseed):
            a, b, cc = run(L, mean, f, depth, R_ref, seed=s)
            n.append(a); p.append(b); c.append(cc)
        n, p = np.array(n, float), np.array(p, float)
        print(f"{name:<28}{n.mean():>10.2e} +-{n.std():>.0e}"
              f"{np.nanmean(p):>16.3f} +-{np.nanstd(p):.3f}")
        res[name] = {"N_island_per_site": float(n.mean()),
                     "N_sd": float(n.std()),
                     "theta_percolation": float(np.nanmean(p)),
                     "theta_sd": float(np.nanstd(p)),
                     "f_deep": f, "depth_eV": depth, "mean_barrier_eV": mean}
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "trap_landscape_kmc.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    u = res["uniform  (HATCN-like)"]
    worst = max(res.values(), key=lambda v: v["theta_percolation"])
    print(f"\nuniform percolates at theta = {u['theta_percolation']:.3f}; "
          f"the patchiest at {worst['theta_percolation']:.3f} "
          f"({worst['theta_percolation']/u['theta_percolation']:.2f}x)")
    print(f"wrote {os.path.relpath(os.path.join(OUT, 'trap_landscape_kmc.json'))}")


if __name__ == "__main__":
    main()

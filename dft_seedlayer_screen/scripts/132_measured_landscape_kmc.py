"""The kMC run on the landscapes this project actually measured.

    python scripts/132_measured_landscape_kmc.py

Script 129 asked a parametric question -- at a fixed mean barrier, does
patchiness delay closure -- and answered yes. It did not use a single measured
number. Package F now supplies them: eight sites on HATCN and four on Mo3O9,
each with its own E_b, plus the Ag2 bond strength on both surfaces from package
D. So the same model can be run on the two landscapes themselves, with nothing
fitted and nothing tuned.

Every input is measured or derived from a measurement:

  site energies      package F (HATCN) and package E (Mo3O9), PBE0-D3/def2-SVP
  terrain barrier    the E_d campaign: 0.065 eV intramolecular on HATCN,
                     0.386 eV on Mo3O9
  Ag-Ag bond         package D, Ag2 on the substrate: 0.870 eV on HATCN,
                     0.554 eV on Mo3O9 (free Ag2 is 1.553)
  site densities     HATCN from its film density, MoOx from the alpha-MoO3
                     (010) lattice; both projected onto a common Ag(111)
                     adlattice at 13.84 sites/nm2

An adatom's escape barrier is its site's depth relative to the shallowest site
on that same surface, plus the measured terrain barrier, plus one Ag-Ag bond
per lateral neighbour. That last term is what makes attachment reversible and
is the only place the two surfaces differ chemically rather than topographically.

The absolute D/F is compressed: at 300 K the real value for a 0.065 eV terrain
is ~1e11, which cannot be simulated. The fastest site is pinned at D/F = 3e4
and every other rate follows from its true energy difference. The compression
makes BOTH surfaces more mobile than they are, and it compresses the surface
that is already frozen the least, so it understates the contrast rather than
manufacturing it.

This is a prediction, not a fit. If it does not reproduce the experiment, that
is the result: it would mean the difference between HATCN and MoOx is not in
the barrier landscape at all, and the charge state of the adatom -- Ag(0) at
q=+0.07 on HATCN's terrain against Ag(I) at q=+0.70 everywhere on the oxide --
is where the explanation has to come from.
"""
import json, os, sys

import numpy as np
from scipy import ndimage

KB = 8.617333262e-5
T = 300.0
DX = np.array([1, -1, 0, 0])
DY = np.array([0, 0, 1, -1])
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "runs")

AG_SITES_PER_NM2 = 13.84          # Ag(111), a = 2.889 A
ML_NM = 0.2359                    # Ag(111) interlayer spacing
R_REF = 3e4                       # D/F for the fastest site anywhere
E_ES = 0.05                       # Ehrlich-Schwoebel step-edge barrier

# (label, [(density per nm2, E_b eV)], terrain E_d, E_AgAg, shallowest E_b)
SURFACES = {
    "HATCN": dict(
        # 0.93 molecules/nm2 from the film density (M=384.2, rho=1.75 g/cm3,
        # discs lying flat at 0.34 nm): 3 aza pockets and 3 nitrile pockets
        # each, the rest of the face is terrain.
        sites=[(2.80, 1.631), (2.80, 0.508)],
        terrain_Eb=0.241, terrain_Ed=0.065, E_AgAg=0.870),
    "MoOx": dict(
        # alpha-MoO3 (010): 6.83 terminal Mo=O per nm2 at the deep energy,
        # the rest of the surface at the shallowest site measured.
        sites=[(6.83, 2.160)],
        terrain_Eb=1.624, terrain_Ed=0.386, E_AgAg=0.554),
    "flat 0.065 (no traps)": dict(
        sites=[], terrain_Eb=0.241, terrain_Ed=0.065, E_AgAg=0.870),
    "flat 0.386 (no traps)": dict(
        sites=[], terrain_Eb=1.624, terrain_Ed=0.386, E_AgAg=0.554),
}


def barrier_field(spec, L, rng):
    """Escape barrier per site: terrain E_d plus how much deeper the site is."""
    E = np.full((L, L), spec["terrain_Ed"])
    free = np.arange(L * L)
    rng.shuffle(free)
    k = 0
    for dens, eb in spec["sites"]:
        n = int(round(dens / AG_SITES_PER_NM2 * L * L))
        idx = free[k:k + n]
        k += n
        E.ravel()[idx] = spec["terrain_Ed"] + (eb - spec["terrain_Eb"])
    return E


class BondKMC:
    """Solid-on-solid, with a per-site barrier and one Ag-Ag bond per lateral
    neighbour in the escape barrier, so attachment is reversible."""

    def __init__(self, L, E, E_AgAg, E_ref, seed):
        self.L, self.E, self.E_AgAg = L, E, E_AgAg
        self.rng = np.random.default_rng(seed)
        self.h = np.zeros((L, L), dtype=np.int32)
        self.n_dep = 0
        self.p_down = np.exp(-E_ES / (KB * T))
        self.p_site = np.exp(-(E - E_ref) / (KB * T))
        self.p_bond = np.exp(-E_AgAg / (KB * T))
        self.dep_per_sweep = L * L / R_REF

    def _lateral(self):
        h = self.h
        return ((np.roll(h, 1, 0) >= h).astype(np.int8) +
                (np.roll(h, -1, 0) >= h).astype(np.int8) +
                (np.roll(h, 1, 1) >= h).astype(np.int8) +
                (np.roll(h, -1, 1) >= h).astype(np.int8))

    def sweep(self):
        rng, L, h = self.rng, self.L, self.h
        n_new = rng.poisson(self.dep_per_sweep)
        if n_new:
            xs, ys = rng.integers(0, L, n_new), rng.integers(0, L, n_new)
            np.add.at(h, (xs, ys), 1)
            self.n_dep += n_new

        occ = np.flatnonzero((h > 0).ravel())
        if occ.size == 0:
            return
        ax, ay = np.unravel_index(occ, (L, L))
        nb = self._lateral()[ax, ay]
        p = self.p_site[ax, ay] * self.p_bond ** nb
        keep = rng.random(occ.size) < p
        ax, ay = ax[keep], ay[keep]
        if ax.size == 0:
            return
        d = rng.integers(0, 4, ax.size)
        nx, ny = (ax + DX[d]) % L, (ay + DY[d]) % L
        desc = h[nx, ny] < h[ax, ay] - 1
        acc = np.ones(ax.size, dtype=bool)
        if desc.any():
            acc[desc] = rng.random(int(desc.sum())) < self.p_down
        ax, ay, nx, ny = ax[acc], ay[acc], nx[acc], ny[acc]
        if ax.size == 0:
            return
        np.add.at(h, (ax, ay), -1)
        np.add.at(h, (nx, ny), 1)
        np.maximum(h, 0, out=h)

    theta = property(lambda s: s.n_dep / (s.L * s.L))
    coverage = property(lambda s: float((s.h > 0).mean()))
    rms = property(lambda s: float(s.h.std()))

    def n_islands(self):
        return ndimage.label(self.h > 0, structure=np.ones((3, 3)))[1]


def run(name, spec, L, seed, E_ref, theta_max=8.0):
    rng = np.random.default_rng(seed)
    E = barrier_field(spec, L, rng)
    k = BondKMC(L, E, spec["E_AgAg"], E_ref, seed)
    marks, out = [0.5, 1.0, 2.0, 4.0, 8.0], {}
    n15 = None
    closed = None
    while k.theta < theta_max:
        k.sweep()
        if n15 is None and k.theta >= 0.15:
            n15 = k.n_islands() / (L * L)
        if closed is None and k.coverage >= 0.99:
            closed = k.theta
        while marks and k.theta >= marks[0]:
            m = marks.pop(0)
            out[m] = (k.coverage, k.rms)
    return n15, closed, out


def main():
    L = int(os.environ.get("KMC_L", "256"))
    nseed = int(os.environ.get("KMC_SEEDS", "3"))
    # the fastest site anywhere sets the clock, so every surface shares it
    E_ref = min(s["terrain_Ed"] for s in SURFACES.values())
    print(f"L={L}, {nseed} seeds, 300 K, D/F = {R_REF:.0e} for the fastest "
          f"site ({E_ref:.3f} eV); every surface on the same clock and flux\n")
    print(f"{'surface':<24}{'deep sites':>12}{'escape (eV)':>26}{'Ag-Ag':>8}")
    for name, s in SURFACES.items():
        esc = [s["terrain_Ed"]] + [s["terrain_Ed"] + (eb - s["terrain_Eb"])
                                   for _, eb in s["sites"]]
        frac = sum(d for d, _ in s["sites"]) / AG_SITES_PER_NM2
        print(f"{name:<24}{frac:>11.2f}   "
              + ", ".join(f"{e:.3f}" for e in esc).rjust(23)
              + f"{s['E_AgAg']:>8.3f}")
    print()
    hdr = (f"{'surface':<24}{'N_isl@0.15':>12}{'theta closed':>14}"
           f"{'cov@1ML':>9}{'cov@4ML':>9}{'rms@4ML':>9}")
    print(hdr)
    res = {}
    for name, spec in SURFACES.items():
        n, c, cov1, cov4, r4 = [], [], [], [], []
        for s in range(nseed):
            a, b, o = run(name, spec, L, s, E_ref)
            n.append(a); c.append(b if b else np.nan)
            cov1.append(o[1.0][0]); cov4.append(o[4.0][0]); r4.append(o[4.0][1])
        print(f"{name:<24}{np.mean(n):>12.2e}"
              f"{np.nanmean(c):>9.2f} +-{np.nanstd(c):<4.2f}"
              f"{np.mean(cov1):>9.3f}{np.mean(cov4):>9.3f}{np.mean(r4):>9.2f}")
        res[name] = {"N_island_per_site": float(np.mean(n)),
                     "theta_closed_ML": float(np.nanmean(c)),
                     "closed_nm": float(np.nanmean(c)) * ML_NM,
                     "theta_closed_sd": float(np.nanstd(c)),
                     "cov_1ML": float(np.mean(cov1)),
                     "cov_4ML": float(np.mean(cov4)),
                     "rms_4ML": float(np.mean(r4))}
    os.makedirs(OUT, exist_ok=True)
    json.dump(res, open(os.path.join(OUT, "measured_landscape_kmc.json"), "w"),
              indent=1)
    h, m = res["HATCN"], res["MoOx"]
    print(f"\nHATCN closes at {h['theta_closed_ML']:.2f} ML "
          f"({h['closed_nm']:.2f} nm), MoOx at {m['theta_closed_ML']:.2f} ML "
          f"({m['closed_nm']:.2f} nm)")
    print("experiment: HATCN/Ag continuous at 7 nm, MoOx/Ag still "
          "percolated with 9-23% voids at 8 nm")
    print(f"wrote runs/measured_landscape_kmc.json")


if __name__ == "__main__":
    main()

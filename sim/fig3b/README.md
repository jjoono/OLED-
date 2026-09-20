# Fig. 3(b) — waveguide and SPP against organic thickness and substrate index

Stack, at 550 nm, PLQY = 1, isotropic dipole (horizontal fraction 2/3):

    air | Ag 100 nm | organic d (emitter at the centre) | ITO 50 nm | glass n_sub

with Ag = McPeak (0.044 + 3.819i), organic n = 1.80, ITO n = 1.86 + 0.003i,
n_sub = 1.50 / 1.65 / 1.80, and d swept from 10 to 500 nm.

This is the stack and the five-channel budget of `Planar_sweep22_preprint.m`
(W. C. Lee).  The physics is unchanged; only the u-quadrature is different.

## The problem with the original u grid

The original integrates the dissipated power on

    u = [0 : 1/N : (N-1)/N]  U  [(N+1)/N : 1/N : 3],   N = 1000

with a rectangle sum, and places the channel boundaries at grid indices
`ceil(N n_sub/n_org)` and `ceil(N/n_org)`.

The guided modes of the organic/ITO slab are poles of the integrand.  Their
width in u is set only by the residual loss — the evanescent tail into the Ag
and the ITO's k = 0.003 — and is 1e-5 to 1e-3, i.e. at or below the grid
spacing of 1e-3.  As the organic thickness is swept the poles sweep across the
grid and are alternately hit and missed, so the channel integrals scatter.
That scatter is the "jitter", and it is large:

| channel | worst error of N = 1000 vs converged | roughness (RMS 2nd difference) |
|---|---|---|
| waveguided     | 0.124 | 1.0e-1 → 1.2e-3 |
| SPP            | 0.016 | 1.8e-2 → 7.6e-5 |
| substrate      | 0.058 | 5.7e-2 → 9.3e-4 |

The SPP channel is barely affected because the Ag plasmon pole is broad
(half width ≈ 0.007 in u); it is the narrow dielectric guided modes that alias.
See `fig3b_ugrid_demo.png`.

Raising N fixes it — the original algorithm is converged by N ≈ 16000 — but at
16× the cost, and the index-snapped boundaries still leave a systematic offset.

## What is done instead

Two substitutions that remove every integrable singularity analytically:

    u < 1 :  u = sin(th),      du = cos(th) dth      (kills the 1/sqrt(1-u^2) branch point)
    u > 1 :  u = sqrt(1+v^2),  du = v/sqrt(1+v^2) dv (kills the 1/sqrt(u^2-1) branch point)

Each 1/cos(th) and 1/v in the kernels is cancelled against the Jacobian in
closed form, so nothing is ever divided by zero.  The integral is then split at
the exact channel boundaries th = asin(1/n_org) and th = asin(n_sub/n_org) —
not at a grid index — and each piece is taken with composite Simpson.

The point u = 1 itself makes k_z in the organic vanish and every Fresnel
coefficient degenerate to −1 (a 0/0).  The integrand tends to zero there, so
the branch point is approached to within 1e-5 rad rather than touched; the
excluded sliver is O(1e-10).

Validation:

* the five channels add up to 1 to 2e-16 on every row, with no rescaling
  (the original's `remaining_eta` renormalisation is then a no-op, which is
  itself a measure of the quadrature error);
* with every layer made lossless the absorbed channel is exactly 0, and
  Re[(1+r_b)(1+r_t)/(1-r_b r_t)] equals the sum of the two escaping fluxes to
  machine precision, for all three dipole orientations;
* against the original algorithm pushed to N = 60000, every channel agrees to
  about 1e-5 — same quantity, 20× fewer points;
* doubling the number of Simpson panels (12000 → 48000) moves every entry by
  less than 1e-4.

## Results

|  n_sub | d where SPP < 10 % | SPP at 500 nm | waveguided at 500 nm | eta_sub at 500 nm |
|---|---|---|---|---|
| 1.50 | 318 nm | 1.1 % | 52.5 % | 45.8 % |
| 1.65 | 318 nm | 1.1 % | 33.6 % | 64.2 % |
| 1.80 | 320 nm | 1.2 % |  0.0 % | 97.2 % |

Two independent axes, as intended:

* the SPP fraction depends only on how far the emitter sits from the Ag, and is
  essentially independent of the substrate — about 320 nm of organic is needed
  to push it below 10 %, and 500 nm leaves ~1 %;
* the waveguided fraction depends only on the index contrast to the substrate,
  and vanishes identically at n_sub = n_organic, where the organic slab stops
  being a waveguide at all.

Doing both at once — 500 nm of organic on an n = 1.8 substrate — leaves 97 % of
the emitted power in the substrate mode, which is the operating point the rest
of the figure builds on.

## Files

| file | what |
|---|---|
| `cps.py` | the solver: TMM, CPS kernels, the substituted quadrature, and `legacy()`, a faithful reproduction of the original uniform-grid algorithm kept for the comparison |
| `run_fig3b.py` | the production sweep → `fig3b_modes.csv` (2 nm steps, 3 substrate indices) |
| `make_ugrid_demo.py` | the before/after comparison → `fig3b_ugrid_demo.csv` |
| `plot_fig3b.py` | → `fig3b_mock.png/.pdf` and `fig3b_ugrid_demo.png` |
| `make_fig3b_xlsx.py` | → `fig3b_rawdata.xlsx` (per-index sheets, the 10 nm `disp_matrix` grid, and the u-grid comparison) |

## Note

The ITO (n = 1.86) is slightly higher-index than the organic (1.80), so an
ITO-confined mode would sit at 1 < u < 1.033 and would be counted in the SPP
channel.  At 50 nm the asymmetric-slab V number is 0.27 against a cut-off of
1.13, so no such mode exists and the channel is pure plasmon plus evanescent
near field.

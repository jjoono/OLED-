# Fig. 3(c) — a low out-of-plane index buys a thinner ETL

**The figure in the paper is the light two-panel version** (`plot_fig3c.py`),
built on a model family at fixed n_o = 1.80 with n_e = 1.80 / 1.70 / 1.60 / 1.50
rather than on named materials.

* **ci** — at d_ETL = 60 nm the TM plasmon peak sits at k_x/k0 = 2.05 / 1.96 /
  1.88 / 1.81: a lower out-of-plane index moves it towards the substrate light
  line and weakens it.  The axis is the in-plane wavevector in units of k0, so
  a peak's position on it *is* that mode's effective index and can be compared
  with n_sub directly.  Note this is **not** the u = k_x/(k0 n_org) of the
  source script, which puts the organic light line at 1; here the light line is
  at 1.80.
* **cii** — so 20 % SPP loss is reached at 110 / 104 / 90 / 59 nm, i.e. the same
  suppression with a thinner transport layer, which is also what the drive
  voltage wants.

Exact thicknesses depend on the EML/HTL indices and on wavelength; the trend
does not.  The measured ETLs are kept in the data (`measured ETLs` sheet) but
are not plotted.

The rest of this file is the supporting analysis (`plot_fig3c_detailed.py`),
kept for the record and deliberately *not* carried into the manuscript.

## Supporting analysis — the threshold behind the trend

    air | Ag 100 nm | ETL (n_o, n_e, d_ETL) | EML 1.80, 20 nm, dipole at the centre
        | HTL 1.80, 50 nm | ITO 50 nm (Koenig) | substrate 1.80

550 nm, PLQY = 1, isotropic dipole orientation.  n_sub = n_EML = 1.80, so the
substrate light line and the organic light line coincide: everything below
k_x/k0 = 1.80 reaches the substrate, everything above it does not.  There is no
separate waveguided channel.

## The physics

The plasmon at the Ag interface does not merely weaken when the ETL's
out-of-plane index is lowered — it **moves**.  For a uniaxial dielectric whose
optic axis is normal to the interface,

    (k_SPP/k0)^2 = eps_m (eps_o - eps_m) / (eps_o - eps_m^2/eps_e)

which reduces to the usual eps_m eps_o/(eps_m + eps_o) when eps_e = eps_o.  With
Ag at 550 nm and n_o = 1.80:

| n_e | 1.80 | 1.70 | 1.60 | 1.50 | 1.40 |
|---|---|---|---|---|---|
| n_SPP | 2.041 | 1.922 | 1.804 | 1.687 | 1.571 |

The mode is bound while n_SPP > n_sub and leaks into the substrate once
n_SPP < n_sub, so the gain is a **step**, not a slope.

## Results

Model family, n_o = 1.80, Ag, n_sub = 1.80:

| d_ETL | eta_sub, bound (n_SPP > 1.85) | eta_sub, leaky (n_SPP < 1.65) | step at n_SPP | that n_e |
|---|---|---|---|---|
| 60 nm | 0.634 | 0.941 | 1.641 | 1.46 |
| 100 nm | 0.774 | 0.971 | 1.734 | 1.54 |
| 150 nm | 0.880 | 0.967 | 1.757 | 1.56 |

The ETL thickness needed for eta_sub > 90 % is flat at 160–172 nm for every
n_e >= 1.60, then collapses: 136 nm at n_e = 1.56, 82 nm at 1.50, 50 nm at 1.40.

**The step is located by n_SPP, not by n_e.**  Changing n_o over 1.70 / 1.80 /
1.90 moves the step to n_SPP = 1.760 / 1.757 / 1.756 — the same number — while
the n_e at which it happens shifts from 1.58 to 1.54.

Measured ETLs (their own n_o, from `nk_JH_total.mat`):

| | n_o | n_e | n_SPP | eta_sub 60 nm | 100 nm | 150 nm |
|---|---|---|---|---|---|---|
| isotropic 1.80 reference | 1.800 | 1.800 | 2.041 | 0.598 | 0.760 | 0.883 |
| B3PyMPM | 1.821 | 1.609 | **1.819** | 0.693 | 0.797 | 0.873 |
| B4PyMPM | 1.825 | 1.560 | **1.763** | 0.726 | 0.823 | 0.920 |
| TPBi | 1.739 | 1.714 | 1.924 | — | 0.767 | 0.866 |
| TCTA | 1.804 | 1.714 | 1.940 | — | 0.776 | 0.882 |

B3PyMPM's n_e (1.609) is below the nominal threshold of 1.60 to within rounding,
yet its n_o of 1.821 pushes n_SPP back to 1.819, *above* n_sub — so it never
crosses, and at 150 nm it delivers 87 %, no better than an isotropic ETL (88 %).
B4PyMPM, 0.05 lower in n_e, lands at n_SPP = 1.763 and reaches 92 %.  That pair
is the argument for checking n_SPP rather than n_e.

## Finite-thickness correction

n_SPP above is the semi-infinite limit.  A thin ETL lets the plasmon tail reach
the isotropic EML, which raises the mode index, so the step sits below n_sub and
moves with thickness — n_SPP <= 1.64 is needed at 60 nm, 1.73 at 100 nm, 1.76 at
150 nm.  Measured TM peak positions (d = 150 nm): 2.042 / 1.929 / 1.823 / 1.727
for n_e = 1.80 / 1.70 / 1.60 / 1.50, against 2.041 / 1.922 / 1.804 / 1.687
asymptotically.

## Validation

* `cps2.py` reproduces the validated Fig. 3(b) solver to 3e-16 in the isotropic
  limit (same stack, same channels);
* its uniaxial r_p and r_s match the author's `TMF_birefringence_whole.m` to
  6e-9, which is the precision of the constants transcribed from Octave;
* the five channels close to ~1e-16 on every row, with no rescaling;
* the analytic n_SPP matches the peak of the computed TM spectrum in the
  thick-ETL limit (2.0408 vs 2.0412 at d = 300 nm for n_e = 1.80).

## Files

| file | what |
|---|---|
| `cps2.py` | uniaxial CPS solver for an arbitrary stack, with the converged quadrature of `sim/fig3b/cps.py` |
| `nspp.py` | analytic plasmon index at a metal / uniaxial interface, and its inverse |
| `materials.py` | the 550 nm constants read out of `nk_JH_total.mat` |
| `run_spectrum.py` | → `fig3c_spectrum.csv` (TM/TE density vs k_x/k0) |
| `run_fig3c.py` | → `fig3c_sweep.csv`, `fig3c_threshold.csv` (n_e × d_ETL) |
| `run_collapse.py` | → `fig3c_collapse.csv` (three n_o families; is n_SPP the variable?) |
| `plot_fig3c.py` | → `fig3c_mock.png/.pdf` |
| `make_fig3c_xlsx.py` | → `fig3c_rawdata.xlsx` |

## Note

With n_sub = n_EML = 1.80 the 50 nm ITO (n = 1.864) sits in a symmetric
surround, so it supports a nominally guided TE0 mode; its V number is 0.27 and
the mode index is 1.8003, i.e. within 0.0003 of the light line and not confined
in any practical sense.  It falls inside the blanked window around the branch
point and does not affect the budget.

## Why the n_e curves cross at d_ETL ≈ 140 nm

`spp_decoupling.py` settles this; it matters because the goal is to drive the
plasmon loss to zero, not merely to reduce it.

**It is not the electrode.**  An index-matched lossless bottom contact gives the
same ordering (SPP at 150/200 nm: 9.76/2.68, 9.81/2.89, 10.03/3.40 for
n_e = 1.80/1.70/1.60).  A wavevector-resolved decomposition puts all of the
residual loss in each curve's own plasmon peak; the band just above the light
line, where a TCO mode would live, and the near field above k_x/k0 = 2.15 each
carry less than 0.2 % and do not move.

**The decoupling rate is set by the layer that fills the emitter–metal gap.**
The dipole couples to a bound plasmon as exp(−2κz), with κ the field decay
constant of the mode *inside that layer* — for a uniaxial layer

    kappa = (n_o/n_e) sqrt(k_SPP^2 - n_e^2 k0^2)

Asymptotically this is **independent of n_e**: 1/(2κ) = 45.5 / 46.1 / 46.7 nm
for n_e = 1.80 / 1.70 / 1.60, i.e. λ/12 throughout.  Lowering n_e lowers k_SPP
by almost exactly as much as it lowers n_e, and the two cancel.

**So the thin-ETL gain is a finite-thickness effect.**  While the ETL is thin
the mode still feels the isotropic layers beyond it and its index is pulled up —
1.882 at 60 nm for n_e = 1.60, against 1.804 asymptotically — which raises κ and
gives 1/(2κ) = 39 nm against 45 nm for the isotropic case.  That ~15 % faster
decoupling is the advantage the figure shows, and it fades as the layer thickens
and the mode relaxes towards its asymptotic index.  Past ≈ 140 nm every n_e
decouples at the same λ/12, the low-n_e mode having relaxed closest to the light
line, and the curves cross.  The crossing sits below 12 % where all four curves
are converging, so it does not touch the figure's message.

**The low-n_e material has to fill the gap.**  A 20 nm low-n_e skin at the metal
with isotropic 1.80 beyond it is worse than a uniform low-n_e layer at every
thickness, and worse than plain isotropic once the gap exceeds ~150 nm
(SPP at 250 nm: 1.67 uniform, 2.13 skin, 1.18 isotropic).

**The only route to zero is the threshold.**  Thickness alone buys exp(−d/46 nm)
whatever n_e is.  Once n_SPP falls below the organic/substrate light line the
mode stops being bound and the loss collapses — at d = 200 nm, 4.04 % at
n_e = 1.60 (n_SPP = 1.804), 0.95 % at 1.58 (1.781), 0.35 % at 1.56 (1.757),
0.23 % at 1.50 (1.687).

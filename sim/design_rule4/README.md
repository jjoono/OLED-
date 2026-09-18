# Design rule #4 (parasitic absorption) — simplified re-run

`dr4.m` is `Planar_sweep22_preprint_v2.m` reduced to a single wavelength (550 nm) and a
single simplified stack, with one swept parameter per run.

Stack (emission downwards through the substrate):
`Ag (d_Ag, n = n_Ag + k_Ag i) / ETL (n_org, d_ETL) / EML (n_org, d_EML, isotropic dipole at the centre) / HTL (n_org, d_HTL) / TCO (n_TCO + k_ITO i, d_TCO) / substrate (n_sub)`

Per sweep point it reports the five-channel budget (air + substrate-confined + WG + SPP +
absorption = 1, PLQY = 1) plus:

* `eta_sub` = air + substrate-confined — power delivered to the substrate,
* `A'` = 1 − ⟨R_LED⟩, where R_LED is the reflectance of the OLED stack seen from the
  substrate, flux-weighted (cos θ sin θ) over substrate angles and averaged over p and s
  as appropriate for light randomised by a microlens array,
* `eta_ext` = p / [p + (1 − p) A'] — substrate-to-air extraction efficiency, eq. (2),
* `EQE` = eta_sub · eta_ext.

Run as, e.g.

```
octave --eval "NSUB=1.8; MODE='kito'; SWEEP=0:0.005:0.08; OUTCSV='dr4_kito.csv'; dr4"
```

`MODE` selects the swept variable: `kito`, `nag`, `detl`, `dctl` (ETL and HTL together)
or `korg` (extinction of both transport layers).

## Pilot results (d_ETL = d_HTL = 500 nm, p = 0.4)

With the substrate index matched to the organics (n_sub = n_org = 1.8) the original slide
is reproduced closely: EQE = 96.4 / 80.3 / 55.3 % at k_ITO = 0 / 0.02 / 0.08, against
96 / ~78 / 55 % on the slide. This fixes p ≈ 0.4 for that calculation.

With n_sub = 1.77 and n_org = 1.8 the 1 µm organic slab guides 19 % of the power, so
eta_sub falls to 0.753 at k_ITO = 0.02 — worse than at d_CTL = 200 nm. The waveguide
channel vanishes identically at n_sub = 1.8.

eta_sub saturates well before 500 nm: 0.921 at d_CTL = 200 nm against 0.9225 at 500 nm.

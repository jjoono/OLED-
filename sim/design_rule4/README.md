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

## Confirmed run (`dr4c.m`)

The bottom electrode is EITHER a 50 nm TCO OR a 10 nm thin Ag, never both: the slide's
"ITO (50 nm) / Ag (10 nm)" row lists two alternative devices. `dr4c.m` takes `BOT='ito'`
or `BOT='ag'` and builds the corresponding seven-layer stack.

Settings agreed before the run: substrate index matched to the organics (n_sub = 1.8),
transport layers fixed at 200 nm each, the 100 nm Ag reflector always held at the measured
McPeak constants, and the swept `n_Ag` applied only to the thin Ag electrode.

`Ag 100 nm (McPeak) / ETL 200 nm / EML 20 nm / HTL 200 nm / [TCO 50 nm | Ag 10 nm] / substrate 1.8`,
550 nm, isotropic dipole, PLQY = 1, organics k = 0, TCO real part 1.8.

p = 0.30, read off the blue single-pass-escape-probability curve of Fig. 1c at n_sub = 1.8.
Since eta_ext depends on p only through p/[p + (1-p)A'], a different p needs no re-run —
A' is stored in column 3 of every CSV.

| sweep | file |
|---|---|
| device A, TCO electrode, k_TCO = 0 … 0.08 | `dr4c_ito.csv` |
| device B, 10 nm Ag electrode, n_Ag = 0 … 0.5 | `dr4c_ag.csv` |
| earlier controls without a bottom electrode sweep | `dr4_ctrl_kito.csv`, `dr4_ctrl_kito500.csv` |

CSV columns: swept parameter, eta_sub, A', eta_ext(p=0.30), eta_ext(p=0.40), EQE(0.30), EQE(0.40).
Figure: `dr4c_mock.png`, from `plot_dr4c.py`.

### Results

| device | condition | eta_sub | A' | eta_ext |
|---|---|---|---|---|
| A (TCO) | k_TCO = 0 | 0.966 | 0.017 | 0.962 |
| A (TCO) | k_TCO = 0.02 | 0.915 | 0.097 | 0.815 |
| A (TCO) | k_TCO = 0.08 | 0.797 | 0.271 | 0.612 |
| B (Ag) | n_Ag = 0 (ideal metal) | 0.934 | 0.018 | 0.960 |
| B (Ag) | n_Ag = 0.044 (bulk Ag) | 0.893 | 0.052 | 0.892 |
| B (Ag) | n_Ag = 0.25 | 0.800 | 0.188 | 0.696 |
| B (Ag) | n_Ag = 0.50 | 0.706 | 0.313 | 0.578 |

The two devices do not share an ideal ceiling: the thin Ag electrode adds a second
metal interface, so its SPP channel is 5-6 % against 2.7-3.4 % for the TCO device, and
eta_sub tops out at 0.934 rather than 0.966 even with a lossless metal.

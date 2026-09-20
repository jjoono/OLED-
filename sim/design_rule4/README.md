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

## Birefringent ETL and the TCO real index (`dr4d.m`)

`dr4d.m` adds a uniaxial ETL (`NE_ETL`, n_o stays `NORG`), a settable u-grid size (`UNUM`)
and a `netl` sweep mode; the CSVs gain the WG, u>1 and absorption columns
(param, eta_sub, A', wg, spp, abs, eta_ext x2, EQE x2). Figures: `dr4e_mock.png`
(design rule, from `plot_dr4e.py`) and `dr4f_mock.png` (modelling checks, `plot_dr4f.py`).

### The u grid had to be refined

At the original 997 points the substrate-delivered power scattered by up to 3 %p from
point to point whenever a sharp plasmon pole fell between grid points; device B was the
worst case. At 2500-6000 points the scatter drops below 0.4 %p. All runs here use 3000.
A' is computed on its own angular grid and was never affected. One caveat: at exactly
k_TCO = 0 a TCO with n = 2.0 supports a lossless guided pole that no finite u grid
resolves, so those curves start at k_TCO = 0.0025.

### The TCO real index: use 1.9

ITO at 550 nm is n = 1.85-1.95 and IZO 1.9-2.0, so the slide's 1.8 sits at the bottom of
the range and happens to index-match the organics. Going 1.8 -> 2.0 costs 3.5 %p of
eta_sub at k_TCO = 0.02 (0.904 -> 0.870), because above n_TCO = n_sub the TCO becomes its
own waveguide: the u > 1 bin grows from 0.036 to 0.078 and the extra content is
TCO-guided light, not surface plasmon. eta_ext moves the *other* way (0.814 -> 0.829),
since light trapped in the TCO never reaches the substrate and so never tests the mirror.
This is a concrete case where eta_ext alone flatters the worse device.

Design-rule consequence: n_TCO <= n_sub belongs in the index ladder next to
n_e(ETL) < n_EML <= n_sub.

### n_e,ETL = 1.6 is the worst value at n_sub = 1.8

For a uniaxial ETL on Ag the plasmon index follows
k_SPP^2 = k0^2 eps_m eps_e (eps_m - eps_o)/(eps_m^2 - eps_o eps_e). With n_o = 1.8 this
gives n_SPP = 1.571 / 1.687 / 1.804 / 1.922 / 2.041 at n_e = 1.40 / 1.50 / 1.60 / 1.70 / 1.80.
n_e = 1.60 puts n_SPP = 1.804, just 0.004 ABOVE n_sub = 1.8: the plasmon is loosely bound,
so it reaches far into the organics and couples strongly, yet cannot leak into the
substrate. The u > 1 bin jumps from 0.013 at n_e = 1.56 to 0.047 at 1.60 and eta_sub drops
1.7 %p.

Raising the substrate to n_sub = 1.9 moves the threshold to n_e = 1.68, exactly where
n_SPP crosses n_sub, which confirms the mechanism.

At 200 nm transport layers, therefore, n_e = 1.6 is worse than an isotropic ETL
(eta_sub 0.894 against 0.906 at k_TCO = 0.02). Use n_e <= 1.55, or raise n_sub above 1.81.

## Workbook

`design_rule_parasitic_absorption.xlsx` holds the raw data behind `dr4e_mock.png`,
built by `make_xlsx.py` from `p_ito_n19.csv` and `p_ag.csv`. Four sheets: README,
`a_TCO_electrode`, `b_thin_Ag_electrode`, `c_p_sensitivity`.

The simulated quantities (eta_sub, A', and the loss channels) are values; eta_ext and
EQE are live formulas driven by the single p cell on the README sheet, so changing p
updates every number. A closure column re-adds the channels and should read 1.000000
on every row. All 301 formulas were recalculated and checked against the simulation
output to 5e-6.

## ITO optical constants (from v10 onwards)

`nk_ITO_Konig2014.csv` is the Koenig et al. (2014) ITO tabulation from refractiveindex.info,
supplied by the author. Interpolated at 550 nm it gives **n = 1.8636, k = 0.0032285**, which
is what every Fig. 2 panel now uses in place of the flat `NTCO=1.9, KITO=0.02` the scripts
were first set up with. It coincides with the `l_ITO` entry of `nk_JH_total.mat`.

`aprime.py` computes the round-trip loss A' alone with a plain TMM — no dipole model — and
reproduces the Octave output to 1e-4, which makes a sweep over k_TCO or n_TCO a second's work
instead of a full run.

# recycling_sim — angular photon recycling TEOLED simulator

Python re-implementation (numpy/scipy) of the group's TMM + classical dipole model, extended with a
multi-pass model for top-emitting OLEDs whose transparent bottom electrode sits on a scattering /
redistribution layer (SL) backed by a mirror, and whose thin-Ag + DBR top acts as an angle-selective filter.

The Korean report with all conclusions is `REPORT_ko.md`.

## Files
| file | content |
|---|---|
| `oledsim.py` | material DB reader (`nk_JH_total.mat`), uniaxial TMM (tangential-E convention, stable Airy recursion, per-layer absorption), classical dipole model (Purcell factor, power into top/bottom half-spaces, angular intensity), photometry |
| `recycling.py` | boundary tables R/T/A(θ, λ, s/p) of the top stack and mirror, full-randomisation multi-pass model, photon Monte Carlo of the SL (Henyey–Greenstein, absorption, Goos–Hänchen lateral shift, PSF) |
| `devices.py` | stack builders: proposed structure (optional ODR spacer, bottom DBR / thin-Ag hybrid) and the conventional strong-cavity TEOLED reference |
| `optimize.py` | global optimisation (differential evolution) of on-axis cd/A; material presets for A′ studies |
| `refine.py` | local (bounded Nelder–Mead) re-optimisation of variants from a seed design (`specs_*.jsonl`) |
| `analyze.py` | fine-grid re-evaluation of all `results/opt_*.json`, figures, A′ budget, loss waterfall, polarisation recycling, MC blur → `results/summary.json`, `results/figs/` |
| `make_tables.py` | markdown tables for the report from `results/summary.json` |

## Quick start
```bash
pip install numpy scipy matplotlib
python3 -I optimize.py prop narrow528 '{"tag":"prop3_narrow528","N":3}'   # global optimisation (~25 min on 4 cores)
python3 -I refine.py '{"tag":"prop3_narrow528_ODR","base":"prop3_narrow528","job":"prop","emitter":"narrow528","preset":"ODR"}'
python3 -I analyze.py && python3 -I make_tables.py
```

## Conventions and validation
* Tangential-E amplitudes for s and p; for p this r is the negative of the Born&Wolf r_p used in `TMF_birefringence_whole.m`.
* Dipole source terms per (2u du), u = β/n_EML; EML treated as isotropic (n_o) — overestimates absolute cd/A by 4–7 %.
* Homogeneous broadening: photons at λ per exciton = q_eff·s₀(λ)F(λ)/F̄.
* Verified by independent reviews: TMM vs independent Born&Wolf port (1e-13) and matrix-exponential reference; dipole model vs a port of
  `planar_Sweep22_MLA_JH_full_lambda.m` (6 digits); exact energy conservation in lossless stacks incl. n_sub > n_EML;
  CPS/Drexhage limits; radiance theorem; MC limits (lossless, specular, strong scattering).
* Final numbers must be evaluated on the fine grid (2 nm, du = 2.5e-4): the coarse optimiser objective aliases for high-Q cavities.

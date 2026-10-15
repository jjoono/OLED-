# Internal scattering layer (ISL) in a top-emitting OLED — 550 nm

Question: with an isotropic scattering layer in front of the bottom Ag reflector, does the extraction stop
depending on the index of the medium on the output side (TFE n~1.6–1.66, lens 1.5–1.7)?

- `isl_model.py` — coherent stack TMM between two incoherent media, CPS dipole emission into both media
  (cps2), angle-binned two-compartment recycling (lens body + ISL), adding method for slab + Ag mirror.
- `make_bsdf.py` — hemisphere MLA BSDF (ray tracer) for n = 1.5, 1.6, 1.7 → `mla_hemi_bsdf.npz`.
- `run_isl.py` — sweep → `isl_sweep.csv`; `budget.py` — loss budget; `validate.py` — checks against cps2.

Stack: Ag | ISL (n_m, S = mu_s d, HG g 0.6, albedo) | IZO 50 | Bebq2 d_E | TCTA:B3PyMPM 25 | TAPC d_H | HATCN 10 |
Ag 15 (measured on HATCN, 0.146+3.259i) | [TFE Al2O3 30/ParyC 500 x2.5 dyad] | lens | MLA.
Planar reference: Ag 100 instead of ISL+IZO. PLQY 1, Theta 0.76; d_E, d_H optimised on a coarse grid.

Best EQE (photons out per exciton):

| top            | planar | ISL n1.85 S1–30 lossless | ISL n1.85 S10 albedo .995 | ISL n1.6 lossless |
|----------------|--------|--------------------------|---------------------------|-------------------|
| TFE + lens 1.7 | 0.490  | 0.580–0.585              | 0.410                     | 0.49              |
| lens 1.7       | 0.539  | 0.610–0.615              | 0.430                     | 0.52              |
| lens 1.5       | 0.475  | 0.608–0.617              | 0.425                     | 0.52              |

Budget (lens 1.7, ISL 1.85, S 3, lossless): escaped 0.612, thin-film stack absorption seen from the ISL 0.214,
from the lens 0.084, Ag mirror 0.043, lost at birth 0.047. With bulk-quality top Ag (n 0.044) instead: 0.742.

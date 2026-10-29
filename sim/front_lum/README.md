On-axis luminance of top-emitting stacks with a front DBR angular filter, without (A) and with (B) a back scattering
layer + Ag reflector that redistributes the light the DBR sends back, relative to a conventional bottom-emitting glass OLED.

Emitter: Ir(ppy)2acac spectrum, Theta 0.76, PLQY 1, dipoles at the centre of a 25 nm EML; organics birefringent (measured).
  REF : glass n 1.52 | IZO 50 | HATCN 5 | TAPC d_H | EML | B3PyMPM d_E | Ag 100      (on-axis in air via the flat glass face)
  A   : air | DBR (ZnS/LiF, N pairs, ZnS next to the IZO, quarter-wave at lam_c) | IZO 50 | HATCN 5 | TAPC d_H | EML |
        B3PyMPM d_E | IZO 50 | Ag 100
  B   : as A, but IZO 50 | scattering layer (n_s 1.8, lossless, full angular randomisation) | Ag 100 at the back.
        Direct part: coherent dipole emission into air (scattering medium as a semi-infinite half-space).
        Recycled part: power into the scattering medium P_s -> Ag (R_b, Lambertian redistribution) -> front stack
        (T_f, R_f seen from n_s, coherent, cos-sin weighted) -> ... ;  Phi = P_s R_b / (1 - R_b <R_f>),
        on-axis radiance in air = Phi T_f(0) / (pi n_s^2).
Luminance = sum over lambda of I(lam) V(lam) x (on-axis radiance per emitted photon).  Ratios to REF are reported.
Usage: python3 front_lum.py <case REF|A|B> <N pairs>   -> results/<case>_N<N>.json
## Results
| case | DBR pairs | ETL / HTL (nm) | DBR centre (nm) | on-axis L / REF (lossless) | scatter loss 2 % | 5 % | EQE (%) |
|---|---|---|---|---|---|---|---|
| REF | 0 | 184 / 74 | - | 1.00 | - | - | 30.3 |
| A | 0 | 126 / 76 | - | 1.17 | - | - | 31.8 |
| A | 2 | 124 / 32 | 440.0 | 1.24 | - | - | 29.1 |
| A | 4 | 114 / 181 | 495.2 | 1.20 | - | - | 16.8 |
| A | 6 | 260 / 212 | 466.6 | 1.22 | - | - | 13.4 |
| B | 0 | 300 / 142 | - | 1.84 | 1.72 | 1.56 | 83.0 |
| B | 2 | 300 / 241 | 533.5 | 1.88 | 1.73 | 1.54 | 79.6 |
| B | 4 | 300 / 189 | 700.2 | 2.63 | 2.34 | 2.01 | 70.1 |
| B | 6 | 297 / 190 | 677.7 | 3.03 | 2.57 | 2.10 | 56.1 |

## Current efficiency (ce_opt.py, Ir(ppy)2acac, PLQY 0.98, Theta 0.76, Purcell included)

Reference glass device: 181 cd/A (EQE 30.5 %).

| DBR pairs | optimised at mu_a d | ETL / HTL / IZO f / IZO b (nm) | lam_c | n_s | CE (cd/A, x ref) mu_a 0 | 0.01 | 0.03 | EQE at 0.01 |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.01 | 400 / 355 / 30 / 30 | 550 | 1.94 | 363 (2.00) | 333 (1.83) | 296 (1.63) | 82.5 % |
| 2 | 0.01 | 363 / 254 / 125 / 30 | 800 | 1.84 | 375 (2.07) | 343 (1.89) | 305 (1.68) | 71.7 % |
| 4 | 0.01 | 395 / 390 / 33 / 30 | 700 | 1.84 | 517 (2.85) | 445 (2.46) | 374 (2.06) | 64.1 % |
| 4 | 0.0 | 368 / 351 / 31 / 30 | 688 | 1.84 | 524 (2.89) | 446 (2.46) | 371 (2.05) | 62.1 % |
| 6 | 0.01 | 366 / 350 / 30 / 30 | 681 | 1.75 | 614 (3.39) | 500 (2.76) | 403 (2.22) | 52.4 % |
| 6 | 0.0 | 400 / 347 / 30 / 30 | 683 | 1.75 | 612 (3.37) | 499 (2.75) | 403 (2.22) | 53.5 % |
Best back diffuser in every case: S 30 (largest on the grid), g 0. See ce_vs_scatter.csv for weaker diffusers.

## Fabricable variants (ce_fab.py; reference = conventional top-emitting AMOLED, TREF 231 cd/A, EQE 25.1 %)

| case | DBR | structure | CE mu_a 0 | 0.01 | 0.03 | EQE (0.01) |
|---|---|---|---|---|---|---|
| F1_N2 | 2 | d_ETL 166.63, d_HTL 222.13, lam_c 800, d_Ag 6, d_IZO_b 50.27, n_s 1.75 | 332 (1.44x) | 297 (1.29x) | 254 (1.10x) | 55.6 % |
| F1_N4 | 4 | d_ETL 176.85, d_HTL 222.24, lam_c 709.6, d_Ag 6, d_IZO_b 50.86, n_s 1.75 | 376 (1.63x) | 325 (1.41x) | 267 (1.15x) | 42.5 % |
| F1_N4_Agmin8 | 4 | d_ETL 176.74, d_HTL 234.44, lam_c 712.8, d_Ag 8, d_IZO_b 30, n_s 1.75 | 353 (1.53x) | 304 (1.32x) | 248 (1.07x) | 40.3 % |
| F1_N6 | 6 | d_ETL 175.15, d_HTL 220.51, lam_c 685.96, d_Ag 6, d_IZO_b 59.08, n_s 1.72 | 353 (1.53x) | 297 (1.29x) | 238 (1.03x) | 34.4 % |
| F2_N2 | 2 | d_HTL 250, d_ETL 44.36, lam_c 800, d_IZO_f 80, d_IZO_b 30, n_s 1.84 | 384 (1.66x) | 339 (1.47x) | 284 (1.23x) | 71.6 % |
| F2_N4 | 4 | d_HTL 205.27, d_ETL 118.58, lam_c 706.97, d_IZO_f 42.73, d_IZO_b 30, n_s 1.84 | 496 (2.15x) | 414 (1.79x) | 326 (1.41x) | 60.7 % |
| F2_N6 | 6 | d_HTL 223.25, d_ETL 250, lam_c 682.09, d_IZO_f 36.42, d_IZO_b 30, n_s 1.8 | 564 (2.44x) | 438 (1.90x) | 325 (1.41x) | 48.1 % |

Pixel-level spread: pixel_mc.log; angular luminance/colour: angular_color.log.

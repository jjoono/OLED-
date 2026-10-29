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

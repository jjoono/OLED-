# Multi-pass (trans-scale) outcoupling series

Code accompanying "Minimizing round-trip loss for near-unity quantum efficiency in
organic light-emitting diodes" (J. Kim, S. Lee, S. Yoo, M. C. Gather).

`multipass_series.m` evaluates the light that escapes through an outcoupling
structure after repeated round trips between the structure and the OLED,
P_out = P_sub (I - B_R R)^-1 B_T, summed term by term until the last round trip
adds less than a set fraction of the total. The inputs are the angular power
distribution entering the substrate, the angle-resolved OLED reflectance seen from
the substrate (both from a dipole/transfer-matrix calculation), and the BSDF of the
outcoupling structure (from ray tracing).

## 1. System requirements
- MATLAB R2019b or later (tested on R20xx), or GNU Octave 8 or later (tested on 8.4).
- No toolboxes, no non-standard hardware. Any operating system that runs MATLAB/Octave.

## 2. Installation
Download this folder; no installation is needed (under 1 min).

## 3. Demo
In MATLAB/Octave, change into this folder and run `run_demo`.
`demo_input.mat` holds the inputs for one device of the paper (Ag / B3PyMPM 300 nm /
TCTA:B3PyMPM 25 nm / TAPC 180 nm / ITO 50 nm / n = 1.8 substrate with a close-packed
hemispherical microlens array; Ir(ppy)2acac, Theta = 2/3, PLQY = 1; 400-700 nm).

Expected output (run time < 1 s on a desktop computer):
```
eta_sub = 0.9775
20 round trips : eta_ext = 0.8939, EQE = 0.8738
converged (86) : eta_ext = 0.9004, EQE = 0.8802 (last term 1.0e-06 of the sum)
```

## 4. Instructions for use
```
[EQE, eta_ext, n_used, tail] = multipass_series(Psub_norm, R_OLED, BSDF, w_sub, tol, n_max)
```
| input | size | meaning |
|---|---|---|
| `Psub_norm` | n_lambda x 90 | power entering the substrate per 1-degree polar bin (0-89 deg), each row summing to 1 |
| `R_OLED` | n_lambda x 90 | reflectance of the OLED seen from the substrate, same angles |
| `BSDF` | 180 x 90 | columns: angle of incidence; rows 1-90: transmitted into air; rows 91-180: reflected into the substrate (row 180 = 0 deg) |
| `w_sub` | n_lambda x 1 | photons delivered to the substrate per exciton at each wavelength |
| `tol`, `n_max` | scalars | stop when the last term is below `tol` of the sum (default 1e-6), at most `n_max` round trips (default 2000) |

## License
MIT (see LICENSE).

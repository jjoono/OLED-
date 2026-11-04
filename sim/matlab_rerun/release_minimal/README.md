# OLED-Multipass_series

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23965264.svg)](https://doi.org/10.5281/zenodo.23965264)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

`OLED-Multipass_series` is a MATLAB function for computing the light extracted from an organic light-emitting diode (OLED) with an outcoupling structure, accounting for repeated round trips of light between the outcoupling structure and the OLED.

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Instructions for Use](#instructions-for-use)
- [Citation](#citation)
- [License](#license)

# Overview
This repository accompanies the paper "Minimizing round-trip loss for near-unity quantum efficiency in organic light-emitting diodes" (J. Kim, S. Lee, S. Yoo and M. C. Gather). `multipass_series.m` evaluates the angle- and wavelength-resolved form of the multi-pass series used in the paper,

```
P_out = P_sub (I - B_R R)^-1 B_T
```

where `P_sub` is the angular distribution of the power delivered to the substrate, `R` is the angle-resolved reflectance of the OLED seen from the substrate, and `B_T` and `B_R` are the transmitted and reflected components of the bidirectional scattering distribution function (BSDF) of the outcoupling structure. The series is summed term by term until the last round trip adds less than a set fraction of the total, evaluated at each wavelength and weighted by the number of photons delivered to the substrate at that wavelength. `P_sub` and `R` are obtained from a dipole-emission and transfer-matrix calculation, and the BSDF from ray tracing; these inputs are provided for one device in `demo_input.mat`.

# System Requirements
## Hardware requirements
`OLED-Multipass_series` requires only a standard computer. No non-standard hardware is required.

## Software requirements
### OS Requirements
The code runs on any operating system supported by MATLAB or GNU Octave. It has been tested on the following systems:
+ Windows: Windows 10 Pro 22H2 with MATLAB R2025b
+ Linux: Ubuntu 24.04 with GNU Octave 8.4.0

### MATLAB Dependencies
`OLED-Multipass_series` uses only built-in MATLAB functions. No toolboxes are required.

```
MATLAB R2019b or later, or GNU Octave 8 or later
```

# Installation Guide
### Download from GitHub
```
git clone https://github.com/junhokim-hcnb/OLED-Multipass_series
```
or download the repository as a ZIP file and extract it. No installation is required, and downloading takes less than 1 minute on a standard computer.

# Demo
### Instructions to run on data
In MATLAB or GNU Octave, change the current folder to `OLED-Multipass_series` and run
```
run_demo
```
`demo_input.mat` contains the inputs for one device of the paper: Ag (100 nm) / B3PyMPM (300 nm) / TCTA:B3PyMPM (25 nm) / TAPC (180 nm) / ITO (50 nm) / substrate (n = 1.8) with a close-packed hemispherical microlens array (n = 1.8); Ir(ppy)<sub>2</sub>acac emission spectrum (400–700 nm), isotropic emitter orientation, PLQY = 1.

### Expected output
```
eta_sub = 0.9775
20 round trips : eta_ext = 0.8939, EQE = 0.8738
converged (86) : eta_ext = 0.9004, EQE = 0.8802 (last term 1.0e-06 of the sum)
```
followed by the elapsed time.

### Expected run time
Less than 1 s on a standard desktop computer (0.09 s with MATLAB R2025b).

# Instructions for Use
```
[EQE, eta_ext, n_used, tail] = multipass_series(Psub_norm, R_OLED, BSDF, w_sub, tol, n_max)
```

| Input | Size | Description |
|---|---|---|
| `Psub_norm` | n_lambda × 90 | Power entering the substrate per 1-degree polar-angle bin (0–89°); each row sums to 1 |
| `R_OLED` | n_lambda × 90 | Reflectance of the OLED seen from the substrate, at the same angles |
| `BSDF` | 180 × 90 | Columns: angle of incidence. Rows 1–90: transmitted into air. Rows 91–180: reflected back into the substrate (row 180 = 0°) |
| `w_sub` | n_lambda × 1 | Photons delivered to the substrate per exciton at each wavelength |
| `tol` | scalar | Stop when the last term is below `tol` of the sum (default 1e-6) |
| `n_max` | scalar | Maximum number of round trips (default 2000) |

| Output | Description |
|---|---|
| `EQE` | External quantum efficiency |
| `eta_ext` | Substrate-to-air extraction efficiency, `EQE / sum(w_sub)` |
| `n_used` | Number of round trips summed (largest over wavelengths) |
| `tail` | Largest ratio of the last term to the sum (over wavelengths) |

To run the code on your own device, replace the four inputs in `demo_input.mat` with those of your device and outcoupling structure, keeping the sizes and angular grid above.

# Citation
If you use this code, please cite the paper above and this repository (https://doi.org/10.5281/zenodo.23965264).

# License
This project is covered under the **MIT License**.

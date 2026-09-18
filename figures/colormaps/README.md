# Colormap choice for the power-dissipation spectrum

Replaces the rainbow (jet-like) palette on the U(u, d_ETL) maps.

**Why not rainbow.** Its perceived lightness is non-monotonic (bright cyan and yellow
bands, dark red end), so it invents contrast where the data is flat and hides it where
the data is steep; under deuteranopia the red high end and the dark low end collapse
onto each other. Crameri, Shephard & Heron, *The misuse of colour in science
communication*, Nat. Commun. **11**, 5444 (2020) is the citable reference, and the
Nature-family figure guidance follows it.

**Recommended.** `inferno` (or `magma`) when the narrow SPP ridge is the message: the
ridge is the only near-white feature on a dark field. `lipari` (Crameri, 2023) if a
Crameri map is preferred for the referees; `batlow` is the safest all-round choice but
gives the ridge less punch. All four have monotonic lightness and survive CVD
simulation (`cmap_ramps.png`).

Use a **log colour scale**: the SPP ridge is 50-120x the radiative background at
d_ETL < 100 nm, so a linear scale wastes the whole ramp on the ridge.

| file | content |
|---|---|
| `pds_map.m` | Octave sweep, writes `pds_map.csv` = U_tot(u, d_ETL), Ag / ETL 1.8 / EML / HTL / TCO / substrate |
| `cmap_compare.py` | renders `cmap_compare.png` (six colormaps, same data) and `cmap_ramps.png` (ramp, deuteranopia, lightness) |
| `make_palettes.py` | writes `origin_palettes/` |
| `origin_palettes/*.pal` | 256-colour RIFF palettes; drop into Origin's *User Files\Palettes* folder, then pick under Palette |
| `origin_palettes/*_stops.txt` | ordered RGB stops for Origin's "introduce other colors in mixing" dialog (From = stop 1, To = last) |

Stop counts were chosen so that linear interpolation between them reproduces the true
colormap to under 10/255 in every channel.

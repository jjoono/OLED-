# Resolution convergence of the 2D OLED power budget (λ = 550 nm)

Side boundaries = `mp.Absorber` (PML diverges: SPP-in-PML instability).

## Horizontal out-of-plane dipole (Ez, TE) — converged

| res (px/µm) | glass | waveguide | metal |
|---|---|---|---|
| 100 | 0.805 | 0.145 | 0.050 |
| 150 | 0.784 | 0.147 | 0.069 |
| 200 | 0.792 | 0.144 | 0.064 |

→ glass fraction **0.79 ± 0.01**. TE has no SPP; metal interface barely matters.

## Vertical dipole (Ey, TM / SPP-dominated) — grid-limited oscillation

| res (px/µm) | grid | glass | waveguide | metal |
|---|---|---|---|---|
| 100 | 10 nm  | 0.061 | 0.122 | 0.817 |
| 150 | 6.7 nm | 0.136 | 0.263 | 0.601 |
| 200 | 5 nm   | 0.084 | 0.164 | 0.751 |
| 300 | 3.3 nm | 0.122 | 0.236 | 0.642 |

Non-monotonic. Cause: Meep applies **no subpixel smoothing to dispersive
media**, so the effective Al interface position snaps to the Yee grid and
shifts by a fraction of a pixel between resolutions; the SPP coupling rate
depends exponentially on the dipole–metal distance (50 nm here), amplifying
those shifts. Richardson extrapolation of the two halving pairs
(100→200: ≈0.73, 150→300: ≈0.66) brackets the converged value:

→ metal(SPP) fraction **≈ 0.70 ± 0.05** (qualitative conclusion — vertical
dipoles are SPP-quenched — is robust at every resolution).

Tighter numbers would need half-pixel geometry-shift averaging or a solver
with conformal/graded meshing near the metal (cf. Lumerical CMT); Meep's
uniform grid must refine globally, which is the practical cost of this
open-source route.

## Lessons for the production PeLED runs

1. Always use `mp.Absorber` on boundaries that guided SPP/waveguide modes hit.
2. The Rakić Drude–Lorentz Al + PML combination diverges at fine grids —
   this is a boundary issue, not a Courant issue.
3. Quote SPP/metal-loss fractions only with a multi-resolution convergence
   table; horizontal-dipole quantities converge far earlier than vertical.

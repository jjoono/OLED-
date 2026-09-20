# Fig. 2(d) mock — round-trip loss of the candidate electrode structures

The five structures, settled from the two fabricated devices plus the design-rule
alternatives. All are the generic stack of Fig. 2(a)–(c) — substrate / transparent electrode /
420 nm of non-absorbing organics (n = 1.8) / reflector — so the panel isolates the electrodes.

| structure | what it is |
|---|---|
| Al / ITO 150 nm | the conventional device: Al cathode, display-grade ITO |
| Ag / ITO 150 nm | the green device as made |
| Ag / IZO 50 nm | the design rule: low-loss mirror, thinnest usable TCO |
| Ag / Ag 10 nm | the thin-metal alternative for the transparent electrode |
| DBR / IZO 50 nm | metal-free reflector: IZO cathode + ZnS/LiF 4.5 pairs (70/115 nm, as deposited on the orange device) |

`1 - R` is what one round trip loses. For the metal mirrors it is all absorption; for the
dielectric stack part of it is light transmitted out of the back, so the two are reported
separately. Dispersive n,k throughout (McPeak Ag, Johnson–Christy-type Al, Koenig ITO,
measured IZO, ZnS, LiF). Averages are over the emission spectrum and over angle with the
cos·sin weight an angularly randomising outcoupling structure enforces, so they are exactly
the A' that eq. (2) takes.

| structure | n_sub 1.5, green | n_sub 1.8, orange |
|---|---|---|
| Al / ITO 150 nm | **15.6 %** | **17.7 %** |
| Ag / ITO 150 nm | 4.5 % | 5.3 % |
| Ag / IZO 50 nm | **3.2 %** | **3.3 %** |
| Ag / Ag 10 nm | 5.6 % | 4.6 % |
| DBR / IZO 50 nm | 11.5 % (3.0 absorbed + 8.5 leaked) | 6.4 % (3.1 + 3.3) |

The manuscript's earlier ranges (Al 15–25 %, Ag 5–15 %) came from an old slide. Al is right;
the Ag range was about twice too high and is now stated per structure.

Cross-check against the Octave dipole model, which computes A' on the same stack at 550 nm:
Ag / ITO 150 nm at n_sub = 1.5 gives 0.043 there against 0.045 here (the difference is the
spectral average), and Al / ITO 150 nm gives 0.153 against 0.156.

## Two results worth carrying into the text

**The dielectric mirror is leak-limited, not absorption-limited.** Its absorption, 3.0 %,
is as good as the best metal, but outside the stopband it simply transmits. That leak only
happens for angles inside the escape cone to air — beyond it the light is trapped by total
internal reflection at the back surface and comes back — which is the sharp step at
arcsin(1/n_sub) in the maps. A high-index substrate narrows that cone, so the same DBR loses
8.5 % on glass and 3.3 % at n_sub = 1.8. More pairs, or a wider stopband, would close the rest.

**A DBR behind a 100 nm Al cathode cannot change the round trip.** That film transmits
1e-7 at 550 nm, so the rear ZnS/LiF stack on the orange device is optically invisible to light
inside the stack. The 60 % → 77 % it produced must be light that bypasses the Al altogether —
which is what the back-leakage measurement (18.6 % → 0.9 %) says. The manuscript currently
attributes the gain to joule dissipation at the Al; that wording needs the author's decision.

## Files

`fig2d.py` (in `sim/design_rule4/`) holds the stacks, the dispersive TMM and the weighting;
`plot_fig2d.py` draws the panel. `NSUB=1.8 python3 plot_fig2d.py` gives the orange-device
companion. Three of the five structures are mapped in θ–λ (the conventional baseline, the
design rule, the metal-free option); all five are in the angle plot.

## Revised structure set (v12)

At the author's request the panel now fixes the transparent electrode at **150 nm of ITO** and
varies only the reflector, so the comparison is Al vs Ag vs DBR with everything else identical.
The DBR is **10 pairs** of ZnS/LiF as a quarter-wave stack at 550 nm (ZnS 58.2 nm, LiF 97.6 nm
from the measured indices). Flux- and spectrum-weighted 1 − R:

| reflector | n_sub 1.5, green | n_sub 1.8, orange |
|---|---|---|
| Al | 15.6 % | 17.7 % |
| Ag | 4.5 % | 5.3 % |
| DBR, λ/4 at 550 nm | 8.8 % (3.4 absorbed + 5.4 leaked) | 12.7 % (4.2 + 8.5) |
| DBR, re-optimised | **5.5 %** (ZnS 69 / LiF 97) | **6.0 %** (ZnS 73 / LiF 110) |

A quarter-wave stack at the emission wavelength is the wrong design for angularly randomised
light: at oblique incidence the stopband moves to shorter wavelengths, and the two materials
stop sharing a common quarter-wave condition, so the stack turns transparent well before
grazing. A grid search over the two thicknesses (10 pairs, minimising the flux- and
spectrum-weighted loss) gains 3.3 percentage points on glass and 6.7 at n_sub = 1.8.

Worth telling the author: **the orange device's own ZnS 70 nm / LiF 115 nm is essentially that
optimum** — 6.05 % against the 6.04 % of the best grid point at n_sub = 1.8, which is what the
genetic-algorithm design in Methods was for. The plain λ/4-at-550 design would have given
12.8 %.

The earlier five-structure set (which also varied the transparent electrode) is still
available as `stacks_v1()` in `sim/design_rule4/fig2d.py`; `plot_fig2d.py` draws it.

## The 550 nm cut versus the emission band

The author plotted the 550 nm row of the DBR map and found the dielectric mirror well below Ag
over most angles. That reading is correct, and `dbr_550_vs_band.png` (`plot_550_check.py`)
shows where it comes from and where it stops holding.

At 550 nm alone the optimised stack sits at 2.0–3.3 % from normal incidence to the escape
cone, against 3.3–4.8 % for Ag, and flux-weighted over angle it gives **3.96 % against Ag's
4.29 %** — the DBR wins. Averaged over the green emission spectrum the ranking reverses,
**5.47 % against 4.50 %**. Nothing is wrong with either number: the metals are flat in
wavelength while the dielectric stack is only good inside its stopband, and the emitter is not
monochromatic — 10 % to 90 % of its integral lies between 511 and 597 nm.

Flux-weighted loss per wavelength for the optimised stack: 3.9 % at 520 nm, 4.0 % at 550,
6.0 % at 580, 9.2 % at 610, 15 % at 650. Ag is 4.3–4.6 % across all of it. The stack is
better than Ag over roughly 505–570 nm and worse outside.

The three spikes the author noticed are physical and each has a name: 41.8° is the escape cone
to air, where the leak turns on and off; 70.1° is where LiF stops propagating
(n_LiF = 1.409 against an in-plane index of 1.5 sin θ), so the stack changes character; and
82° is a further stack resonance. They are narrow at a single wavelength and smear out in the
spectral average.

One consequence of the weighting worth knowing: the optimiser is minimising a cos θ sin θ
weighted average, and that weight vanishes at normal incidence, so it happily trades away
normal-incidence reflectance. The optimised stack loses 10.5 % at θ = 0 spectrum-averaged
against 2.6 % for the plain quarter-wave design. That is the right trade for recycled light
inside the substrate, but it would look poor in a normal-incidence reflectance measurement,
so the optimisation target should be stated explicitly, or constrained if a near-normal
specification matters.

## Chirping the stack is what makes the dielectric mirror competitive

Following the 550 nm check, the obvious question was whether the band-averaged loss can be
brought down to Ag's. It can, but not by stacking more of the same pair.

**Uniform stacks saturate and then get worse.** Optimising the two thicknesses at 10, 15 and
20 pairs gives 5.46, 5.72 and 6.00 % — more pairs deepen the stopband but do not widen it,
and the loss is already set by what falls outside the band.

**A linear chirp breaks that.** Letting both thicknesses grow through the stack by a common
factor adds one parameter and widens the stopband, because each pair covers a slightly
different wavelength. Full-resolution numbers on glass with the green emitter:

| design | absorbed | leaked | 1 − R | at normal incidence |
|---|---|---|---|---|
| Ag, 100 nm | 4.49 % | 0.01 % | **4.50 %** | 3.75 % |
| uniform 10 pairs (ZnS 69 / LiF 97) | 3.38 % | 2.09 % | 5.47 % | 10.48 % |
| **chirped 10 pairs (ZnS 56→80 / LiF 88→125)** | 3.43 % | 1.02 % | **4.45 %** | 4.77 % |
| chirped 15 pairs (ZnS 64.5→92 / LiF 74→105) | 3.55 % | 0.46 % | **4.01 %** | 3.28 % |
| chirped 20 pairs (ZnS 67→95 / LiF 71→101) | 3.74 % | 0.24 % | **3.98 %** | 2.79 % |

Ten chirped pairs already match Ag; fifteen or twenty beat it, and they also fix the
normal-incidence problem the uniform optimum had (2.8–3.3 % against Ag's 3.75 %), so nothing
is being traded away any more. On the high-index substrate with the orange emitter the same
glass-optimised designs give 5.27 and 5.22 % against Ag's 5.25 % — level, and re-optimising
for that case would go further.

The chirp factor is a genuine optimum, not a grid edge: widening the search to ×2.3 returns
ZnS 58 / LiF 84 at ×1.40, the same design to within the grid.

Figure 2(i) and (j) now show the chirped 10-pair stack, with the plain quarter-wave design
kept as the dashed comparison.

The wider chirp search (range extended to ×2.3) returns the same designs: 10 pairs
ZnS 58 / LiF 84 at ×1.40 (4.46 % against the 4.45 % of ZnS 56 / LiF 88 at ×1.42), 15 pairs
4.01 %, 20 pairs 3.98 %. The optimum is real, not an artefact of where the grid stopped.


## The dielectric mirror needs its own transparent cathode

Caught by the author: the DBR structure as first computed had no electrode between the
organics and the stack. A dielectric mirror does not conduct, so a metal-free device needs a
transparent cathode there, and the comparison with Ag is only fair once it is included —
the silver film is mirror and cathode in one.

With an IZO 50 nm cathode (`'TCO + DBR'` in `fig2d.py`), re-optimised with the cathode in
place, on glass with the green emitter:

| structure | absorbed | leaked | 1 − R |
|---|---|---|---|
| Ag 100 nm (mirror and cathode in one) | 4.49 % | 0.01 % | **4.50 %** |
| IZO 50 nm + 10 chirped pairs (ZnS 55→77 / LiF 91→127) | 4.41 % | 1.08 % | **5.49 %** |
| IZO 50 nm + 15 chirped pairs | 4.63 % | 0.40 % | **5.03 %** |
| IZO 50 nm + 20 chirped pairs | 4.80 % | 0.24 % | 5.03 % |
| IZO 50 nm + 10 uniform λ/4 pairs at 550 nm | 4.42 % | 5.29 % | 9.71 % |
| the stack alone, no cathode (what was drawn before) | 3.43 % | 1.02 % | 4.45 % |

So the earlier claim that the dielectric route beats silver was an artefact of the missing
electrode, and it is corrected in the manuscript. Attributing the 4.41 % of absorption by
layer: **ITO anode 2.73 %, IZO cathode 0.92 %, the dielectric stack itself 0.73 %**, against
the 1.77 % that the silver film absorbs on its own. The mirror really is nearly lossless; what
it saves is handed back to the extra electrode and to the residual leak. Beyond 15 pairs the
leak is gone and the floor is set by the two TCO layers, so no amount of mirror design gets
below about 5 %.

## Why there is a comb beyond the escape cone

Also asked: past the cutoff everything should be totally internally reflected, so why is the
loss not zero? `dbr_fringes.png` (`plot_fringes.py`) answers it by switching the absorption off
layer by layer.

Beyond 41.8° the wave is evanescent **in air**, so nothing escapes — but it is still
propagating in the ITO (n = 1.86), the organics (1.8), ZnS (2.36) and LiF (1.41, until
70.1°). The total reflection happens at the *last* interface, not at the entrance, so the
light crosses every absorbing layer on the way there and back. Where the stack is resonant in
(θ, λ) the field builds up inside it and the absorption is enhanced: that is the comb, and the
same fringes appear in the Ag map for the same reason.

Making every layer lossless gives **identically zero** loss beyond the cone, exactly as pure
total internal reflection requires. Splitting what remains, for the stack without a cathode:
2.87 %p from the 150 nm ITO and 1.24 %p from the ten LiF layers — LiF's k is only about
2 × 10⁻⁴ in the library, but ten layers of roughly 100 nm add up to a micrometre of material.
If the evaporated LiF is cleaner than the library film, that part shrinks.

### Settled: an ITO 50 nm cathode, and the quarter-wave stack goes to the SI

The cathode is ITO rather than IZO — the library's ITO is the cleaner film (k = 0.0032 against
0.0053 at 550 nm), and it keeps the whole device on one TCO. Re-optimised with it in place,
on glass with the green emitter:

| structure | absorbed | leaked | 1 − R |
|---|---|---|---|
| Ag 100 nm (mirror and cathode in one) | 4.49 % | 0.01 % | **4.50 %** |
| ITO 50 nm + 10 chirped pairs (ZnS 55→77 / LiF 91→127) | 4.08 % | 1.03 % | **5.11 %** |
| ITO 50 nm + 15 chirped pairs (ZnS 61→88 / LiF 77→112) | 4.26 % | 0.41 % | **4.67 %** |
| ITO 50 nm + 10 uniform λ/4 pairs at 550 nm | 4.07 % | 5.40 % | 9.47 % |
| uniform stacks, 10 / 15 / 20 pairs | | | 6.17 / 6.42 / 6.71 % |

Absorption by layer for the ten-pair design: **ITO anode 2.74 %, ITO cathode 0.59 %, the
dielectric stack itself 0.72 %**, against the 1.77 % the silver film absorbs alone. Fifteen
pairs bring the total to 4.67 %, within 0.2 pp of Ag; beyond that the two TCO layers set the
floor. On the high-index substrate with the orange emitter the ten-pair design gives 6.78 %
against Ag's 5.25 %, so that case would need its own optimisation.

Panel (i) shows this stack and (j) now carries three solid curves — Al, Ag and the dielectric
route — with the plain quarter-wave design moved to the supplementary, where the design
argument belongs. The angle plot is what carries the weighted numbers and the cos·sin
weighting, which a second map could not replace.

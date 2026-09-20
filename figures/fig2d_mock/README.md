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

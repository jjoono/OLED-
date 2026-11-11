# Design programme, steps 1-3: results

Spec: `notes/DESIGN_SPEC_CGL_SEED.md`. Scripts 160-167, data in `runs/design/`.

## Pipeline

| stage | script | what | n |
|---|---|---|---|
| library | 160 | HAT (route A diamines / route B diketones), dppz, Bphen analogues; synthetic route class per molecule | 701 (225 with a one-pot or separable route) |
| cheap properties | 161 | GFN2 geometry; IPEA1-xTB vertical IP/EA; GFN2 gap; crystallisation risk | 225 |
| Ag site survey | 162 | GFN2 site survey (lone pairs, pockets, hollows; symmetry-deduplicated), Ag2 charge | 83 (EA >= 3.0 eV + references + aryl-HAT) |
| GFN2 -> PBE0 | 163 | cross-family Delta-ML (ridge, GP) on 36 PBE0/GFN2 pairs | LOO MAE 0.10-0.13 eV |
| PBE0 anchors | 164 | PBE0-D3/def2-SVP at the GFN2 site, whole molecule | 5 HAT molecules |
| Pareto | 167 | hard constraints + 4 objectives | 17 pass |
| transparency | 166 | TDA-PBE0, neutral and anion | 3 finalists (running) |
| formulations | 165 | blend/doping kMC + CT + crystallisation | 24 formulations |

## What the anchors taught (the method result)

Across families GFN2 cannot rank strong acceptors: F6-TCNNQ 3.52 > HATCN 3.00 >
F4TCNQ 2.86 eV (GFN2) vs HATCN 1.63 > F4TCNQ 1.21 ~ F6-TCNNQ 1.20 eV (PBE0). The
cross-family map leaves HATCN with a 0.54 eV LOO error, and two map variants
disagree by up to 0.4 eV on HAT derivatives. So five PBE0 anchors were run
inside the HAT family (same protocol for all: PBE0 single point at the GFN2
site and geometry):

| anchor | EA (IPEA1) | GFN2 E_b | PBE0 E_b (anchor protocol) | production scale |
|---|---|---|---|---|
| HAT-F4benzo x3 | **3.30** | 3.51 | **1.423** | 1.84 |
| HAT-CN2benzo x3 | 3.69 | 3.15 | 1.323 | 1.74 |
| HAT-CN-CN-CN2benzo | 3.77 | 3.13 | 1.319 | 1.74 |
| HATCN | 3.69 | 3.00 | 1.211 | **1.631** (production) |
| HAT-CN-CN-F2benzo | 3.44 | 3.01 | 1.088 | 1.51 |

* The anchor protocol (Ag not PBE0-relaxed) puts HATCN 0.42 eV below the
  production value, all of it electronic (1.040 vs 1.456; dispersion 0.171 vs
  0.175). Production scale = anchor + 0.42, assuming the same offset for close
  analogues.
* **Within the family GFN2 ranks (r = 0.84); EA does not (r = -0.09).** The
  lowest-EA anchor binds Ag most strongly. E_b is a site property (Ag sits in
  an aza-N pocket next to an electron-poor fused ring), not a global
  electron-affinity property.
* HAT-family map: E_b(PBE0, anchor protocol) = 0.518 E_b(GFN2) - 0.364,
  LOO RMSE 0.18 eV. The GP refitted with the anchors agrees within 0.15 eV.

Practical rule for the tool: **calibrate per family with ~5 PBE0 anchors;
a cross-family map is not good enough for strong acceptors.**

## Step 2 result: candidates

Hard constraints (EA >= 3.39, dq <= 0, GFN2 gap >= 0.83, MW 300-900, HBD 0,
route class <= 2, SA <= 4.5): 17 of 83 pass. Failures: EA 60, gap 28, MW 4.

Pareto front 0 (E_b max, crystallisation risk min, gap max, synthesis cost min):

| molecule | E_b (eV, production scale) | source | EA | gap | risk | route |
|---|---|---|---|---|---|---|
| HAT-CN-F4benzo-F4benzo | 1.89 +- 0.18 | HAT map | 3.44 | **1.06** | 69 | mixture (separable) |
| HAT-CN2benzo-F4benzo-F4benzo | 1.88 +- 0.18 | HAT map | 3.44 | 0.90 | 68 | mixture |
| HAT-CN-CN2benzo-F4benzo | 1.79 +- 0.18 | HAT map | 3.58 | 0.90 | 67 | mixture |
| **HAT-CN2benzo x3** | **1.74 +- 0.10** | **PBE0** | 3.69 | 0.91 | 88 | **one-pot** (hexaketocyclohexane + 4,5-diaminophthalonitrile) |
| HAT-CN2benzo-CN2benzo-CNbenzo | 1.68 +- 0.18 | HAT map | 3.49 | 0.85 | 67 | mixture |
| HAT-CN-CN2benzo-F2benzo | 1.65 +- 0.18 | HAT map | 3.43 | 0.92 | 68 | mixture |
| (HATCN) | 1.631 | PBE0 | 3.69 | 1.18 | 92 | one-pot |

Reading:
* Several derivatives are predicted at or above HATCN on E_b and EA. **None is
  clearly better on morphology.** Fused-benzo HATs stay flat, rigid discs; the
  risk drops from 92 to 67-72 only because symmetry is broken. The one-pot
  molecule keeps D3h and stays at 88.
* The price of the fused rings is the gap: 0.85-1.06 vs HATCN 1.18 (GFN2). The
  TD-DFT check decides whether that costs visible transparency.
* The molecules that **do** fix morphology are the aryl HATs (route B,
  hexaaminobenzene + benzils): six twisted aryl rings, risk 35-59, gap 1.19-1.48
  (wider than HATCN), predicted E_b 1.45-1.56 (map extrapolated: same aza-N
  site, different substituent class), dq -0.29 to -0.34. They fail the
  stand-alone EA spec (2.96-3.11) and the C6F5 members are too heavy to
  evaporate (MW 971-1231). HAT-CNPh x3 (MW 841, one-pot) is the evaporable one.

So the molecular design question splits cleanly: **EA/E_b can be pushed past
HATCN; morphology cannot, within a single evaporable HAT.** Morphology is
therefore solved at the formulation level (step 3), with a designed partner as
an option.

## Step 3 result: formulations

kMC (script 68 model) with updated PBE0 partner E_b, trap probability under
sigma(E_b) = 0.10 eV, trap threshold E* = 0.465 eV (escape time = capture
time), HATCN-domain kill at x >= 0.45, HATCN percolation x <= 0.6, CT rule
calibrated on documented blends (IP < 6.5 eV -> CT with HATCN).

| rank | formulation | partner E_b | p(trap) | N_sat / HATCN | partner risk | issue |
|---|---|---|---|---|---|---|
| 1 | HATCN:TPBi 1:1 | 0.46 | 0.48 | 0.88 | 17 | - |
| 2 | HATCN:PO-T2T 1:1 | 0.48 | 0.55 | 0.88 | 29 | - |
| 3 | HATCN:p-bPPhenB 1:1 | 0.79 | 1.00 | 0.92 | 56 | - |
| 4 | HATCN neat | 1.63 | 1 | 1.00 | (92) | HATCN domains span the film |
| - | HATCN:HAT-C6F5-CNPh-CNPh 1:1 (designed) | 1.49 | 1.00 | 0.92 | 39 | MW 971 |
| - | HATCN:HAT-CNPh x3 1:1 (designed) | 1.45 | 1.00 | 0.92 | 59 | partner crystallises |
| low | HATCN:NPB / TAPC / TPD / TCTA | 0.25-0.28 | ~0 | 0.80 | | CT absorption |
| low | Bphen:Cs2CO3 alone | 0.89 (Cs) | | 0.82-0.89 | | n-side, not a p-seed |

* TPBi stays first, but for a different reason than before: its multi-site
  PBE0 E_b (0.46) sits on the trap/diluent boundary, so it is a trap only half
  the time. It wins on morphology (risk 17), not on anchoring; the nucleation
  cost of 1:1 dilution is ~12 %.
* The amine HTLs are excluded by CT with HATCN (IPEA1 IP 6.29-6.31 eV vs >= 7.08
  for every ETL; the split is clean).
* A designed HAT partner makes every partner site a permanent trap (p = 1) at a
  lower risk than p-bPPhenB. Today's candidates are either too heavy or
  crystallise themselves; a lighter twisted-aryl HAT is the design target that
  follows from this.

## Transparency (TD-DFT, TDA-PBE0)

PBE0 places CT-like states lower than a range-separated hybrid, so "above 3.1 eV
at PBE0" is a conservative pass.

| molecule | neutral lowest (eV) | neutral bright below 3.1 eV | anion visible sum f |
|---|---|---|---|
| HATCN | 3.71 (dark; 3.96, 4.17 dark) | none | running |
| HAT-CN2benzo x3 | running | | |
| HAT-CN-F4benzo-F4benzo | running | | |

## What is solid and what is not

Solid: the library and routes; the cheap properties; the GFN2 site survey; the
finding that GFN2 ranks within a family but not across, and EA does not
predict E_b; the five PBE0 anchors (same protocol); the formulation logic with
current E_b.

Not solid: production-scale E_b assumes HATCN's 0.42 eV offset for all
analogues (one PBE0-relaxed derivative would test it); the aryl-HAT E_b is an
extrapolation of the HAT map; crystallisation risk is a heuristic score, not a
simulation; the CT rule is calibrated on known blends, not computed; the
TD-DFT is PBE0-TDA with a capped Davidson.

Next calculations that would change a decision, cheapest first:
1. PBE-D3 Ag relaxation + PBE0 for HAT-CN2benzo x3 (tests the 0.42 eV offset on
   the one-pot finalist).
2. PBE0 anchor for HAT-CNPh x3 (the aryl class; tests the map extrapolation).
3. Glass-formation check beyond the heuristic (MD of a 50-molecule film or
   literature Tg of HAT aryl derivatives).

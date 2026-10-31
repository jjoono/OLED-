# Common evaporable OLED materials ranked as Ag seed layers (final, 2026-10-13)

Question: with the screening rule fixed in `LIT_SEED_SYNTHESIS_20261011.md`
(gates first, then rank by E_b), where do the everyday HIL/HTL/host/ETL
materials fall, and does 3TPYMB work as a seed under thin Ag?

Data: `runs/common_ranking.json` (script 150), built from
`runs/common_survey_xtb.json` (script 148, GFN2 multi-site survey with
lone-pair starts), `runs/common_pbe0_fragment.json` (scripts 149 and 151) and
the earlier full-molecule PBE0 values.

## Method

1. **Sites (GFN2-xTB).** Each molecule built from SMILES, optimised with
   GFN2, and Ag relaxed (molecule frozen) from every heteroatom, ring hollow,
   lone-pair direction and heteroatom pocket.
2. **Score (PBE0-D3BJ/def2-SVP).** GFN2 E_b is not transferable across
   elements (PBE0/GFN2 = 0.19–0.68 on N sites, 0.7–0.85 on carbon sites), so
   the ranking value is PBE0: electronic part on a fragment (heavy atoms within
   5 Å, ring systems completed, substituents whole, H caps) plus D3BJ on the
   whole molecule. Validated on HATCN, Bphen, B3PyMPM, p-bPPhenB against
   full-molecule PBE0: largest error 0.013 eV (script 149).
3. **Crowded sites re-placed (script 151).** On sites between ortho
   hydrogens GFN2 pulls Ag to 1.7–2.2 Å from H, which PBE0 finds repulsive
   (T2T 0.07, Alq3 0.14 eV at the GFN2 position). For T2T, NBphen, DPEPO,
   PO-T2T and Alq3 the Ag position was re-found: start at the open point
   nearest the N/O anchors, relax Ag only at PBE-D3 on a 4 Å fragment, score
   at PBE0 as in step 2. Every crowded site moved to ≥ 2.6 Å from H except
   NBphen (2.19 Å, see below).
4. **Veto (GFN2).** Δq(Ag₁→Ag₂), sign validated against PBE0 earlier.

Every value in the table is PBE0 except four π-only materials (CBP, TCTA,
TPD, TAPC), which carry their motif's PBE0 value. The two π-only materials
computed directly confirm the motif values: mCP 0.267 (carbazole motif
0.279), NPB 0.254 (triarylamine motif 0.259).

## Ranking

Gates: condensation E_b ≥ ~0.5 eV; veto Δq ≥ +0.3 e. **No material fails the
veto** — every organic hands Ag₂ charge back or moves none — so the order is
set by E_b. E_b is split into its electronic part (fragment PBE0, no D3) and
dispersion (whole-molecule D3).

| # | material | role | E_b | electronic | dispersion | site | Δq (e) | class |
|---|---|---|---|---|---|---|---|---|
| 1 | HATCN | HIL / CGL | **1.63** | 1.46 | 0.18 | aza N | −0.20 | 1 acceptor |
| 2 | F4TCNQ | p-dopant | **1.21** | 1.15 | 0.06 | nitrile N | −0.27 | 1 |
| 3 | F6-TCNNQ | p-dopant | **1.20** | 1.10 | 0.10 | nitrile N | −0.26 | 1 |
| 4 | p-bPPhenB | ETL | **0.79** | 0.59 | 0.20 | phen chelate | −0.31 | 2 chelate |
| 5 | Bphen | ETL | **0.62** | 0.47 | 0.13 | phen chelate | −0.35 | 2 |
| — | *0.5 eV gate* | | | | | | | |
| 6 | PO-T2T | ETL / exciplex acceptor | 0.48 | 0.06 | 0.42 | groove at triazine N | −0.23 | 3* |
| 7 | BCP | HBL / ETL | 0.47 | 0.29 | 0.18 | phen N, off-chelate | −0.37 | 3 |
| 8 | TPBi | ETL / HBL | 0.46 | 0.22 | 0.24 | benzimidazole N3 | −0.21 | 3 |
| 9 | NBphen | ETL | 0.46 | 0.18 | 0.28 | half-blocked chelate | −0.32 | 3 |
| 10 | DPEPO | blue host / HBL | 0.45 | 0.17 | 0.28 | P=O | −0.51 | 3 |
| 11 | Alq3 | ETL / green host | 0.43 | 0.19 | 0.24 | quinolinolate O | −0.28 | 3 |
| 12 | B3PyMPM | ETL | 0.42 | 0.19 | 0.22 | pyridyl N | −0.24 | 3 |
| 13 | Liq | EIL | 0.40 | — | — | Li–O | −0.07 | 3 |
| 14 | TSPO1 | blue host / HBL | 0.38 | 0.27 | 0.11 | P=O | −0.02 | 3 |
| 15 | **3TPYMB** | ETL / HBL | **0.37** | 0.28 | 0.09 | pyridyl N | −0.13 | 3 |
| 16 | TmPyPB | ETL / HBL | 0.37 | 0.28 | 0.09 | pyridyl N | −0.11 | 3 |
| 17 | T2T | ETL / host | 0.32 | 0.08 | 0.25 | groove at triazine N | −0.23 | 4* |
| 18 | CBP | host | (0.28) | — | — | ring | −0.00 | 4 π-only |
| 19 | TCTA | HTL / host | (0.28) | — | — | ring | −0.03 | 4 |
| 20 | mCP | host / EBL | 0.27 | 0.04 | 0.22 | ring | +0.01 | 4 |
| 21 | TPD | HTL | (0.26) | — | — | ring | −0.02 | 4 |
| 22 | TAPC | HTL / EBL | (0.26) | — | — | ring | +0.02 | 4 |
| 23 | NPB | HTL | 0.25 | −0.00 | 0.26 | ring | −0.02 | 4 |

Parentheses: motif value (carbazole / triarylamine), confirmed by mCP and NPB.
\* see "dispersion-held sites".

## Reading the table

**Only five materials clear the gate, and they are two motifs.** Nitrile/aza
acceptors (1.2–1.6 eV; Ag donates into a deep LUMO, charge bounded at one
electron per molecule) and unhindered phenanthroline chelates (0.6–0.8 eV;
two N lone pairs at the right bite). Of these only HATCN, Bphen and
p-bPPhenB are practical neat seed films.

**Most ETLs and blue hosts sit in one band, 0.37–0.48 eV, just below the
gate.** Their spread (0.11 eV) is close to what the method resolves; their
order inside the band should not be read as a prediction until it is tested.

**HTLs and carbazole hosts are π-only (0.25–0.28 eV), benzene-like.** The
amine N of a triarylamine is buried; NPB's electronic binding is zero and all
of its E_b is dispersion. Expected to dewet like inert aromatics.

**Device role does not predict the class.** "ETL" spans Bphen 0.62 to T2T
0.32. Electron transport tracks the LUMO; Ag binding tracks whether a lone
pair is exposed.

### Substituents next to the binding atom decide

- **2,9-substituted phenanthrolines lose the chelate.** Bphen 0.62 → BCP 0.47
  (methyls), NBphen 0.46 (naphthyls). In NBphen Ag gets to 2.42 Å from one N
  but 2.85 Å from the other and touches a naphthyl H at 2.19 Å; the electronic
  part halves (0.47 → 0.18).
- **Aryl-triazines hide their N.** In T2T and PO-T2T each triazine N sits
  between two ortho hydrogens of the attached phenyls. From two different
  starts T2T's Ag ends 2.7–2.9 Å from N, next to H; electronic binding
  0.06–0.08 eV. The triazine is not a site in these molecules.
- **3TPYMB's boron is unreachable.** Ortho-methyls shield it; Ag never comes
  within 5 Å. The pyridyl N binds like a bare pyridine (0.374 vs 0.373 eV).

### Dispersion-held sites

PO-T2T ranks 6th on E_b but 0.42 of its 0.48 eV is dispersion: Ag sits in a
groove with three H at 2.7 Å. The electronic part (0.06) is the lowest of
any heteroatom material. Ranked by the electronic part alone, class 3 orders
differently:

BCP 0.29 · TmPyPB 0.28 · 3TPYMB 0.28 · TSPO1 0.27 · TPBi 0.22 · B3PyMPM 0.19
· Alq3 0.19 · NBphen 0.18 · DPEPO 0.17 · T2T 0.08 · PO-T2T 0.06

Which order matters on a real film is open. On a single molecule a groove
offers more dispersion than an open lone pair; on an amorphous film every
surface site has neighbours, so dispersion is supplied everywhere and the
electronic part is what differs between sites. The rule keeps E_b as the
ranking variable (it is what was validated against the literature), but
inside class 3 the two orders disagree, and **class 3 should be treated as
one tier rather than ranked internally.** The electronic-part order is a
specific prediction to test (e.g. BCP vs PO-T2T: same E_b, 5x different
electronic binding).

## 3TPYMB

- Site: pyridyl N, all three equivalent; boron inaccessible.
- E_b = 0.37 eV (electronic 0.28, dispersion 0.09). Same as TmPyPB and bare
  pyridine; below B3PyMPM (0.42), TPBi (0.46), BCP (0.47); well below Bphen
  (0.62) and the gate.
- Δq Ag₁→Ag₂ = −0.13 e: passes the veto; Ag₁ nearly neutral — a weak dative
  bond, no charge transfer.
- Expected: Ag condenses, nucleates sparsely, Volmer–Weber islands, later
  closure than on Bphen or HATCN. Usable as an ETL under a thicker Ag or with
  a seed on top (HATCN, Bphen, or metal-doped Ag); not a seed by itself.

## Method lessons from this round

- **GFN2 positions are not safe on crowded sites.** Single points at GFN2 Ag
  positions are within 0.03 eV of PBE0 on open sites (DMABN, Liq, Cu4I4) but
  are repulsive at PBE0 where GFN2 has packed Ag against hydrogens. Check the
  shortest Ag–H distance (< 2.4 Å → re-place) before trusting a score.
- **The GFN2-deepest site is not always the PBE0-deepest.** TmPyPB: GFN2
  preferred a bridge near N (PBE0 0.31); the lone-pair site is 0.37.
- **Locate at PBE-D3, score at PBE0.** PBE0 gradients on 40-atom fragments
  took 1.5–2 h per step; PBE-D3 took 5–40 min.

## Caveats

- One molecule, one Ag atom. Real films are amorphous; how many binding
  motifs face the surface matters as much as E_b, and bulky molecules (3TPYMB,
  TAPC) may bury theirs.
- Absolute E_b is PBE0/def2-SVP; for HATCN it drops 0.14 eV at
  TZVP+counterpoise. Differences < 0.05 eV between neighbours are not
  resolved, and class 3 as a whole spans only 0.11 eV.
- The 0.5 eV gate is a heuristic placed between the literature's good and
  poor seeds; BCP, TPBi, NBphen and PO-T2T at 0.46–0.48 are borderline.
- n-doped ETLs (Bphen:Cs, Liq co-deposited, Yb) change the surface chemistry
  and are not covered by this neat-film ranking.

## n-doped Bphen (Bphen:Cs2CO3) — script 153 (2026-10-28)

Model: Cs gives its electron to the phenanthroline and sits as Cs⁺ in the
chelate pocket (GFN2: q(Cs) +0.57, Cs–N 2.92 Å). At 20 wt% Cs2CO3 that is
about one Cs per two Bphen, so the doped surface mixes this complex with
neat Bphen. The whole-molecule PBE0 on Bphen–Cs (44 atoms with Ag) did not
finish one SCF in 2 h on 4 cores, so PBE0 was run on phen–Cs (4,7-phenyls →
H; GFN2 reproduces Bphen–Cs: q(Cs), Cs–N, site order, Δq) with a neat-phen
chelate baseline in the same model. Whole-Bphen Gaussian jobs with PBE0-
optimised Ag positions are in `bphencs_jobs/` (script 154).

| site (phen model) | PBE0 E_b | electronic | dispersion | state | q(Ag) |
|---|---|---|---|---|---|
| neat phen, chelate (baseline) | 0.45 | 0.31 | 0.13 | doublet | +0.01 |
| **phen⁻Cs⁺, Ag next to Cs** | **0.87** | 0.84 | 0.03 | singlet (triplet +0.93 eV) | **−0.45** |
| phen⁻Cs⁺, Ag at the N side | −0.04 (unbound at the GFN2 position) | −0.25 | 0.20 | triplet | −0.07 |

GFN2 Δq(Ag₁→Ag₂) at the Cs site: −0.01 (phen–Cs), −0.02 (Bphen–Cs).

Reading:
- **Doping adds a new, stronger site and removes the old one.** The radical
  anion's extra electron pairs with Ag's 5s electron next to Cs⁺ (closed-shell
  singlet, 0.93 eV below the triplet): a real chemical bond, +0.43 eV over
  the neat chelate in the same model. Cs⁺ occupies the chelate pocket, and
  the N side is no longer a site.
- **Charge flows the other way from HATCN/MoOx.** Ag takes electrons
  (−0.45 e at PBE0); Ag₂ keeps the same charge — no veto; the Ag stays
  metallic/electron-rich.
- **Expected place in the ranking.** The Cs site carries almost no
  dispersion (0.03 eV), so the 4,7-phenyls should change it little: about
  0.9 eV on Bphen–Cs (to be confirmed by `bphencs_jobs`). The doped surface
  then offers ~0.9 eV sites on the doped half of the molecules and the 0.62 eV
  chelate on the rest: above neat Bphen (0.62), Cs2CO3 alone (0.70) and
  p-bPPhenB (0.79), below F4TCNQ (1.21) and HATCN (1.63).
- **Prediction:** Ag on Bphen:Cs2CO3 closes at a smaller thickness than on
  neat Bphen, but not as early as on HATCN. In the stack
  Bphen:Cs2CO3 / HATCN / Ag the HATCN remains the better seed by E_b.

Caveats: phen model; Ag positions are GFN2 single points (the N-side value
in particular may improve if relaxed at PBE0 — the Gaussian jobs do this);
Cs2CO3 partly decomposes on evaporation, and Cs is mobile.

### Whole-Bphen confirmation (script 155, 2026-10-31)

Cs-keeping fragment PBE0 (phen core + Cs⁺ + caps) with whole-molecule D3,
Ag located at PBE-D3 on a 4 Å fragment (replaces the Gaussian package, which
did not run on the workstation).

| site on Bphen⁻Cs⁺ | PBE0 E_b | electronic | dispersion | state | q(Ag) | Ag–Cs |
|---|---|---|---|---|---|---|
| **next to Cs (Cs4)** | **0.89** | 0.86 | 0.03 | singlet | −0.44 | 3.57 Å |
| N side (hollow8), relaxed | 0.51 | 0.32 | 0.19 | singlet | −0.39 | 3.70 Å |
| π face, far from Cs (hollow3)* | 0.25 | 0.09 | 0.16 | triplet | −0.05 | 8.82 Å |

\* relaxation stopped at its lowest point (|g| 0.14 eV/Å) after the line
search looped.

PBE0 Ag₂ at the Cs site: total charge −0.46 against −0.44 for Ag₁ —
Δq(Ag₁→Ag₂) = −0.015 e, passes the veto (GFN2 gave −0.02).

- The phen-model estimate is confirmed: 0.89 eV (phen 0.87). The phenyls add
  0.02 eV to the Cs site.
- Every bound site draws its strength from Cs: the N-side minimum moves
  towards Cs (3.70 Å) and takes electrons; far from Cs the surface is
  π-only (0.25, unpaired). So on a doped film the useful sites are the Cs
  neighbourhoods, and their density is set by the Cs content.
- Place in the ranking: neat Bphen 0.62 < Cs₂CO₃ 0.70 < p-bPPhenB 0.79 <
  **Bphen⁻Cs⁺ 0.89** < F4TCNQ 1.21 < HATCN 1.63.

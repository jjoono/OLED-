# Common evaporable OLED materials ranked as Ag seed layers (2026-10-12)

Question: with the screening rule fixed in `LIT_SEED_SYNTHESIS_20261011.md`
(gates first, then rank by E_b), where do the everyday HIL/HTL/host/ETL
materials fall, and does 3TPYMB work as a seed under thin Ag?

Data: `runs/common_ranking.json` (script 150), built from
`runs/common_survey_xtb.json` (script 148, GFN2 multi-site survey with
lone-pair starts) and `runs/common_pbe0_fragment.json` (script 149, fragment
PBE0 at the GFN2-deepest site).

## Method in one paragraph

Each molecule was built from SMILES, optimised with GFN2-xTB, and searched for
Ag sites (every heteroatom, every ring hollow, lone-pair directions, pockets
between heteroatoms). GFN2 finds the sites well but its E_b is not
transferable across elements (PBE0/GFN2 = 0.19–0.68 on N sites, 0.7–0.85 on
all-carbon sites), so the ranking value is PBE0: the electronic part on a
fragment (heavy atoms within 5 Å of Ag, ring systems completed, substituents
kept whole, H caps) plus D3BJ on the whole molecule. Validated against
full-molecule PBE0 on HATCN, Bphen, B3PyMPM and p-bPPhenB: largest error
0.013 eV. The charge-capacity veto uses GFN2 Δq(Ag₁→Ag₂), whose sign was
validated against PBE0 earlier. Materials not yet run at PBE0 carry the PBE0
value of their closest computed motif and are marked *analog*.

## Ranking

Gates: condensation E_b ≥ ~0.5 eV; veto Δq ≥ +0.3 e. **No material in the list
fails the veto** — every organic hands Ag₂ charge back or moves none — so the
order is set by E_b alone.

| # | material | role | E_b (eV) | site | Δq Ag₁→Ag₂ (e) | source |
|---|---|---|---|---|---|---|
| 1 | HATCN | HIL / CGL | 1.63 | aza N | −0.20 | PBE0 |
| 2 | F4TCNQ | p-dopant | 1.21 | nitrile N | −0.27 | PBE0 |
| 3 | F6-TCNNQ | p-dopant | 1.20 | nitrile N | −0.26 | PBE0 fragment |
| 4 | p-bPPhenB | ETL | 0.79 | phen chelate | −0.31 | PBE0 |
| 5 | Bphen | ETL | 0.62 | phen chelate | −0.35 | PBE0 |
| — | *0.5 eV condensation gate* | | | | | |
| 6 | NBphen | ETL | (0.47) | phen pocket | −0.32 | analog — running |
| 7 | BCP | HBL / ETL | 0.47 | N, off-chelate | −0.37 | PBE0 fragment |
| 8 | TPBi | ETL / HBL | 0.46 | benzimidazole N3 | −0.21 | PBE0 fragment |
| 9 | B3PyMPM | ETL | 0.42 | pyridyl N | −0.24 | PBE0 |
| 10 | TmPyPB | ETL / HBL | (0.42) | pyridyl N | −0.11 | analog (B3PyMPM) |
| 11 | DPEPO | blue host / HBL | (0.41) | P=O | −0.51 | analog (Ph₃P=O) |
| 12 | TSPO1 | blue host / HBL | (0.41) | P=O / ring | −0.02 | analog (Ph₃P=O) |
| 13 | PO-T2T | ETL / exciplex acceptor | (0.41) | P=O / triazine | −0.23 | analog — running |
| 14 | Liq | EIL | 0.40 | Li–O | −0.07 | PBE0 |
| 15 | Alq3 | ETL / host | (0.40) | quinolinolate O | −0.28 | analog (Liq) |
| 16 | 3TPYMB | ETL / HBL | 0.37 | pyridyl N | −0.13 | PBE0 fragment |
| 17 | CBP | host | (0.28) | ring hollow | −0.00 | analog (PhCz) |
| 18 | mCP | host / EBL | (0.28) | ring hollow | +0.01 | analog (PhCz) |
| 19 | TCTA | HTL / host | (0.28) | ring hollow | −0.03 | analog (PhCz) |
| 20 | T2T | ETL / host | (0.27) | triazine N | −0.23 | analog (triazine) |
| 21 | NPB | HTL | (0.26) | ring hollow | −0.02 | analog (TPA) |
| 22 | TPD | HTL | (0.26) | ring hollow | −0.02 | analog (TPA) |
| 23 | TAPC | HTL / EBL | (0.26) | ring hollow | +0.02 | analog (TPA) |

Values in parentheses are analogs and will be replaced as the PBE0 queue
finishes (NBphen, PO-T2T, T2T, TmPyPB, DPEPO, Alq3, TSPO1).

## Reading the table

**Four classes, separated by the binding motif, not by the role in the device.**

1. *Nitrile/aza acceptors (1.2–1.6 eV)* — HATCN, F4TCNQ, F6-TCNNQ. Ag donates
   into a deep LUMO; charge stays bounded at one electron per molecule. Only
   HATCN is a practical neat seed (the TCNQs are dopants, rarely used neat, and
   F4TCNQ is volatile and diffuses in films).
2. *Phenanthroline chelates (0.6–0.8 eV)* — Bphen, p-bPPhenB. Two N lone pairs
   held at the right bite angle. The only ETL-class materials above the gate.
3. *Single lone-pair donors (0.37–0.47 eV)* — pyridyl, benzimidazole,
   phosphine oxide, metal-quinolinolate O. Ag binds one lone pair; the 2–3
   rings around it add 0.05–0.1 eV through dispersion and polarisation. This is
   where most ETLs and blue hosts sit, just below the gate.
4. *π-only surfaces (0.26–0.28 eV)* — carbazole and triarylamine HTLs/hosts,
   triazine without exposed N. Ag sits over a ring; benzene-like. Expected to
   dewet like the inert aromatics in the literature.

**Substituents decide inside a class.** BCP is Bphen with 2,9-methyls; the
methyls block the chelate pocket and cost 0.15 eV (0.62 → 0.47), dropping it a
class. The same logic predicts NBphen (2,9-naphthyl) below Bphen; the PBE0
result will test it. In 3TPYMB the boron is shielded by ortho-methyls and Ag
never gets within 5 Å of it — the binding is that of a bare pyridine (0.374 vs
pyridine 0.373 eV).

**Role in the device does not predict the class.** "ETL" spans classes 2–4
(Bphen 0.62 … T2T ~0.27). Electron-transport ability tracks LUMO energy;
Ag binding tracks lone-pair geometry.

## 3TPYMB

* Site: pyridyl N (Ag–N 2.6 Å, GFN2), all three pyridyls equivalent.
  Boron inaccessible.
* E_b = 0.37 eV (PBE0 fragment). Below B3PyMPM (0.42), TPBi (0.46), BCP (0.47),
  well below Bphen (0.62) and the 0.5 eV gate.
* Δq Ag₁→Ag₂ = −0.13 e: passes the veto; Ag₁ is nearly neutral (−0.03 e GFN2,
  −0.11 Mulliken PBE0) — a weak dative bond, no charge transfer.
* Expected behaviour: Ag condenses but nucleates sparsely, Volmer–Weber islands,
  later closure than on Bphen or HATCN. Usable as an ETL under a thicker Ag or
  with a seed (HATCN, Bphen, or metal-doped Ag) on top; not a seed by itself.

## Caveats

* Single isolated molecule, one Ag atom. Real films are amorphous; the density
  of exposed binding motifs at the surface matters as much as E_b, and bulky
  molecules (3TPYMB, TAPC) bury their heteroatoms partially.
* Absolute E_b is at PBE0/def2-SVP. For HATCN it drops by 0.14 eV at
  TZVP+counterpoise (1.63 → 1.50); differences of < 0.05 eV between
  neighbours are not resolved.
* The 0.5 eV gate is a heuristic placed between the literature's good and poor
  seeds; materials at 0.45–0.5 eV (BCP, TPBi) are borderline, not failures.
* n-doped ETLs (Bphen:Cs, Liq co-deposited, Yb) change the surface chemistry
  and are not covered by this neat-film ranking.

# Seed screening re-audit, 2026-08-21

Prompted by a challenge to the Al2O3 entry. The challenge was correct, and
checking it turned up four more entries that should not have been in a single
ranked list.

## Retracted or reclassified

| entry | was | why it is wrong |
|---|---|---|
| **Al2O3 0.910 eV** | ranked 4th | Al-rich **Al10O10** cluster. `scripts/23_al4o6_ag.py` already records it as overbinding — undercoordinated Al with dangling bonds. The stoichiometric **Al4O6** cluster gives **0.422 eV**, which is MoOx territory, not top-five. |
| **MoO3 −0.718 eV** | quoted as repulsive | Unrefined: xTB geometry with a DFT single point. Relaxing the Ag position at DFT level gives **+0.284 eV**. Superseded by `psi4_binding_refined_eV.json`. |
| **CuI 0.819 / CuSCN 0.567** | ranked as materials | These are **Ag@Cu** sites — a silver-copper metallic contact, not the surface a seed layer presents. The relevant sites are Ag@I (0.222) and Ag@S (0.136). |
| **ZnS 1.60 eV** | in the kMC table | Not our calculation. A literature prior from Kim AFM 2015, carried for the growth model only. |
| **PhCN·ZnS 1.218 eV** | ambiguous | This is **ZnS binding to a nitrile**, i.e. whether ZnS sticks to HATCN. It is not an Ag binding energy at all. |

## The deeper problem — three model classes were being ranked together

**A. Whole-molecule DFT.** Real materials, mutually comparable.
HATCN (CN) 1.029 · F4TCNQ 0.966 · Cs2CO3 0.900 · TPBi 0.889 · p-bPPhenB 0.870 ·
B3PyMPM 0.626 · Bphen 0.490 · Liq 0.167 · HATCN (face) −0.024

**B. Inorganic cluster models.** Strongly model-dependent — composition and
coordination number move these by more than the spread between candidates.
MoOx Mo3O8 0.445 · Al2O3 Al4O6 0.422 · MoO3 Mo3O9 0.284 · LiF32 0.249 ·
CuI Cu4I4@I 0.222 · CuSCN@S 0.136

**C. Functional-group model compounds.** Not materials. Proxies for a binding
motif, useful for ranking chemistries, meaningless as a material property.
pyridine 0.342 · triazine 0.291 · PhCz→TCTA 0.279 · BTD 0.277 ·
TPA→TAPC 0.254 · Me3P=O→DPEPO 0.253 · thiophene 0.169 · benzene 0.167

Note that TPA and PhCz were being read as materials. `scripts/12_htl_sites.py`
declares them as proxies: "Ag on triphenylamine (TAPC proxy), N-phenylcarbazole
(TCTA proxy)". The HTL cores, not the HTLs.

## What survives

HATCN's position is unchanged, and the correction strengthens rather than
weakens it: removing the spurious Al2O3 0.910 takes away the only inorganic
that appeared to come close.

Among materials that can sit on the anode side, HATCN at 1.029 eV leads the
runner-up (oxygen-deficient MoOx, 0.445 eV) by **2.3x**. Everything ranked
between them is an electron-transport material that cannot be used there.

The pre-existing caveat in `notes/EB_EA_VALIDATION_AUDIT.md` still holds and
should be quoted with any of these numbers: single-adatom adsorption energies
have no experimental counterpart, so the paper should claim **ordering** with a
−0.2 to +0.3 eV model error, never absolute values. The cluster-versus-slab gap
for HATCN itself (1.03 vs 1.346 eV) is the scale of that error.


---

# Second pass — binding mode, 2026-08-21

Prompted by a challenge to TPBi at 0.889 eV. Checking it meant opening the
optimised geometries and asking where the silver actually sits, which turned
out to be a different question from which site it was placed at.

## TPBi is not a single-nitrogen site

It was **placed** at one benzimidazole nitrogen. After the GFN2 optimisation
that the DFT single point runs on, the silver has moved into a pocket and
coordinates **two** nitrogens:

| system | placed at | relaxed to |
|---|---|---|
| TPBi | N 2.30 Å | **N 2.33 / N 2.34 Å** — bidentate |
| p-bPPhenB | N 1.39 Å (bad guess) | **N 2.28 / N 2.39 Å** — bidentate |
| HATCN | N 2.30 Å | **N 2.30 Å, next atom C at 3.46 Å** — monodentate |

That resolves the puzzle: a single pyridine nitrogen gives 0.342 eV, so 0.889
for one nitrogen made no sense. Two coordination bonds do.

## Binding mode, read off the relaxed geometries

**Monodentate** — HATCN 1.029 (N 2.30) · F4TCNQ 0.966 (N 2.04) ·
B3PyMPM 0.626 (N 2.31) · pyridine 0.342 (N 2.59) · triazine 0.291 (N 2.24) ·
BTD 0.277 (N 2.09) · Me3P=O 0.253 (O 2.64) · thiophene 0.169 (S 2.90)

**Bidentate chelate** — Cs2CO3 0.900 (O 2.44/2.53) · TPBi 0.889 (N 2.33/2.34) ·
p-bPPhenB 0.870 (N 2.28/2.39) · Bphen 0.490 (N 2.28/2.28)

**Dispersion only** — benzene 0.167 (C 3.36)

## The comparison this makes possible

HATCN is monodentate and still beats every bidentate chelate in the set. One
nitrile nitrogen out-binds two benzimidazole nitrogens. And HATCN carries six
of those sites on a planar molecule, so all of them are surface-accessible,
against three buried pockets on TPBi's propeller.

Site energy and site density both favour HATCN, and they are the two terms that
set nucleation density. This is a stronger statement than "HATCN has the
highest E_b" and should replace it.

## Three more values that do not mean what they claim

| entry | what the geometry shows |
|---|---|
| **PhCz → TCTA proxy, 0.279** | Silver relaxes to C 2.74 / H 2.76 — it is on the π face, not the nitrogen. Carbazole nitrogen is trisubstituted and has no available lone pair, so this is correct chemistry but the number is a dispersion floor, not a site energy. |
| **TPA → TAPC proxy, 0.254** | Silver sits over three hydrogens at 2.89 Å. Tertiary amine, same reason. |
| **Liq, 0.167** | Silver ends up beside the lithium at 2.46 Å, not in the intended O,N chelate. |

The first two are worth keeping as a **result** rather than a ranking entry:
hole-transport cores cannot anchor silver at all, which is why an HTL by itself
is useless as a seed. The Liq number should not be quoted.

## Checked and sound

B3PyMPM's `structures/` file has silver 1.85 Å from a hydrogen, but that is the
pre-optimisation guess; the geometry the DFT actually used has N at 2.31 Å.
Bphen, Cs2CO3, F4TCNQ, benzene, pyridine, triazine, BTD, thiophene and Me3P=O
all relax to chemically sensible contacts.

There is no size artefact: correlation between E_b and molecule size is +0.22,
and against the number of atoms within 5 Å it is **negative** (−0.39), so
dispersion from bulk is not inflating the larger molecules.

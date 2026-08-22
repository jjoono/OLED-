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

# HATCN co-deposition partner ranking (from lab logbook, 2026-08-14)

Source: `data/systems_logbook_Cologne_260804.xlsx` — sheets 'Materials in stock'
(80 entries), 'Materials Database' (82), 'Creaphys Materials' (19 lots).
Only IN-STOCK materials are ranked; database-only entries are listed at the end.

## Selection criteria (absorption-first, per project priority)

The mixed seed layer must suppress HATCN crystallization WITHOUT adding an
absorption channel. That fixes four filters:

1. **No donors.** HATCN is a strong p-dopant (EA ~4.8-5.2 eV). Any partner with
   HOMO shallower than ~ -6.0 eV gets p-doped in the mixed film → CT/polaron
   bands in the visible-NIR. This kills every arylamine/carbazole-donor host.
2. **UV-only absorber.** Optical gap > 3.1 eV so the partner itself is
   transparent at 400-800 nm. Kills all emitters and dyes.
3. **Glass former with 3D shape.** The whole point is disrupting HATCN's
   planar stacking; bulky tetrahedral/star-shaped molecules with high Tg win.
4. **Evaporable & stable** at normal rates (all listed materials qualify —
   they are OLED evaporants by definition).

Bonus criterion from the DFT screen: pyridine-N and P=O groups bind Ag
(pocket effect +0.53 eV), so a partner carrying them can ADD anchoring
density instead of merely diluting HATCN.

## RANKING — in stock, best first

| # | Material | HOMO (eV) | Gap | Tg (°C) | Why this rank |
|---|----------|-----------|-----|---------|---------------|
| 1 | **TPBi** | -6.2..-6.7 | 3.5 | ~122 | Our own crystallization screen scored it 17/100 (HATCN 92) — best glass former on the list. Benzimidazole N binds Ag. Deep HOMO → no CT with HATCN. Massive stock (7 Creaphys lots). The obvious first experiment. |
| 2 | **PO-T2T** | ~ -7.5 | 3.5 | ~108 | Deepest HOMO in the lab → physically zero CT risk. Strongest pure acceptor. Three P=O groups = strong Ag anchors (same chemistry as our TSPO1 finding). Star-shaped, good glass. |
| 3 | **3TPYMB** | ~ -6.8 | 3.5 | ~108 | Tetrahedral boron center — maximal shape disruption of HATCN stacking. Three pyridines for Ag anchoring. Transparent, deep HOMO. |
| 4 | **B3PYMPM** | ~ -6.8 | 3.6 | ~107 | Pyridine/pyrimidine acceptor, transparent, deep HOMO. Slightly lower rank: tends to pack/orient (H-bonded layers) rather than form a fully isotropic glass. |
| 5 | **TmPPPyTz** | ~ -6.8 | ~3.6 | high | Triazine+pyridine ETL, same class as PO-T2T/B3PYMPM. Less literature on film morphology → rank below the known quantities. |
| 6 | **NBPhen** | ~ -6.4 | 3.2 | ~105 | Phenanthroline N binds Ag well (BPhen/Ag is a classic contact). NBPhen fixes BPhen's crystallization problem. Gap 3.2 eV is the narrowest in this tier — absorption edge nears 400 nm. |
| 7 | **TmPyPB** | ~ -6.7 | 3.9 | ~79 | Chemically ideal (deep HOMO, very wide gap, 3 pyridines) but Tg 79 °C is low — a marginal glass former is a weak crystallization suppressant. |
| 8 | **DPEPO** | ~ -6.1..-6.6 | 4.0 | ~93 | Widest gap available; two P=O Ag anchors. Downsides: thermally fragile (decomposes near evaporation T if run hot), borderline HOMO. |
| 9 | **BPhen** | ~ -6.4 | 3.2 | ~62 | Good Ag chemistry but itself crystallizes notoriously (Tg 62). Only if NBPhen runs out. |

Borderline, use only with a control: **mCPPO1, mCPCN, PCzAc, CzSi-class** —
the P=O/nitrile side is fine but the carbazole HOMO (~ -5.9..-6.1 eV) sits at
the edge of HATCN's doping window; a weak CT tail in the blue is possible.
Check a mixed-film absorption spectrum before committing.

## EXCLUDED (in stock), with reasons

- **Donors — CT absorption with HATCN guaranteed** (HATCN p-dopes them; this
  is literally how HATCN is used as HIL): NPB, TAPC, X-F6-TAPC, TCTA, CBP,
  mCBP, mCP, Spiro-TAD, Spiro-TTB, TSBF, Rubrene.
- **Visible absorbers / emitters** (violates absorption-first): 4CzIPN,
  5TCzBN, C545T, DCM, DBP, TBPe, TBRb, TTPA, BDAVBi, v-DABNA, BSBCz, CzDBA,
  DMAC-TRz, SpiroAC-Trz, SubPc, SubNc, Cl6SubPc, C60, C70, DCNP, Super
  Yellow, F8BT, Ir(MDQ)2(acac), Ir(ppy)2(acac), Pt(TPBP), MADN (anthracene
  edge ~400 nm + crystalline).
- **F6-TNAP**: F4TCNQ-class strong dopant — Ag-salt risk at the metal
  interface plus NIR polaron bands when it dopes anything nearby. Same logic
  that excluded F4TCNQ as a seed.
- **Alq3 / BAlq**: absorption edge extends past 400 nm (green PL); also Al
  chelates offer no acceptor character advantage.
- **Inorganics** (not organics, listed for completeness): MoO3, WO3, LiF,
  LiQ, Cs2CO3, Cs.

## In DATABASE but not confirmed in stock

PPT (would rank ~top-3: P=O, very wide gap), TpPyPB, DPy/2DPyS/3PyS-class,
OXD family (oxadiazole acceptors — would rank mid-tier), F6-TCNNQ (project
anchor molecule, but excluded for co-deposition by the same dopant logic as
F6-TNAP). If PPT turns out to be physically present, slot it at #2-3.

## Suggested first experiment

HATCN:TPBi co-deposition at 3:1 and 1:1 (volume), 15-20 nm total, then Ag
5 nm at 0.5 nm/s (or two-step). Readout: relative T + absorption estimate +
sheet resistance vs the pure-HATCN sample 5 baseline. TPBi first because it
is the only partner our own risk screen has already quantified (17 vs
HATCN's 92), and stock is deepest. PO-T2T 1:1 as the second arm — it tests
whether P=O anchoring adds nucleation density on top of dilution.

## ANODE-SIDE variant: HTL/HIL-only ranking (user request, 2026-08-14)

When the mixed seed sits on the ANODE side, the partner must come from the
hole-side materials. That makes some CT absorption unavoidable (HATCN p-dopes
donors -- that is its HIL mechanism), so the ranking metric flips to
"deepest HOMO = weakest CT", then glass formation:

1. mCBP   (HOMO ~ -6.0, Tg ~95)  -- edge of HATCN's doping window, weak CT
   tail at worst; biggest Creaphys stock. Best absorption compromise.
2. TCTA   (HOMO ~ -5.7, Tg ~151) -- clear CT expected, but the strongest
   glass former / crystallization suppressant in the lab. Control arm.
3. mCPCN / mCPPO1 (HOMO ~ -6.0)  -- mCBP-class HOMO + CN/P=O Ag anchors.
4. mCP    (HOMO ok, Tg ~60: itself marginal against crystallization)
5. CBP    (HOMO ok, but planar + notorious crystallizer -- self-defeating)
6. TSBF   (spiro shape ideal; properties unverified -- promote if confirmed)
7. TAPC, 8. NPB, 9. Spiro-TAD/TTB -- HOMO -5.5..-5.3, textbook HATCN CT
   pairs, visible-NIR polaron/CT absorption. Not recommended.

Injection note: HATCN charge generation (electron extraction from the
adjacent HTL's HOMO into HATCN's LUMO) survives dilution; keep the mix
HATCN-rich (3:1) rather than 1:1 on the anode side.

Cheap pre-check: absorption spectrum of glass / HATCN:mCBP 20 nm mixed film
BEFORE any Ag work -- a CT band shows as a broad visible-NIR hump.

# Design spec: Ag seed on the p-side of a charge-generation junction

Step 1 of the design programme. The role is fixed first; the ranking (scripts
160-165) only asks which molecule does this one job better than HATCN.

## The role

```
 ... / n-ETL (Bphen:Cs2CO3, or Bphen:Liq) / SEED 3-5 nm / Ag 8 nm (top, semi-transparent)
```

The seed does three jobs at once:

1. **Electrical**: the p-side of the n/p charge-generation junction. Its LUMO
   has to sit deep, close to the Fermi level of Ag and of the n-doped ETL,
   so that electrons cross the Ag/seed/n-ETL stack without a barrier. HATCN
   does this because its electron affinity is about 5 eV (solid). A seed whose
   LUMO is 0.5 eV shallower puts a barrier into the cathode.
2. **Growth**: Ag on it has to nucleate densely and close early. From the
   screen: E_b(PBE0) at the deepest site ≥ 1.0 eV (HATCN 1.63, F4TCNQ 1.21;
   ETL-class ≤ 0.8), and the Ag1→Ag2 charge must not grow (Δq ≤ 0, the
   veto that MoOx fails).
3. **Optical**: it sits inside the top electrode of a transparent/top-emitting
   device, so it may not absorb in the visible, **including as the anion**:
   the first monolayer is reduced by Ag (HATCN takes ~0.3-0.5 e), and an
   anion band in the visible is a loss that a neutral spectrum hides. The
   F4TCNQ anion's 2.4 eV band (f = 0.20) is the warning case.

And two practical constraints: it is **thermally evaporated** (sublimes
intact) and it has to be **made** in a chemistry lab within months.

## Thresholds (all relative to HATCN computed the same way)

| # | property | how computed | threshold | HATCN | why |
|---|---|---|---|---|---|
| 1 | EA | IPEA1-xTB vertical (script 161) | ≥ HATCN − 0.3 eV = **3.39 eV** | 3.69 | CGL p-side alignment; 0.3 eV ≈ what a dipole/doping change can absorb |
| 2 | E_b | GFN2 site survey (162) → PBE0 via Δ-ML (163), confirmed by fragment PBE0 for the finalists | **≥ 1.0 eV** (predicted − 1σ for "safe") | 1.63 | dense nucleation; ETL-class materials sit ≤ 0.8 |
| 3 | Δq(Ag1→Ag2) | GFN2, Ag2 at the deepest site | **≤ 0** | −0.23 | the veto |
| 4 | optical (screen) | GFN2 HOMO-LUMO gap | ≥ HATCN − 0.35 eV (= 0.83 eV) | 1.18 | cheap prefilter only |
| 4b | optical (finalists) | TDA-ωB97X/def2-SVP, neutral and anion | neutral: no state with f > 0.02 below 3.1 eV (400 nm); anion: Σf(1.65-3.1 eV) ≤ HATCN anion's | neutral 3.92 eV; anion 1.96 eV f 0.15×2 | transparency inside the electrode |
| 5 | sublimability | MW, HBD | 300 ≤ MW ≤ 900, HBD = 0 | 384 | evaporable without decomposition |
| 6 | synthesis | route class + SA score | route class 1 (homoleptic one-pot) or 2 (mixture, separable); SA ≤ 4.5 | class 1 | made in-house |
| 7 | morphology | crystallisation risk (script 53: flat + rigid + symmetric + small) | **objective, lower is better** | 92 | the Park & Suh failure mode: HATCN crystallises at 7 nm |

Hard constraints: 1, 3, 5, 6 (and 4 as prefilter). Objectives for the Pareto
front: maximise E_b, minimise crystallisation risk, maximise gap, minimise
synthesis cost. 4b, the fragment PBE0 E_b and the TD-DFT go to the finalists
only.

## What HATCN loses on

Only on axis 7. HATCN is the best on E_b, it is EA-aligned and transparent;
its weakness is that it is a small, flat, D3h disc that crystallises. So the
design question is sharp: **keep EA and E_b at HATCN level and break the
crystal packing**, without opening a visible absorption. That is what the
periphery substitutions in the library (fused benzo with CN/F/Cl, route-B
aryls, mixed units that break D3h) are for.

Machine-readable copy: `runs/design/spec.json`.

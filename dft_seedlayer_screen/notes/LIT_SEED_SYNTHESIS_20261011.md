# Seed-layer literature synthesis and the priority of screening descriptors (2026-10-11)

Question: putting the lab's own result (HATCN/Ag continuous at 7 nm; MoOx/Ag
percolated with 9–23 % voids at 8 nm) beside everything reported for Ag seed
layers, which descriptor should rank candidates, and in what order?

Sources are web abstracts, open full texts and theses found on 2026-10-11 (list
at the end). **Numbers come from different groups, deposition methods and
definitions of "percolation" (optical, conductance, mass-equivalent); they are
not head-to-head.** Ordering within a class is reliable only where one paper
compares directly.

## 1. What the literature reports

| class | seed (thickness) | reported Ag result | mechanism offered | direct comparison? |
|---|---|---|---|---|
| none | glass / SiO₂ | evaporated: not closed, not conductive < 10 nm; sputtered: percolation ~6 nm, continuous ~12 nm | Volmer–Weber, weak Ag–oxide adhesion | baseline |
| metal | Cu (1 nm) | percolation 3 nm (vs 6 nm without), RMS < 0.5 nm, layer-like growth | Ag–Cu bonding, smaller grains | yes (vs no seed) |
| metal | Ge (1–2 nm) | RMS improved >10×; 10 nm Ag 22 vs 51 Ω/sq; FM-like growth | Ag–Ge > Ag–Ag | yes (vs SiO₂) |
| metal | Ca, Al, Au (1 nm) | wetting improves *with seed surface energy* | surface energy | yes (three seeds) |
| metal dopant | Al in Ag (co-dep.) | 6–7 nm continuous, sub-nm roughness, 80 % T at 7 nm | Al suppresses 3D islanding | yes (vs pure Ag) |
| metal dopant | Cu in Ag | percolation 6 nm, RMS 0.19 nm; 6.5 nm DMD at 100.3 % relative T | same | yes |
| gas dopant | O₂ / N₂ in early Ag | suppresses VW growth; 13–18 Ω/sq at 86–87 % T | surfactant | yes |
| oxide, weak acceptor | **CuO (0.25–1 nm)** | **percolated near 1 nm** on SiO₂ and CaF₂; transparent | not stated | yes (vs bare) |
| oxide, weak acceptor | ZnO (5 nm on TiO₂) | 8 nm Ag continuous where a-TiO₂ is not; 5.68 vs 7.56 Ω/sq | ZnO(0001)‖Ag(111) texture, fewer grain boundaries | yes (vs TiO₂) |
| oxide | TiO₂ (amorphous) | 8 nm Ag discontinuous, abnormal grains | — | yes (vs ZnO) |
| oxide, deep acceptor | **MoO₃ (1 nm)** | **7 nm percolates but porous, not islanded**; 9 nm: 20 Ω/sq, 74 % T | nucleation layer | yes (vs bare) |
| chalcogenide | ZnS | continuous where TiO₂ underlayer is not | Ag–S | yes (vs TiO₂) |
| polymer donor | PEI (+Au) | 8 nm Ag 9 Ω/sq, 93 % T (400–800 nm) | amine nucleation inducer | partial |
| SAM | thiol, primary amine | promote Ag nucleation | high Ag affinity | qualitative |
| organic, inert | perfluorinated polymers/SAMs | condensation coefficient C ≈ 0 above ~5–10 nm | no reactive moiety, low polarisability | yes |
| organic | Bphen / BCP | Ag forms [Ag(Bphen)]⁺, [Ag(Bphen)₂]⁺ — binds Ag⁺, not Ag⁰ | coordination + n-doping | — |
| organic | HATCN (7 nm) | 15 nm Ag still has voids/cracks (Park & Suh 2018) | used as HIL | — |
| **this lab** | HATCN | **continuous at 7 nm** | — | yes (vs MoOx) |
| **this lab** | MoOx | **percolated with 9–23 % voids at 8 nm** | — | yes (vs HATCN) |

Two independent observations of MoO₃ — the literature's "porous instead of
islanded at 7 nm" and this lab's "percolated with voids at 8 nm" — describe the
same morphology: **many nuclei, connected early, never closing.**

## 2. Which descriptor explains what

| descriptor | explains | fails on |
|---|---|---|
| **A. Ag–seed bond strength** (seed surface energy; molecular E_b) | metals ≫ oxides ≫ inert organics; Ca/Al/Au ranking; fluorinated C ≈ 0; thiol/amine/PEI work; Ge/Cu > bare | **MoO₃ vs HATCN** (MoO₃ binds harder by 0.45 eV at TZVP+CP, yet loses); MoO₃'s early percolation but porosity |
| **B. Electron capacity under a growing cluster** (computed: total Ag charge bounded vs extensive) | MoO₃ porous (literature + lab); HATCN, Bphen, F4TCNQ return charge from Ag₂; Bphen binds a *lone* Ag⁺ (= our n = 1 result); CuO and ZnO — weak acceptors (EA ~4.1–4.3 eV vs MoO₃ ~6.7) — seed well | says nothing about metals (trivially metallic) or about inert surfaces |
| **C. Structural templating** (lattice/texture match) | ZnO > TiO₂ at similar EA; low-resistivity CSL boundaries | molecular/amorphous seeds |
| **D. Alloying / adatom immobilisation** (dopant or reactive metal pins Ag) | Al-, Cu-, O-, N-doped Ag; Cu, Ge seeds | needs a metal or reactive gas in the electrode |
| **E. Condensation gate** (any reactive moiety at all) | fluorinated/perfluoro surfaces | everything that binds |
| F. Diffusion barrier E_d, Venables | qualitative only | absolute numbers unusable (±0.1 eV → 1.7× d_c) |

Read across: **A is the descriptor with the broadest support** — it orders
classes correctly across metals, oxides, polymers and fluorinated surfaces, and
the literature's own explanation (seed surface energy) is its macroscopic form.
**B is narrow but decisive where A fails**: the only cases where a stronger
binder makes a worse film are deep-acceptor oxides, and B is what distinguishes
MoO₃ from CuO, ZnO and HATCN. Neither A nor B alone covers the literature; the
two together do, except for the ZnO/TiO₂ structural effect (C).

An important nuance for B: an electron-affinity proxy is not enough. HATCN's
LUMO (EA ~5.7 eV) is nearly as deep as MoO₃'s conduction band, yet its charge is
*bounded* — a molecule accepts about one electron and then charges up, whereas a
reducible oxide offers a band of Mo⁶⁺/Mo⁵⁺ sites at high density. Capacity, not
depth, is the quantity; it has to be computed (cluster or slab charge growth),
not read off a band diagram.

## 3. Revised priority for screening (thermally evaporated top electrode)

**Gates (pass/fail, applied first):**

1. **Process.** Vacuum-evaporable onto a finished organic stack at room
   temperature, no solvent, no sputter damage. Removes PEI, PEDOT:PSS, SAMs,
   sputtered ZnO/CuO for this application (they stay valid as mechanistic
   evidence).
2. **Optical cost at working thickness.** The paper's thesis is electrode
   absorption; a seed that absorbs (Ge, Au, Ni, Cr; Cu after oxidation) must pay
   for itself in thinner Ag. Sub-nm seeds (Cu, CuO 0.25–1 nm) may pass.
3. **Condensation (descriptor E).** A reactive moiety must exist — computed E_b
   above roughly 0.5 eV at the deepest site. Removes inert aromatics and
   fluorinated layers (benzene-like surfaces condense but dewet).
4. **Electron capacity (descriptor B).** Veto if total Ag charge grows with
   cluster size (MoO₃-like). Pass if bounded (≤ ~0.6–0.7 e through Ag₄) or if
   the seed donates (amines, Cs₂CO₃, Ca).

**Ranking among survivors:**

5. **E_b at the deepest site (descriptor A)** — the literature's best-supported
   correlate. Higher = denser nucleation = earlier closure. Multi-site value,
   never single-site.

**Tie-breakers:**

6. Templating / seed smoothness (C), measured RMS of the bare seed.
7. E_d, only as an ordering sanity check.

**Compared with the previous ordering.** The 2D criterion (E_b × metallicity)
stands, but its two axes play different roles: metallicity is a **veto**, not a
weight, and E_b is the **ranking** variable. Treating metallicity as a ranking
weight would put benzene at the top — the literature says inert hydrocarbons
are among the worst. Treating E_b alone as the ranking puts MoO₃ above HATCN —
the literature and this lab both say no.

## 4. Where this leaves the candidates

| candidate | gate 1 process | gate 2 optics | gate 3 condensation | gate 4 capacity | E_b (multi-site) |
|---|---|---|---|---|---|
| HATCN | ✓ | ✓ (UV-blue edge) | ✓ | bounded | 1.63 (1.50 TZVP+CP) |
| F4TCNQ | ✓ | ✓ | ✓ | bounded | 1.21 |
| Bphen | ✓ | ✓ | ✓ | bounded (xTB) | 0.62 |
| Cs₂CO₃ | ✓ | ✓ | ✓ | donor (xTB) | 0.70 (single site, confirmed) |
| Mo3O9 / MoOx | ✓ | ~ | ✓ | **extensive → veto** | 2.16 (1.94 TZVP+CP) |
| Al4O6 / AlOx | ✓ | ✓ | ✓ | partly retained (xTB +0.27) | 0.70 (corrected) |
| Cu4I4 / CuI | ✓ | ✓ | ✓ | bounded (xTB) | 0.74 (corrected) |
| benzene-like hosts | ✓ | ✓ | ✗ (0.21) | — | — |
| Cu seed 1 nm | ✓ | ~ (oxidises) | ✓ | metallic | not computed |
| Al/Cu/Yb co-doped Ag | ✓ | ~ | ✓ | metallic | not computed |

HATCN remains the top evaporable transparent candidate under the revised
priority. **The competitors the paper must address are not other organics but
metal co-doping of Ag (Al, Cu, Yb)** — process-compatible, closing at 6–7 nm —
whose price is absorption from the dopant. That is the comparison the paper's
own thesis (lower electrode absorption) is built to win or lose.

## 5. Predictions that would test descriptor B (no new experiment needed to state them)

- **V₂O₅ and WO₃** (deep acceptors, layered/reducible like MoO₃) should give the
  MoO₃ morphology — early percolation, persistent porosity. No literature found
  that tests either as an ultrathin-Ag seed.
- **CuO and ZnO** should return charge from small Ag clusters (bounded). The
  literature's excellent CuO result is consistent; not yet computed here.
- **A charge-capacity calculation on ZnO(10-10), TiO₂(110), CuO and V₂O₅ slabs**
  at GFN2 (same protocol as the MoO₃ slab, script 137) would show whether the
  veto separates the oxides the literature already ranks — ZnO and CuO good,
  MoO₃ porous — before any new film is grown. That is the cheapest decisive test
  of the revised priority.

## Sources

- Logeeswaran et al., Nano Lett. 9, 178 (2009) — Ge nucleation layer.
- Formica et al., ACS AMI 5, 3048 (2013) — 1 nm Cu seed, Ag percolation 3 nm.
- Zhang, …, Guo, Adv. Mater. (2014) — Al-doped Ag.
- Ji, Liu, Zhang, Guo, Nat. Commun. 11, 3367 (2020) — Cu-doped Ag, 6.5 nm.
- Huang et al., Sol. Energy Mater. Sol. Cells (2018) — Cu-doped Ag percolation 6 nm.
- Martínez-Cercós (ICFO thesis); ACS AMI (2021) "Ultrathin metals on a transparent seed" — CuO seed, ~1 nm.
- Nanomaterials 8, 473 (2018), PMC6071051 — 1 nm MoO₃, 7 nm Ag porous.
- Dannenberg, Stach et al., OSTI 821432 — Ag on a-TiO₂ vs ZnO.
- Cueva & Carretero (Zaragoza) — Ag(O₂)/Ag(N₂) seeds.
- Frontiers in Materials 6, 18 (2019) — Au/PEI seed, 8 nm Ag.
- ACS Appl. Energy Mater. 7, 7140 (2024), PMC11412282 — condensation-coefficient review.
- Nat. Commun. 10, 866 (2019), PMC6382909 — Ag⁺–Bphen coordination n-doping.
- Park & Suh, Opt. Express 26, 4979 (2018) — HATCN/Ag 15–25 nm.

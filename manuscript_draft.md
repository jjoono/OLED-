# Limits of external light outcoupling in thin-film light-emitting devices: why shape is an exhausted degree of freedom

**[제목 후보]**
- *Limits of external light outcoupling in thin-film LEDs: shape is an exhausted degree of freedom*
- *Étendue bounds on directional outcoupling: why inverse design saturates at a hemisphere*
- *Efficiency and directionality cannot be increased together by shape design*

**Authors** — [TBD]

---

## Abstract

Microlens arrays (MLAs) and scattering layers are the standard external structures used to improve light
extraction from thin-film light-emitting devices, and freeform design and numerical inverse design are
increasingly applied to them in pursuit of tailored angular emission. Here we establish the boundary of
what such structures can achieve. Any passive external layer is completely described by its bidirectional
scattering distribution function, and energy conservation together with reciprocity then force the emitted
radiance to satisfy *L*<sub>out</sub> ≤ *L*<sub>in</sub>/*n*<sub>s</sub>², **irrespective of internal
multiple reflections, total-internal-reflection cascades or volume scattering**. The power deliverable
into a prescribed solid angle in a single pass is therefore bounded by the product of that solid angle's
étendue and the maximum attainable radiance; applying this one relation with different target solid
angles yields limits for a range of figures of merit within a single framework. We show that a **bare flat
interface already saturates this single-pass bound**, so no microlens or freeform profile can exceed a
flat surface on a single pass, and that efficiency and angular selectivity can only be **exchanged, not
simultaneously increased**, by shape design. Angular compression is not forbidden in general — it is paid
for in area — but an array coextensive with the emitter has an output-to-source area ratio of unity, which
removes that currency; this simultaneously explains why an index-matched hemispherical macroextractor
works and why it cannot be tiled. Three routes lie beyond the single-pass bound: radiance build-up by
recycling, microcavity engineering of the source distribution, and aperture expansion. Optimizing over all
physically admissible layers, the recycling bound reaches 62% of generated power in a 40–60° band at 10%
loss per round trip, whereas refractive structures are found to reach only 29%; we report this factor-two
gap explicitly and identify the reason refractive structures cannot close it as an open question.
Combined with a zero-mean-slope identity that forbids net azimuthal steering in any periodic tiling, these
results define a map of achievable and forbidden angular targets that depends only on substrate index and
geometry, and therefore applies across OLED, PeLED, QLED and micro-LED platforms.

> **[한글 주석]** 초록에서 (i) BSDF 유도로 일반성을 확보, (ii) "패스당" 한정 명시,
> (iii) 62% vs 29% 간극을 결과로 승격, (iv) 집광=면적 대가로 재서술했습니다.

---

## 1. Introduction

Thin-film light-emitting devices such as OLEDs, PeLEDs and QLEDs offer large-area surface emission and
compatibility with flexible and conformable form factors, and are consequently regarded as key light
sources for both lighting and display applications. Their external quantum efficiency, however, remains
well below that of inorganic LEDs, which typically reach 70–80% extraction, and this gap continues to
limit high-luminance operation and device lifetime [ref].

Internal quantum efficiency is no longer the bottleneck. With phosphorescent and TADF emitters, the
internal radiative efficiency approaches unity [ref]. What remains limiting is light extraction: in planar
stacks, waveguided modes in the organic layers, total internal reflection at the substrate/air interface,
and coupling to surface plasmon polaritons at the metal electrode together trap a large fraction of the
generated photons, leaving only about 20–25% extracted [ref].

Among the various strategies developed to address this, **external outcoupling structures** such as
microlens arrays and scattering layers occupy a distinctive position. Because they are formed outside the
device, they do not perturb its electrical operation. They are also largely agnostic to the emitter
material and to the fabrication route, which makes them broadly applicable across OLED, PeLED and QLED
platforms, and they provide substantial efficiency gains, which has made them attractive for industrial
adoption [ref].

The limits of extraction **efficiency** are already well understood. Yablonovitch's statistical ray optics
established that the path-length enhancement available in a high-index medium with a randomizing texture
is bounded by 4*n*² [Yablonovitch 1982], and the decomposition of optical loss channels in planar OLEDs
has been analysed in detail [Brütting 2013]. Brütting and co-workers further showed that horizontally
oriented phosphorescent emitters with near-unity internal radiative efficiency could reach external
quantum efficiencies of up to 70% **if all of the substrate emission could be extracted** — a condition
realized in practice with an index-matched hemispherical macroextractor, which is not applicable to
large-area devices [Brütting 2013].

> **[한글 주석]** Brütting 인용에 "반구 macroextractor 전제" 단서를 유지. 2.2절에서 이 예외가
> 에텐듀 식으로부터 유도되므로, 서론의 이 문장과 이론이 연결됩니다.

What recent applications increasingly demand, however, is not the total amount of light but the **shape of
its angular distribution**. Augmented-reality and waveguide displays require coupling into a specific
angular range [ref]; privacy and dual-view displays require azimuthal asymmetry [ref]; inspection and
machine-vision illumination requires uniformity over a target region [ref]. Accordingly, a large body of
work has sought to improve figures of merit beyond EQE alone, including on-axis luminance and angular
colour stability [ref].

Shaping an angular distribution is itself a mature subject. Freeform design methods based on the
Monge–Ampère equation, ray mapping and supporting quadrics have been developed over roughly two decades
and are now well established in illumination optics [ref]. More recently, these design and inverse-design
techniques have been extended to the outcoupling structures of thin-film emitters, with numerical
optimization used to target both angular distribution and efficiency [ref]. Crucially, however, these
methods were largely developed under a **point-source or near-zero-étendue assumption**, whereas a
thin-film emitter presents an **extended source** comparable in size to the lenslets, and its behaviour is
dominated by **total internal reflection and recycling** inside the substrate. The operating conditions
are therefore fundamentally different.

Under these conditions it has not been established how far an external outcoupling structure can go with
respect to **directional figures of merit**: which angular distributions are attainable, which are
forbidden, what sets the boundary, and whether adding design freedom can move it. This work addresses that
question.

---

## 2. Theory

### 2.1 The master relation, and why internal complexity does not matter

**Radiance bound from the BSDF.** Whatever its internal construction, a passive external layer is
completely characterized by its bidirectional scattering distribution function *f*(θ<sub>i</sub> →
θ<sub>o</sub>), which by definition already contains every internal multiple reflection, total-internal-
reflection cascade and volume-scattering event. Any physically admissible BSDF satisfies energy
conservation,

$$\int f(\theta_o\!\to\!\theta_i)\cos\theta_i \, d\Omega_i \le 1 ,$$

and reciprocity across an index boundary,

$$n_i^{2} f(\theta_i\!\to\!\theta_o) = n_o^{2} f(\theta_o\!\to\!\theta_i) .$$

For Lambertian incidence of radiance *L* from the substrate, the emergent radiance is

$$L_{out}(\theta_o) = \int f(\theta_i\!\to\!\theta_o)\,L\cos\theta_i \, d\Omega_i
= L\,\frac{n_o^{2}}{n_i^{2}}\int f(\theta_o\!\to\!\theta_i)\cos\theta_i \, d\Omega_i
\;\le\; \frac{L}{n_s^{2}} .$$

The radiance bound is thus a direct consequence of **energy conservation and reciprocity alone**. No
assumption about the structure's geometry, feature size or internal ray history is required, so every
microlens array, freeform profile and scattering layer is automatically subject to it. This is what allows
the analysis below to cover a structural class that cannot be enumerated.

> **[한글 주석 — 중요]** 이 유도가 "모든 MLA·산란층을 어떻게 다 확인했나"라는 질문을 원천 차단합니다.
> 확인할 필요 없이 BSDF 제약이 보증합니다. 앞서 우려하셨던 지점의 답입니다.

**Consistency across descriptions.** The same bound follows from three independent routes, which differ
only in how tightly they can be approached:

| Description | Origin of the bound | Attainability |
|---|---|---|
| Single refraction | Snell's law conserves étendue exactly | **Saturated** by a flat interface |
| Multiple internal reflections | Liouville's theorem: phase-space volume is conserved | Generally **below** the bound — mixing raises étendue |
| BSDF | Energy conservation + reciprocity | Depends on the particular BSDF |

Multiple internal reflections therefore never help: each scattering event preserves or increases étendue,
so internal complexity moves a structure away from, not towards, the bound.

**The master relation.** Combining the radiance bound with the definition of the power delivered into a
target solid angle Ω from emitting area *A*, and normalizing by the total power crossing into the
substrate hemisphere *P*<sub>sub</sub> = *L A* π, gives

$$\boxed{\;\frac{P_\Omega}{P_{\mathrm{sub}}}\;\le\;\frac{1}{\pi n_s^{2}}\int_\Omega \cos\theta\, d\Omega\;}
\qquad \text{(single pass)}$$

Different figures of merit correspond only to different choices of Ω. For an off-axis band
θ ∈ [θ<sub>lo</sub>, θ<sub>hi</sub>] the integral gives π(sin²θ<sub>hi</sub> − sin²θ<sub>lo</sub>), so

$$\frac{P_{\mathrm{band}}}{P_{\mathrm{sub}}}\;\le\;\frac{\sin^{2}\theta_{hi}-\sin^{2}\theta_{lo}}{n_s^{2}} ,$$

which for the 40–60° band with *n*<sub>s</sub> = 1.51 evaluates to **14.8%** per pass.

### 2.2 A flat interface saturates the single-pass bound

A planar dielectric interface maps substrate angles onto air angles one to one through Snell's law without
scrambling étendue, and therefore saturates the master relation: the fraction of escaping light in the
40–60° band is exactly sin²60° − sin²40° = 0.337, which combined with the escape probability
1/*n*<sub>s</sub>² gives 14.8% of the substrate power. An **ideal angular filter** transmitting only the
target band reaches the same single-pass value, transmitting 14.8% with unit selectivity; both extremes
are constrained by the same étendue–radiance product.

**Corollary.** No microlens array or freeform profile can exceed a flat interface in single-pass band
power. The saturation of large-scale inverse design at a near-hemispherical profile is therefore a
consequence of a conservation law, not a failure of the optimizer.

**Angular compression costs area, and tiling removes the currency.** Angular compression is not forbidden
in general. Étendue conservation, *A*<sub>src</sub>Ω<sub>src</sub> = *A*<sub>out</sub>Ω<sub>out</sub>,
states that the output solid angle can be reduced only in proportion to an increase in output area. A
single lens collimates precisely because its aperture greatly exceeds the source it covers. What matters
is the **ratio** *A*<sub>out</sub>/*A*<sub>src</sub>, not the absolute lenslet size: enlarging the lenslets
enlarges the source patch beneath them in the same proportion. An array that tiles a large-area emitter is
coextensive with it, so the ratio is fixed at unity and the currency required for compression is simply
absent. This yields the exception noted in the Introduction as a consequence of the same relation: an
index-matched hemispherical macroextractor achieves *A*<sub>out</sub>/*A*<sub>src</sub> ≫ 1 over a small
emitter, which is why it extracts the substrate modes and equally why it cannot be tiled over a large
area.

A secondary, size-dependent penalty applies to small lenslets. Light spreads laterally in the substrate
before reaching the array, so each lenslet of radius *a* draws from an effective radius
*r*<sub>eff</sub> = *a* + *d*<sub>sub</sub> tan θ<sub>c</sub>. For *a* ≪ *d*<sub>sub</sub> tan θ<sub>c</sub>
neighbouring lenslets share light heavily and each sees a nearly Lambertian angular input, removing any
local compression as well.

### 2.3 Efficiency and directionality can only be exchanged

Writing the band power as the product of total extraction efficiency and angular selectivity,
*P*<sub>band</sub> = η<sub>ext</sub> × *S*, the master relation bounds the **product**:

$$\eta_{\mathrm{ext}}\times S \;\le\; \text{const}.$$

Shape design can move a structure along this constraint — trading total extraction for angular
concentration or the reverse — but cannot move it outward. We therefore predict that designs obtained from
independent optimization runs, whatever their individual shapes, populate a common iso-product contour in
the (η<sub>ext</sub>, *S*) plane.

> **[한글 주석 — 중요]** 이 예측은 **아직 검증 전**입니다. BO 결과를 (EQE_total, 선택성) 평면에 찍어
> 등곱 곡선 위에 놓이는지 확인해야 합니다.
> - **성립** → 이 절이 논문의 간판 결과, 제목 반영
> - **불성립** → "곱이 유계된다"는 진술만 남기고 iso-product 예측은 삭제
>
> **투고 전 반드시 검증할 것.**

### 2.4 Azimuthal steering is forbidden in periodic tilings

For a surface *z*(*x*, *y*) that tiles periodically with height continuous across cell boundaries, the
mean surface slope over one cell vanishes,

$$\oint \nabla z \, dA = 0 ,$$

because contributions from opposite cell edges cancel. In the thin-element limit the prismatic deviation
imparted to a ray is proportional to the local slope, so the mean angular deviation is zero: a uniform
periodic array cannot impart net beam steering. The azimuthal fraction reaching a window of width Δφ is
therefore limited to approximately Δφ/360°, i.e. 11% for a 40° window.

This statement applies to **refractive, height-continuous periodic arrays**. Metasurfaces impart a phase
gradient and are not subject to this identity; they lie outside the class considered here [ref].

> **[한글 주석]** 이 한정 문구가 없으면 메타표면 논문 하나로 반박당합니다. 반드시 유지.

### 2.5 Beyond the single-pass bound: three routes

The master relation is stated per pass, for a passive layer acting on a given substrate radiance. Each of
its premises corresponds to a distinct route past it.

**(i) Radiance build-up by recycling.** Light that fails to escape can be returned and given further
attempts, raising the steady-state substrate radiance above its single-pass value and pushing more power
through the same étendue. Build-up requires that out-of-band light be **selectively rejected**, and it is
paid for in loss, since each additional round trip incurs absorption at the reflective electrode.

**(ii) Microcavity engineering of the source distribution.** The bound above assumes Lambertian substrate
radiance. A microcavity redistributes the emission itself, rendering *I*<sub>sub</sub>(θ) non-Lambertian.
Because this modifies the source rather than acting as a passive layer upon it, it is not constrained by
the radiance bound and constitutes a lever independent of the external structure.

**(iii) Aperture expansion.** As shown in Section 2.2, raising *A*<sub>out</sub>/*A*<sub>src</sub> above
unity permits angular compression. This is the operating principle of the hemispherical macroextractor,
and it is unavailable to any structure that tiles the emitter.

Three further routes lie outside the scope of this work, since each violates a premise of the radiance
bound: **non-reciprocal** elements such as magneto-optical media, **frequency conversion** as used in
luminescent concentrators, and **active or gain** media.

**The gap between the recycling bound and refractive structures.** Optimizing over all physically
admissible BSDFs subject only to energy conservation and reciprocity, the optimum is an angular filter,
and the corresponding recycling bound reaches 62% of generated power in the 40–60° band at 10% loss per
round trip. Refractive structures — flat interfaces, microlens arrays, freeform profiles and scattering
layers — reach only about 29% under the same conditions. **This factor-two gap is a result of this work,
not an artefact**: the 29% figure is an achievable frontier for non-selective layers, whereas the 62%
figure is a bound over the full admissible class.

Whether refractive structures can generate the angular selectivity required to close that gap is, in our
assessment, **not settled**. The reflectance of a smooth dielectric surface is governed by Fresnel
coefficients and total internal reflection, which vary smoothly and monotonically with local incidence
angle, so placing a pass-band at an arbitrary angle with rejection on both sides appears to require
interference. Consistent with this, we find numerically that the selectivity of a flat interface is
completely insensitive to the asymmetry parameter of an adjacent scattering layer (Section 4.4). We
nevertheless present this as an **open question** rather than a proof, since energy conservation and
reciprocity alone do not exclude a selective refractive BSDF.

> **[한글 주석]** 이 "열린 문제" 선언이 논문의 정직성이자 방어막입니다. 그리고 62% vs 29% 간극을
> 숨기지 않고 명시적 결과로 올린 것이 v7의 핵심 변경입니다.

---

## 3. Methods

### 3.1 Trans-scale optical model
Dipole emission inside the microcavity is computed with the classical-power-spectrum (CPS) formalism,
yielding the substrate-side angular and spectral intensity *I*<sub>sub</sub>(θ, λ), which then serves as
the source for macroscopic ray tracing through the external structure [ref: 자기 ACS Photonics 2023,
Nat Commun 2025]. The framework has previously been validated experimentally for OLED microlens
structures, where simulated and measured EQE and angular profiles agreed within experimental uncertainty
[refs].

### 3.2 Inverse design
Freeform lens profiles are parameterized by [스플라인 제어점 / RBF …] and optimized with a
surrogate-based global optimizer (RBF surrogate with multi-start and pattern-search polishing) under
geometric feasibility constraints. [세부: 변수 수, 평가 예산, multi-start 수 기입]

### 3.3 Recycling model
The external layer is represented by its angular transmission and reflection, and the round trip through
the substrate to the reflective electrode by an absorption probability *a*. Summing the resulting series
gives the total power delivered to the target band. For an angularly selective layer this is evaluated
using transfer-matrix calculations of a dielectric multilayer, which serves here as a **representative
realization of the selective limit** rather than as a proposed device.

> **[한글 주석]** 3.3절이 "산란층+DBR 구조를 제안한다"로 읽히지 않도록 명시했습니다.

---

## 4. Results

### 4.1 The inverse-design frontier saturates at the single-pass bound
[Fig 2] Across [N] independently optimized freeform profiles, the band power does not exceed the value
attained by a flat interface, and the best designs converge to near-hemispherical profiles. We refer to
this saturation curve as an **achievable frontier**; it is an outcome of finite search and is conceptually
distinct from the bounds of Section 2, which follow from conservation laws.

### 4.2 Efficiency–selectivity trade-off
[Fig 2b] Designs populate a common contour in the (η<sub>ext</sub>, *S*) plane, consistent with the
bounded product of Section 2.3. **[검증 후 확정]**

### 4.3 Map of figures of merit
[Fig 3] Applying the master relation with different target solid angles gives bounds for on-axis
enhancement, off-axis band concentration and flat-top uniformity, while the zero-mean-slope identity
renders net azimuthal steering and dark-split dual-view distributions unattainable. We further find that
angular mixing, which dilutes viewing-angle colour variation, directly opposes angular concentration, so
colour stability and directionality constitute a trade-off rather than independent targets. [단색 기준]

### 4.4 Recycling: the bound, the frontier, and the gap
[Fig 4a,b] At *a* = 0.1 per round trip, the bound over all admissible layers — realized by an ideal
angular filter — delivers 62.1% of generated power into the 40–60° band, whereas non-selective refractive
layers reach 29.1%, a factor of 2.1 below. The selectivity of the latter is found to be completely
independent of the scattering asymmetry parameter over *g* ∈ [−0.5, +0.5], indicating that angular
redistribution by non-selective scattering does not contribute to build-up. For a realizable dielectric
multilayer the attainable value degrades with source bandwidth, from 48.2% for a monochromatic source to
31.0% over a 100 nm band, because the stopband edge shifts with wavelength. This reproduces the behaviour
reported for angle-filtered OLEDs, where angular selection is accompanied by increased absorption arising
from repeated interaction with absorbing layers [ref].

### 4.5 Consistency and generality
[Fig 4c,d] The analytic recycling model and full three-dimensional ray tracing agree to within [X]% across
[N] conditions. Predicted and simulated directional EQE for the studied geometry are 3.2% and 2.93%
respectively. Because the bounds involve only substrate index and geometry, they transfer directly to
PeLED, QLED and micro-LED platforms; sweeping *n*<sub>s</sub> [Fig 4c] shows the expected scaling.

---

## 5. Discussion

The central practical implication is that **shape is an exhausted degree of freedom** for external
outcoupling. Further refinement of lenslet profiles, whether by analytic freeform methods or by
large-scale numerical inverse design, redistributes a bounded product rather than increasing it. Effort is
better directed at the levers identified in Section 2.5 — radiance build-up through angular rejection,
microcavity engineering of the source distribution, and, where planarity can be sacrificed, aperture
expansion — and at reducing the loss that caps the first of these.

The factor-two gap between the recycling bound and what refractive structures achieve deserves emphasis,
because it is the quantity a designer can still act upon. Our results show that the gap cannot be closed
by shape optimization, and identify angular selectivity as the property that would close it, but we do not
establish whether a refractive structure can possess that property.

Our claims are restricted to refractive, height-continuous periodic arrays. Metasurfaces and resonant
diffractive structures impart phase gradients or support guided-mode resonances and are not subject to the
zero-mean-slope identity; non-periodic, per-pixel structures likewise fall outside the class, at the cost
of alignment and patterning requirements that periodic arrays avoid.

By reciprocity, the geometric half of these results transfers to the absorption problem: the set of
attainable angular distributions for emission through a passive layer is isomorphic to its angular
acceptance when the same layer is used as an absorber. The efficiency ceiling does not transfer, because
the useful channel differs — escape into air for emission, absorption in the active layer for detection.

Finally we note the limitations of the analysis: the non-selectivity of refractive structures is argued
rather than proven; the scattering layer is treated in a two-stream approximation; and the treatment is
in the geometrical-optics limit.

---

## 6. Conclusion

We have shown that energy conservation and reciprocity alone bound the radiance emerging from any passive
external outcoupling layer, and hence the power it can deliver into a prescribed angular target in a
single pass, independently of the layer's internal complexity. A flat interface already saturates this
bound, which explains the saturation of inverse design at near-hemispherical profiles as a consequence of
conservation rather than of optimization, and implies that efficiency and directionality can only be
exchanged by shape design. Angular compression requires an output-to-source area ratio above unity, a
currency that tiling removes — the same relation that accounts for the success of the hemispherical
macroextractor and for its non-scalability. Beyond the single-pass bound lie radiance build-up by
recycling, microcavity engineering of the source, and aperture expansion; of these, recycling admits a
bound a factor of two above what refractive structures achieve, and closing that gap requires angular
selectivity whose availability to refractive structures remains open.

---

## Figures

| Fig | Content | Status |
|---|---|---|
| 1 | Platform, étendue picture, **three levers** (shape exhausted / recycling / cavity) | 작성 필요 |
| 2 | Single-pass bound, flat saturation, iso-product plane, frontier | **iso-product 검증 필요** |
| 3 | FoM map (Ω별 상한, 금지 영역, 색–방향성 트레이드오프) | 계산 일부 필요 |
| 4 | Recycling **bound vs frontier gap**, bandwidth, generality, cross-validation | (a)(b) 계산 완료 |

## References — 우선 확보 목록
- Yablonovitch, *J. Opt. Soc. Am.* **72**, 899 (1982)
- Brütting et al., *Phys. Status Solidi A* **210**, 44 (2013) — 70% / macroextractor 단서 포함
- Winston, Jiang, Ricketts, *Adv. Opt. Photon.* **10**, 484 (2018)
- 자기 그룹: ACS Photonics **10**, 1775 (2023); Nat. Commun. (2025)
- Freeform illumination 리뷰 / Monge–Ampère
- 각도필터 OLED (Adv. Photonics Res. 2023), angle-selective film (2022)
- 메타표면 방향성 microLED (Comm. Eng. 2025)
- BSDF 상반성 / radiative transfer 표준 문헌

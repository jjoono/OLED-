# Limits of external light outcoupling in thin-film light-emitting devices: why shape is an exhausted degree of freedom

**[제목 후보]**
- *Limits of external light outcoupling in thin-film LEDs: shape is an exhausted degree of freedom*
- *Étendue bounds on directional outcoupling: why inverse design saturates at a hemisphere*
- *Efficiency and directionality cannot be increased together in external outcoupling structures*

**Authors** — [TBD]

---

## Abstract

Microlens arrays (MLAs) and scattering layers are the standard external structures used to improve
light extraction from thin-film light-emitting devices, and freeform design and numerical inverse
design are increasingly applied to them in pursuit of tailored angular emission. Here we establish the
boundary of what such structures can achieve. Because a passive layer cannot increase basic radiance,
the power that can be delivered into a prescribed solid angle in a single pass is bounded by the product
of that solid angle's étendue and the maximum attainable radiance. Applying this single relation with
different target solid angles yields limits for a range of figures of merit within one framework. We
show that a **bare flat interface already saturates this bound** — it maps angles one to one without
scrambling étendue — so no microlens or freeform profile can exceed a flat surface on a single pass.
This explains, as a consequence of a conservation law rather than as an optimization outcome, why
large-scale inverse design saturates at a near-hemispherical profile. It further implies that total
extraction efficiency and angular selectivity can only be **exchanged, not simultaneously increased**,
by shape design. The only remaining route past the bound is radiance build-up through recycling, which
requires angular rejection and whose gain is capped by mirror loss: at 10% loss per bounce, a
non-selective layer reaches 29% of generated power in a 40–60° band whereas an ideal angular filter
reaches 62%. Combined with a zero-mean-slope identity that forbids net azimuthal steering in any
periodic tiling, these results define a map of achievable and forbidden angular targets. The bounds
depend only on substrate index and geometry, and therefore apply across OLED, PeLED, QLED and
micro-LED platforms.

> **[한글 주석]** 초록에서 "우리가 새 구조를 제안한다"는 뉘앙스를 철저히 배제했습니다.
> 모든 상한은 보존법칙 유도이고, 시뮬레이션은 frontier 확인용으로만 등장합니다.

---

## 1. Introduction

Thin-film light-emitting devices such as OLEDs, PeLEDs and QLEDs offer large-area surface emission and
compatibility with flexible and conformable form factors, and are consequently regarded as key light
sources for both lighting and display applications. Their external quantum efficiency, however, remains
well below that of inorganic LEDs, which typically reach 70–80% extraction, and this gap continues to
limit high-luminance operation and device lifetime [ref].

Internal quantum efficiency is no longer the bottleneck. With phosphorescent and TADF emitters, the
internal radiative efficiency approaches unity [ref]. What remains limiting is light extraction: in
planar stacks, waveguided modes in the organic layers, total internal reflection at the substrate/air
interface, and coupling to surface plasmon polaritons at the metal electrode together trap a large
fraction of the generated photons, leaving only about 20–25% extracted [ref].

Among the various strategies developed to address this, **external outcoupling structures** such as
microlens arrays and scattering layers occupy a distinctive position. Because they are formed outside
the device, they do not perturb the electrical operation of the stack. They are also largely agnostic
to the emitter material and to the fabrication route, which makes them broadly applicable across OLED,
PeLED and QLED platforms, and they provide substantial efficiency gains, which has made them attractive
for industrial adoption [ref].

The limits of extraction **efficiency** are already well understood. Yablonovitch's statistical ray
optics established that the path-length enhancement available in a high-index medium with a randomizing
texture is bounded by 4n² [Yablonovitch 1982], and the decomposition of optical loss channels in planar
OLEDs has been analysed in detail [Brütting 2013]. Brütting and co-workers further showed that
horizontally oriented phosphorescent emitters with near-unity internal radiative efficiency could reach
external quantum efficiencies of up to 70% **if all of the substrate emission could be extracted** — a
condition realized in practice with an index-matched hemispherical macroextractor, which is not
applicable to large-area devices [Brütting 2013].

> **[한글 주석]** ← Brütting 인용에 "반구 macroextractor 전제"라는 단서를 반드시 붙였습니다.
> 이 단서가 곧 본 논문의 동기로 이어지며, 동시에 별도 준비 중인 논문과의 경계도 지킵니다.

What recent applications increasingly demand, however, is not the total amount of light but the **shape
of its angular distribution**. Augmented-reality and waveguide displays require coupling into a specific
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
forbidden, what sets the boundary, and whether adding design freedom can move it. This work addresses
that question.

We show that a single relation — the étendue of the target solid angle multiplied by the maximum
attainable radiance — bounds the power deliverable into any prescribed angular target, and that a flat
interface already saturates it. Consequently, shape is an exhausted degree of freedom: microlens and
freeform profiles redistribute the bounded product of efficiency and selectivity rather than increasing
it. We quantify the only remaining lever, radiance build-up by recycling, and its cost in loss and
bandwidth, and we combine these with a geometric identity that forbids net azimuthal steering in
periodic tilings. The resulting bounds depend only on substrate index and geometry and thus transfer
across emitter platforms.

---

## 2. Theory

### 2.1 The master relation

Consider a passive external layer on a substrate of refractive index *n*<sub>s</sub>, with a Lambertian
radiance *L* incident from the substrate side. A passive, linear, reciprocal optical element cannot
increase the basic radiance *L*/*n*², so the radiance emerging into air satisfies
*L*<sub>air</sub> ≤ *L*/*n*<sub>s</sub>².

The power delivered into a prescribed target solid angle Ω from emitting area *A* is therefore bounded by

$$P_\Omega \;\le\; \frac{L}{n_s^2}\,A\int_\Omega \cos\theta\, d\Omega ,$$

while the total power crossing into the substrate hemisphere is *P*<sub>sub</sub> = *L A* π. Their ratio
gives the **master relation**

$$\boxed{\;\frac{P_\Omega}{P_{\mathrm{sub}}}\;\le\;\frac{1}{\pi n_s^{2}}\int_\Omega \cos\theta\, d\Omega\;}
\qquad \text{(single pass)}$$

This expression uses only passivity, linearity and reciprocity of the layer. It is therefore
**independent of the internal complexity of the structure** — multiple internal reflections, total
internal reflection cascades and volume scattering are all already accounted for, because the constraint
is applied at the level of the layer as a whole.

Its practical value is that **different figures of merit correspond only to different choices of Ω**. For
an off-axis band θ ∈ [θ<sub>lo</sub>, θ<sub>hi</sub>], the integral evaluates to
π(sin²θ<sub>hi</sub> − sin²θ<sub>lo</sub>), giving

$$\frac{P_{\mathrm{band}}}{P_{\mathrm{sub}}}\;\le\;\frac{\sin^{2}\theta_{hi}-\sin^{2}\theta_{lo}}{n_s^{2}} .$$

For the 40–60° band with *n*<sub>s</sub> = 1.51 this yields **14.8%** per pass.

> **[한글 주석]** 이 절이 논문의 심장입니다. "FoM마다 새 이론"이 아니라 "Ω만 바꾼다"는 점을
> 여기서 명시적으로 선언해야 3절의 지도가 설득력을 갖습니다.

### 2.2 A flat interface saturates the bound

A planar dielectric interface maps substrate angles onto air angles one to one through Snell's law
without scrambling étendue. It is therefore an étendue-preserving element, and it saturates the master
relation: of the light escaping a flat interface, the fraction falling in the 40–60° band is exactly
sin²60° − sin²40° = 0.337, and combined with the escape probability 1/*n*<sub>s</sub>² this gives 14.8%
of the substrate power — precisely the bound.

Notably, an **ideal angular filter** that transmits only the target band and rejects everything else
reaches the same value in a single pass: it transmits 14.8% of substrate power with 100% selectivity.
Both extremes are constrained by the same étendue–radiance product.

**Corollary.** No microlens array or freeform profile can exceed a flat interface in single-pass band
power. The saturation of large-scale inverse design at a near-hemispherical profile is therefore a
consequence of a conservation law, not a failure of the optimizer.

### 2.3 Efficiency and directionality can only be exchanged

Writing the band power as the product of total extraction efficiency and angular selectivity,
*P*<sub>band</sub> = η<sub>ext</sub> × *S*, the master relation bounds the **product**:

$$\eta_{\mathrm{ext}}\times S \;\le\; \text{const}.$$

Shape design can move a structure along this constraint — trading total extraction for angular
concentration or vice versa — but cannot move it outward. We therefore predict that designs obtained by
independent optimization runs, whatever their individual shapes, should populate a common iso-product
contour in the (η<sub>ext</sub>, *S*) plane.

> **[한글 주석 — 중요]** 이 예측은 **아직 검증 전**입니다. BO 결과들을 (EQE_total, 선택성)
> 평면에 찍어 등곱 곡선 위에 놓이는지 확인해야 합니다. 확인되면 이 절이 논문의 간판이 되고,
> 확인되지 않으면 "bounded product" 진술만 남기고 iso-product 예측은 빼야 합니다.
> **초안 제출 전 반드시 검증할 것.**

### 2.4 Azimuthal steering is forbidden in periodic tilings

For a surface *z*(*x*, *y*) that tiles periodically with height continuous across cell boundaries, the
mean surface slope over one cell vanishes,

$$\oint \nabla z \, dA = 0 ,$$

because contributions from opposite cell edges cancel. In the thin-element limit the prismatic deviation
imparted to a ray is proportional to the local slope, so the **mean angular deviation is zero**: a
uniform periodic array cannot impart net beam steering. The azimuthal fraction reaching a window of
width Δφ is therefore limited to approximately Δφ/360°, i.e. 11% for a 40° window.

We stress that this statement applies to **refractive, height-continuous periodic arrays**. Metasurfaces
impart a phase gradient and are not subject to this identity; they lie outside the class considered here
[ref].

> **[한글 주석]** 이 한정 문구가 없으면 메타표면 논문 하나로 반박당합니다. 반드시 유지.

### 2.5 The remaining lever: radiance build-up by recycling

The master relation applies per pass. Light that fails to escape can be returned, re-randomized and
given further attempts, and in steady state this raises the substrate radiance above its single-pass
value, allowing more power through the same étendue. Build-up is therefore the only mechanism by which
the single-pass bound can be exceeded.

Build-up requires that out-of-band light be **selectively rejected**; a layer that transmits broadly
returns little and accumulates little. It is also paid for in loss: each additional round trip incurs
absorption at the reflective electrode. Writing *a* for the absorption probability per round trip, the
attainable band power follows from summing the recycling series (Section 3.3).

Whether refractive structures can generate the required selectivity is, in our assessment, **not
settled**. The reflectance of a smooth dielectric surface is governed by Fresnel coefficients and total
internal reflection, which vary smoothly and monotonically with local incidence angle; placing a
pass-band at an arbitrary angle with rejection on both sides appears to require interference. Consistent
with this, we find numerically that the selectivity of a flat interface is completely insensitive to the
asymmetry parameter of an adjacent scattering layer (Section 4.3). We nevertheless present this as an
**open question** rather than a proof, since a general treatment including multiple internal reflections
is not available.

> **[한글 주석]** 이 "열린 문제" 선언이 논문의 정직성이자 방어막입니다. 증명인 척하면 잡힙니다.

---

## 3. Methods

### 3.1 Trans-scale optical model
Dipole emission inside the microcavity is computed with the classical-power-spectrum (CPS) formalism,
yielding the substrate-side angular and spectral intensity *I*<sub>sub</sub>(θ, λ), which is then used as
the source for macroscopic ray tracing through the external structure [ref: 자기 ACS Photonics 2023,
Nat Commun 2025 — 프레임워크의 실측 검증 근거로 인용]. The framework has previously been validated
experimentally for OLED microlens structures, where simulated and measured EQE and angular profiles
agreed within experimental uncertainty [refs].

### 3.2 Inverse design
Freeform lens profiles are parameterized by [spline control points / RBF …] and optimized with a
surrogate-based global optimizer (RBF surrogate with multi-start and pattern-search polishing) under
geometric feasibility constraints. [세부: 변수 수, 예산, multi-start 수 기입]

### 3.3 Recycling model
The external layer is represented by its angular transmission and reflection, and the round trip through
the substrate to the reflective electrode by an absorption probability *a*. Summing the resulting series
gives the total power delivered to the target band. For an angularly filtering layer this is evaluated
using transfer-matrix calculations of a dielectric multilayer.

> **[한글 주석]** 3.3에서 "우리가 산란층+DBR 구조를 제안한다"로 읽히지 않게, 어디까지나
> **모델 검증용 대표 구조**로 서술할 것.

---

## 4. Results

### 4.1 The inverse-design frontier saturates at the bound
[Fig 2] Across [N] independently optimized freeform profiles, the band power does not exceed the value
attained by a flat interface, and the best designs converge to near-hemispherical profiles. We refer to
this saturation curve as an **achievable frontier**; it is an outcome of finite search and is
conceptually distinct from the bounds of Section 2, which follow from conservation laws.

### 4.2 Efficiency–selectivity trade-off
[Fig 2b] Designs populate a common contour in the (η<sub>ext</sub>, *S*) plane, consistent with the
bounded product of Section 2.3. **[검증 후 확정]**

### 4.3 Map of figures of merit
[Fig 3] Applying the master relation with different target solid angles gives bounds for on-axis
enhancement, off-axis band concentration and flat-top uniformity, while the zero-mean-slope identity
renders net azimuthal steering and dark-split dual-view distributions unattainable. We further find that
angular mixing, which dilutes viewing-angle colour variation, directly opposes angular concentration, so
colour stability and directionality constitute a trade-off rather than independent targets. [단색 기준]

### 4.4 Recycling: loss and bandwidth
[Fig 4a,b] At *a* = 0.1 per round trip, a non-selective layer delivers 29.1% of generated power into the
40–60° band, whereas an ideal angular filter delivers 62.1%; the gain vanishes as *a* increases. For a
realizable dielectric multilayer, the attainable value degrades with source bandwidth, from 48.2% for a
monochromatic source to 31.0% over a 100 nm band, because the stopband edge shifts with wavelength. This
reproduces the behaviour reported for angle-filtered OLEDs, where angular selection is accompanied by
increased absorption from repeated interaction with absorbing layers [ref].

### 4.5 Consistency and generality
[Fig 4c,d] The analytic recycling model and full three-dimensional ray tracing agree to within [X]%
across [N] conditions. Predicted and simulated directional EQE for the studied geometry are 3.23% and
2.93% respectively. Because the bounds involve only substrate index and geometry, they transfer directly
to PeLED, QLED and micro-LED platforms; sweeping *n*<sub>s</sub> [Fig 4c] shows the expected scaling.

---

## 5. Discussion

The central practical implication is that **shape is an exhausted degree of freedom** for external
outcoupling. Further refinement of lenslet profiles, whether by analytic freeform methods or by
large-scale numerical inverse design, redistributes a bounded product rather than increasing it. Effort
is better directed at the one remaining lever identified here — radiance build-up through angular
rejection — and at reducing the loss that caps it.

Our claims are restricted to refractive, height-continuous periodic arrays. Metasurfaces and resonant
diffractive structures impart phase gradients or support guided-mode resonances and are not subject to
the zero-mean-slope identity; non-periodic, per-pixel structures likewise fall outside the class, at the
cost of alignment and patterning requirements that periodic arrays avoid.

By reciprocity, the geometric half of these results transfers to the absorption problem: the set of
attainable angular distributions for emission through a passive layer is isomorphic to its angular
acceptance when used as an absorber. The efficiency ceiling does not transfer, because the useful channel
differs — escape into air for emission, absorption in the active layer for detection.

Finally, we note the limitations of the present analysis. [열린 문제 재확인: 굴절형 구조의 비선별성에
대한 일반 증명 부재; 산란층의 2-state 근사; 기하광학 극한.]

---

## 6. Conclusion

We have derived a single relation bounding the power that an external outcoupling structure can deliver
into a prescribed angular target, shown that a flat interface already saturates it, and thereby explained
the saturation of inverse design at near-hemispherical profiles as a consequence of conservation rather
than of optimization. Applying the relation with different target solid angles yields a map of attainable
and forbidden angular distributions for thin-film emitters, and identifies radiance build-up by recycling
— not lens shape — as the only remaining route beyond the bound.

---

## Figures

| Fig | Content | Status |
|---|---|---|
| 1 | Platform, étendue picture, two levers | 작성 필요 |
| 2 | Per-pass bound, flat saturation, iso-product plane, frontier | **iso-product 검증 필요** |
| 3 | FoM map (Ω별 상한, 금지 영역, 색-방향성 트레이드오프) | 계산 일부 필요 |
| 4 | Recycling (loss, bandwidth), generality, cross-validation | (a)(b) 계산 완료 |

## References — 우선 확보 목록
- Yablonovitch, *J. Opt. Soc. Am.* **72**, 899 (1982)
- Brütting et al., *Phys. Status Solidi A* **210**, 44 (2013) — 70% / macroextractor 단서 포함
- Winston, Jiang, Ricketts, *Adv. Opt. Photon.* **10**, 484 (2018)
- 자기 그룹: ACS Photonics **10**, 1775 (2023); Nat. Commun. (2025)
- Freeform illumination 리뷰 / Monge–Ampère
- 각도필터 OLED (Adv. Photonics Res. 2023), angle-selective film (2022)
- 메타표면 방향성 microLED (Comm. Eng. 2025)

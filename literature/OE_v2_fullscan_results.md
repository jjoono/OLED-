# OE 원고 v2 전체 점검: 문헌 조사 결과 (2차)

- **조사일:** 2026-11-10
- **대상:** `manuscript_draft.md` (v2, 제목 *Reciprocity Constrains Angular Mixing in Passive Light-Extraction Structures*), 커밋 5aa261a 기준
- **지침:** `literature/SEARCH_BRIEF_OE_v2_fullscan.md`
- **근거 표기:**
  - **[원문]** 이번 조사에서 본문·부록을 직접 확인
  - **[원문\*]** 보조 조사에서 원문을 확인했고, 핵심 문장은 다시 확인하지 않음
  - **[초록만]** 초록·서지 정보만 확인. 수치가 있으면 원문에서 확인하지 못한 값임
- **관계 표기:** 같음 / 유사 / 보완 / 무관
- **서지:** 모든 DOI는 Crossref에서 확인했습니다(Kelemen 2001은 Crossref 미등록).

---

## 0. 요약

### 새로움 위협 판정

| 항목 | 위협도 | 한 줄 판정 |
|---|---|---|
| A1 상호반사 고유모드 | 중 | 대칭 상호반사 핵의 고유모드 분해와 적분구(고유값 1 하나)는 고전 결과. 각도 공간, 탈출, λ₂ 혼합 지표는 없음 |
| A2 광수송 연산자 스펙트럼 | **상** (주장 1) | "상반성 ⇒ 투영 입체각 내적에서 자기수반 ⇒ 실수 스펙트럼"은 그래픽스 광수송 이론의 기정 사실(Veach 1997, Soler & Subr 2025) |
| A3 q ∝ w(1−p) | **상** (분야 간) | 같은 대칭 분리형 "결합 확산 항"이 상반성·에너지 보존을 이유로 그래픽스 BRDF 모델에 있음(Ashikhmin & Shirley 2000, Kelemen 2001). 태양전지·LED·OLED에는 없음 |
| A4 Barchiesi 2026 *JOSA A* | 하 (잠정) | 원문 접근 불가. 공개 정보상 결맞음 S-행렬 방법이고 마르코프·고유값 언급 없음 |
| A5 적분구 이론 | 하 | 적분구 배율 공식은 λ₂ = 0인 무기억 극한. 비Lambertian 벽의 각도 기억은 몬테카를로로만 다룸 |
| A6 Green 2002 | 하 | 이상적 Lambertian 해석해와 재무작위화 없는 극한만 다룸 |
| A7 LED·OLED 재순환 고유값 | 중 | LED·OLED 반례 없음("OLED 광추출 분야 최초" 가능). 다만 "후기 탈출 = 스펙트럼 반경"은 복사 전달에 있음(Xu 2026) |
| A8 용어 | "trapping ratio" 중 | 다른 분야 용례가 있고, ξ가 클수록 포획이 약해져 이름과 뜻이 반대 |

### 원고에서 고쳐야 할 주장

1. **상반성 ⇒ 상세균형 ⇒ 실수 스펙트럼 (초록, 서론 마지막 문단, 2.2절).** 결과로 내세우지 말고, 알려진 성질을 B_R에 적용한 것으로 쓰고 Veach 1997, Soler & Subr 2025를 인용합니다(A2).
2. **"Reciprocity fixes q" (2.3절).** 같은 형태가 그래픽스 반사 모델에 있으므로 Ashikhmin & Shirley 2000, Kelemen & Szirmay-Kalos 2001을 인용합니다. 새로움은 "다중 패스 광추출 모델에 이 제약을 적용하고 정확도 향상을 정량화한 것"으로 한정합니다(A3).
3. **p∞의 범위 (2.4절).** 하한 p∞ ≥ min p(θ)는 비음수 행렬의 행합(열합) 상하한으로, 상호반사 이론에 같은 형태가 있습니다(Langer 1999, 명제 3). 상한은 레일리–리츠 변분 원리의 직접 적용입니다. 부등식 자체는 새롭지 않고, 준정상 분포의 탈출 확률이라는 해석과 지표화가 새롭다고 써야 합니다(A1, A2).
4. **서론 "What is missing is a common way to compare the structures themselves".** Greiner 2007이 구조별 "특성화 도표"를 이미 제시했습니다. 출발 각도와 소자 반사율에 대한 다중 패스 탈출 확률 표입니다. 이를 인정하고, 탈출과 혼합을 분리하지 않는 수치 표라는 차이를 밝혀야 합니다(B6).
5. **서론 "or not at all, as noted for partial randomization [22]".** Battaglia 2012는 단일 산란에 대한 Lambertionality factor를 정의했으므로 부정확합니다(B4).
6. **"trapping ratio" 이름 (A8).**
7. **"55–60% 정체"의 근거 (B1).** 외부 구조 55% 이상 기록이 대부분 한 그룹에서 나왔습니다. 근거 문헌 추가를 권합니다.
8. **5절 첫 문단의 가장 강한 주장 (B13).** "같은 넓이의 굴절면은 내부 휘도의 밝은 부분을 원하는 대역으로 보낼 수 없다"는 일반 명제는 Lenef 2018, Bortz 1997과 어긋납니다. 같은 넓이의 표면도 각도 투과가 이진형이면 radiance 정렬 상한에 도달할 수 있습니다. "이번에 계산한 구조 계열에서 관찰된 결과"로 범위를 좁혀야 합니다(B13).
9. **"only the slope distribution matters" (3.2절, B10).** 한 번만 상호작용하는 미세면 모델에서만 엄밀합니다. "높이장을 균일하게 확대해도 결과가 같으므로 상관 형태가 같으면 rms 기울기로만 정해진다"로 바꾸기를 권합니다.
10. **B11 후보 중 Slootsky & Forrest 2009는 MLA를 다루지 않아 근거로 쓸 수 없습니다.** Sun & Forrest 2006(원고 [9])의 "약 5 µm 이상에서 광선 추적과 일치"는 원문 확인 전에는 수치 없이 인용해야 합니다.

### 권장 서론 문장 (새로움 범위 재설정, 영문)

> Both spectral quantities, and the reciprocity properties that make them real and bounded, are standard objects of Markov-chain theory, of open dynamical systems, and of light-transport theory [24–30; Veach 1997; Soler & Subr 2025]; the reciprocal coupling of a diffuse lobe to angle-dependent reflection is likewise standard in reflectance modelling [Ashikhmin & Shirley 2000]. What is new here is their use as structure metrics for passive light extraction, the reciprocal reduced model of multi-pass extraction built on them, and the resulting map of OLED outcoupling structures.

---

## 1. 새로움 위협 확인 (A1–A8)

### A1. 상호반사 연산자의 고유모드 — 위협도 **중**

- **Koenderink, J. J.; van Doorn, A. J., "Geometrical modes as a general method to treat diffuse interreflections in radiometry," *J. Opt. Soc. Am.* 73, 843–850 (1983), doi:10.1364/JOSA.73.000843** [원문]
  - 내용: 확산(Lambertian) 상호반사의 핵이 대칭·양정치임을 보이고, 고유함수("pseudofacet")로 분해해 모드별 이득 1/(1−ρλ_k)로 해를 씀. 적분구는 0이 아닌 고유값이 1 하나뿐인 경우(식 12–17과 적분구 예).
  - 관계: **유사.** 연산자를 고유모드로 분해하고 적분구를 rank-one 경우로 보는 구조가 같음. 다만 공간 영역의 확산면 문제이고, 각도 분해·탈출·λ₂ 혼합 지표는 없음.
  - 위치·권장 문장 (2.2절): "Eigenmode expansions of interreflection operators are classical in radiometry, where an integrating sphere is the rank-one case with a single non-zero eigenvalue [Koenderink & van Doorn 1983]; here the same construction is applied to the angular return operator of an extraction structure, which escape makes sub-stochastic."
- **Langer, M. S., "When shadows become interreflections," *Int. J. Comput. Vis.* 34, 193–204 (1999), doi:10.1023/A:1008131719047** [원문]
  - 내용: 주 고유값을 핵의 최소·최대 행 적분으로 상하한함(명제 3). 구형 오목면은 0이 아닌 고유값이 (1 + cos α)/2 하나뿐(명제 1, Moon 1940 인용). 주 모드는 비음수.
  - 관계: **유사.** 원고의 하한 p∞ ≥ min p(θ)와 같은 꼴의 상하한.
  - 위치·권장 문장 (2.4절): "The lower bound p∞ ≥ min p(θ) is the angular counterpart of the row-sum bounds on the principal eigenvalue of interreflection kernels [Langer 1999]; the upper bound follows from the Rayleigh quotient in the étendue-weighted inner product."
- **Moon, P., "On interreflections," *J. Opt. Soc. Am.* 30, 195–205 (1940), doi:10.1364/JOSA.30.000195** [초록만] — 보완. 고유함수 접근의 원조. Langer와 함께 인용 가능.
- **Forsyth, D.; Zisserman, A., "Shape from shading in the light of mutual illumination," Proc. Alvey Vision Conf. 1989, 33.1–33.6, doi:10.5244/C.3.33** [원문\*] — 보완. 상호반사 핵의 대칭성과 노이만 급수 수렴. 고유값 해석 없음.

### A2. 광수송 연산자의 스펙트럼 — 위협도 **상** (주장: 상반성 ⇒ 상세균형 ⇒ 실수 스펙트럼)

- **Veach, E., *Robust Monte Carlo Methods for Light Transport Simulation*, PhD thesis, Stanford University (1997)** [원문]
  - 내용: 산란 연산자 K는 BSDF가 대칭일 때, 그리고 그때만 투영 입체각 측도에서 자기수반(정리 4.6). 일반 BSDF의 상반성을 상세균형과 Kirchhoff 법칙에서 유도(6.1–6.2절, 정리 6.1 [원문\*]). 굴절이 있으면 basic radiance L/n²로 대칭화(7장).
  - 관계: **같음.** 원고 2.2절의 상세균형 w_i B_R(i→j) = w_j B_R(j→i)는 이 성질의 이산형.
  - 위치·권장 문장 (2.2절): "Reciprocity makes the scattering operator of any passive surface self-adjoint in the projected-solid-angle (étendue) inner product [Veach 1997]; the detailed balance of B_R is the discrete form of this property."
- **Soler, C.; Subr, K., "Spectral theory of light transport operators," *ACM Trans. Graph.* 45, 1–18 (2025), doi:10.1145/3774756** [원문]
  - 내용: 국소 BRDF 연산자 K_x의 스펙트럼은 (−1, 1) 안의 실수이고 고유함수가 완비(정리 5.2.1). 실측 재질(MERL)의 "BRDF 고유값"을 계산(Fig. 8). Lambertian radiosity 연산자도 대칭화 가능(정리 5.1.1). 탈출이나 혼합 해석은 없음.
  - 관계: **같음** (실수 스펙트럼), **유사** (λ₂와 같은 대상의 스펙트럼).
  - 위치·권장 문장 (2.2절): "The spectra of reciprocal reflection operators are real [Soler & Subr 2025]; we use the second eigenvalue of the column-normalized return matrix as a measure of how fast a structure forgets the angle."
  - 주의: 반환 행렬에 대각인 소자 반사율 R(θ)만 곱하면 실수 스펙트럼이 유지됩니다. 비대각 상반 연산자 두 개를 곱하면 복소 고유값이 생길 수 있습니다(Soler & Subr 5.4절 [원문\*]). 소자 반사를 정반사로 둔 원고 가정의 범위를 밝혀야 합니다.
- **Arvo, J.; Torrance, K.; Smits, B., "A framework for the analysis of error in global illumination algorithms," Proc. SIGGRAPH '94, 75–84 (1994), doi:10.1145/192161.192179** [원문\*]
  - 내용: 반사·전파 연산자가 자기수반. 상반성 때문에 ‖K‖ ≤ 최대 반사율이고 (I−K)⁻¹의 노름은 1/(1−최대 반사율) 이하(5.3절).
  - 관계: **유사** (p∞ ≥ min p 쪽).
  - 같이 볼 문헌: Arvo, J., "The role of functional analysis in global illumination," Rendering Techniques '95, 115–126 (1995), doi:10.1007/978-3-7091-9430-0_12 [원문\*]; Arvo, PhD thesis, Yale (1995) 6.4.1절 [원문\*].
- **Baranoski, G. V. G.; Bramley, R.; Rokne, J. G., "Eigen-analysis for radiosity systems," Univ. Calgary TR 97/601/03 (1997)** [원문\*]
  - 내용: radiosity 행렬을 면적·반사율 대각 행렬로 대칭화해 실수 고유값을 얻음. 고유벡터의 물리적 의미는 "미해결 문제"라고 씀.
  - 관계: **같음** (이산 대칭화). 참고로 인용 가능.
- **Cohen, M. F.; Chen, S. E.; Wallace, J. R.; Greenberg, D. P., "A progressive refinement approach to fast radiosity image generation," *Comput. Graph.* (SIGGRAPH '88) 22, 75–84 (1988), doi:10.1145/378456.378487** [원문\*]
  - 내용: 면적 평균 반사율로 만든 "ambient" 항 1/(1−ρ_avg).
  - 관계: **유사.** 무기억 추정과 같은 역할.
  - 위치·권장 문장 (토의): "The étendue-averaged escape probability plays the role of the well-mixed 'ambient' estimate of radiosity [Cohen et al. 1988], and the trapping ratio quantifies its failure."
- **Kajiya, J. T., "The rendering equation," Proc. SIGGRAPH '86, 143–150 (1986), doi:10.1145/15922.15902** [원문\*] — 보완. 노이만 급수 수렴 조건이 스펙트럼 반경 < 1이고 항이 반사 횟수에 대응.

**고칠 곳 (A2):**
- 초록 "Reciprocity makes the return matrix of a passive structure obey detailed balance with respect to the étendue, which yields two spectral metrics"는 사실이지만 새 결과처럼 읽힙니다. 2.2절에 Veach 1997과 Soler & Subr 2025를 추가하고, 서론의 새로움 문장에서 "the constraints that reciprocity places on them"을 위 0절의 권장 문장처럼 바꾸기를 권합니다.

### A3. 축약 모델의 상반적 재분배 q ∝ w(1−p) — 위협도 **상** (분야 간), 태양전지·LED·OLED 안에서는 하

- **Ashikhmin, M.; Shirley, P., "An anisotropic Phong BRDF model," *J. Graph. Tools* 5, 25–32 (2000), doi:10.1080/10867651.2000.10487522** [원문]
  - 내용: Fresnel 가중 반사가 있으면 Lambertian 확산 항은 에너지 보존을 위해 0이어야 하므로, 확산 항을 ρ_d(k₁,k₂) = c R_d (1−R(k₁))(1−R(k₂))로 둠(식 12). R(k)는 반사 성분의 반구 반사율. 이 형태가 "에너지 보존과 상반성을 함께 지킨다"고 쓰고, Shirley et al. (Pacific Graphics 1997)에 공을 돌림.
  - 관계: **같음.** étendue 구간으로 옮기면 B_R ∝ (1−p_i)(1−p_r) w_r, 즉 q ∝ w(1−p).
  - 위치·권장 문장 (2.3절): "The reciprocal redistribution q ∝ w(1−p) is the étendue-binned form of the coupled diffuse lobe used in reciprocal, energy-conserving reflectance models [Ashikhmin & Shirley 2000; Kelemen & Szirmay-Kalos 2001]; to our knowledge it has not been imposed in multi-pass light-extraction or light-trapping models."
- **Kelemen, C.; Szirmay-Kalos, L., "A microfacet based coupled specular-matte BRDF model with importance sampling," Eurographics 2001 Short Presentations (2001)** [원문] (DOI는 Crossref에 미등록)
  - 내용: "대칭을 보장하려고" matte BRDF를 대칭 분리형 f = s (1−a_spec(θ′))(1−a_spec(θ))로 둠(2.2절). 정규화 상수 s = 1/[π(1−a_spec 평균)]는 q = w(1−p)/Σw(1−p)와 같음.
  - 관계: **같음** (정규화까지).
- **Nicodemus, F. E. et al., *Geometrical Considerations and Nomenclature for Reflectance*, NBS Monograph 160 (1977), doi:10.6028/NBS.MONO.160** [원문\*]
  - 내용: 방향–반구 반사율과 반구–방향 반사율 인자가 같음(Kirchhoff, 상반성).
  - 관계: **보완** (물리적 근거). 등방 조명에서 구조가 돌려보내는 빛의 각도 분포가 w(θ)[1−p(θ)]에 비례한다는 직관의 근거.
- **Neyts, K.; Ullan Nieto, A., "Importance of scattering and absorption for the outcoupling efficiency in organic light-emitting devices," *J. Opt. Soc. Am. A* 23, 1201–1206 (2006), doi:10.1364/JOSAA.23.001201** [원문]
  - 내용: n sinθ 구간 클래스별로 탈출을 정하고(임계각 아래 투과, 위 TIR), 산란된 빛은 "다시 등방 분포가 된다"고 가정. 금속 거울 반사율 R은 각도 무관.
  - 관계: **유사.** 각도 의존 탈출과 등방 재분배를 쓰는 OLED 축약 모델. 다만 등방 재분배가 별도 산란층 안에서 일어나므로 상반성을 위반하지 않음.
- **McTague, J.; Ratnaweera, R.; Sheldon, M., arXiv:2607.20783 (2026)** [원문\*] — 보완. 국소 전이의 미시적 가역성을 강제한 흡수 마르코프 연쇄. 상태가 공간 층이고 부분 혼합 q는 없음.
- **Rau, U.; Paetzold, U. W.; Kirchartz, T., "Thermodynamics of light management in photovoltaic devices," *Phys. Rev. B* 90, 035211 (2014), doi:10.1103/PhysRevB.90.035211** [초록만] — 무관. 각도 평균 형식(후속 논문 Abebe et al. 2018 [원문\*]의 설명).

**고칠 곳 (A3):**
- 2.3절 "Reciprocity fixes q"를 "As in reciprocal coupled diffuse–specular reflectance models [Ashikhmin & Shirley 2000; Kelemen & Szirmay-Kalos 2001], reciprocity fixes q"로 바꿉니다.
- "The common choice q ∝ w ... violates reciprocity whenever p depends on angle"에는 "when it is applied to the combined return of a structure whose escape depends on angle"를 덧붙입니다. Neyts & Nieto처럼 별도 산란층의 등방 재분배는 상반성을 지키기 때문입니다.

### A4. Barchiesi, *JOSA A* 2026 — 위협도 **하** (잠정)

- **Barchiesi, D., "Height-map S-matrix formalism for light trapping in textured multilayer stacks: polarization coupling, energy conservation, and photovoltaic application," *J. Opt. Soc. Am. A* (early posting, 2026), doi:10.1364/JOSAA.600340** [서지만]
  - 원문: Optica 로그인·보안 확인 때문에 열지 못함. Crossref·OpenAlex·Semantic Scholar 모두 초록 없음. 5절 링크 참조.
  - 동반 소프트웨어 SMaK 2 (Zenodo, doi:10.5281/zenodo.21276131) 설명에는 면내 파수 모드별 편광 S-행렬과 텍스처 회절 모드가 있고, 마르코프 연쇄·고유값·부분 무작위화·상반성 언급은 없음 [초록만].
  - 관계: **보완** (추정). 결맞음 파동 광학으로 단일 만남 연산자를 공급할 수 있는 방법.
  - 원문을 열면 확인할 단어: redistribution, eigen, spectral radius, Lambertian, reciprocity, detailed balance, multiple passes.

### A5. 적분구 이론 — 위협도 **하**

- **Jacquez, J. A.; Kuppenheim, H. F., "Theory of the integrating sphere," *J. Opt. Soc. Am.* 45, 460 (1955), doi:10.1364/JOSA.45.000460** [초록만]; **Goebel, D. G., "Generalized integrating-sphere theory," *Appl. Opt.* 6, 125 (1967), doi:10.1364/AO.6.000125** [초록만]
  - 관계: **보완.** 배율 M = ρ/[1−ρ(1−f)]는 완전 확산 벽을 가정한 무기억(λ₂ = 0) 극한. 공식 자체는 Tang et al., *Appl. Opt.* 57, 1581 (2018), doi:10.1364/AO.57.001581, 식 (1)에서 확인 [원문\*].
  - 위치·권장 문장 (2.4절): "For a Lambertian integrating sphere the redistribution operator has rank one [Koenderink & van Doorn 1983], so N_mix = 1 and p∞ = p, which recovers the classical sphere multiplier [Jacquez & Kuppenheim 1955; Goebel 1967]."
- **Hanssen, L., "Effects of non-Lambertian surfaces on integrating sphere measurements," *Appl. Opt.* 35, 3597 (1996), doi:10.1364/AO.35.003597** [초록만]; **Prokhorov, A. V.; Mekhontsev, S. N.; Hanssen, L. M., "Monte Carlo modeling of an integrating sphere reflectometer," *Appl. Opt.* 42, 3832 (2003), doi:10.1364/AO.42.003832** [초록만]
  - 관계: **보완.** 비Lambertian 벽 때문에 무기억 공식에서 벗어나는 정도를 몬테카를로로만 다룸. 고유값 해석 없음.
- **Tardy, H. L., "Matrix method for integrating-sphere calculations," *J. Opt. Soc. Am. A* 8, 1411 (1991), doi:10.1364/JOSAA.8.001411** [초록만] — 보완. BRDF를 다루는 구역 행렬법, 스펙트럼 분석 언급 없음.
- **열복사 교환:** Naraghi, M. H. N.; Chung, B. T. F., *J. Heat Transfer* 106, 690–698 (1984), doi:10.1115/1.3246740 [초록만]; Billings, R. L. et al., *Numer. Heat Transfer A* 19, 313–326 (1991), doi:10.1080/10407789108944852 [초록만]
  - 관계: **보완.** 흡수 상태가 있는 마르코프 연쇄로 밀폐 공간의 복사 교환을 계산(양방향 반사 포함). 스펙트럼 반경·준정상 분포 언급은 초록에 없음.
  - 위치·권장 문장 (서론, 선택): "Markov chains with absorbing states are established for radiative exchange in enclosures, including bidirectionally reflecting walls [Naraghi & Chung 1984; Billings et al. 1991]."

### A6. Green 2002 — 위협도 **하**

- **Green, M. A., "Lambertian light trapping in textured solar cells and light-emitting diodes: analytical solutions," *Prog. Photovolt.* 10, 235–241 (2002), doi:10.1002/pip.404** [원문]
  - 내용: 이상적 Lambertian 광포획의 해석해. 뒷면에서 매번 재무작위화하는 구조(Fig. 1a)와 재무작위화 없이 각도 분포가 유지되는 구조(Fig. 1b)를 구분. Kirchhoff 상반성을 근거로 LED에도 적용.
  - 관계: **보완.** 부분 무작위화, 각도 의존 탈출, 고유값은 없음.
  - 위치·권장 문장 (서론 2문단): "Closed-form multi-pass solutions exist for ideal Lambertian randomization, including light-emitting devices [Green 2002]; N_mix interpolates between this limit and complete angle retention."

### A7. LED·OLED 광자 재순환 모델의 고유값·스펙트럼 반경 — 위협도 **중**

- **검색:** OpenAlex 제목·초록 논리 검색 3종과 웹 검색 약 10회.
  - 조합: (spectral radius OR eigenvalue OR Perron) × (light extraction OR outcoupling) × (OLED OR LED); "photon recycling" × (eigenvalue OR spectral radius OR Markov OR Neumann series); "Markov chain" × LED/OLED.
  - 결과: LED·OLED 반례는 찾지 못했습니다. "OLED·LED 광추출에서 반환 연산자의 스펙트럼을 처음 쓴다"는 표현은 유지할 수 있습니다.
- **Xu, H., "Photon escape from slab Thomson media: a scattering-order-resolved recursive formalism for Comptonization applications," *Mon. Not. R. Astron. Soc.* 550, stag1392 (2026), doi:10.1093/mnras/stag1392** [원문]
  - 내용: 높은 산란 차수에서 탈출 비율 P_n/P_{n−1}이 재귀 연산자의 스펙트럼 반경으로 수렴하고, 같은 고유모드가 극한 각도 분포를 정함(초록, 식 102–103).
  - 관계: **유사** (p∞와 같은 발상). 산란 차수 단위이고 구조의 통과 단위가 아님.
  - 위치·권장 문장 (2.4절): "That late-order escape is set by the spectral radius of the transport recursion, with a limiting angular distribution, is known in radiative transfer [Xu 2026; Alves-Pereira et al. 2007]."
- **Alves-Pereira, A. R. et al., "Radiation trapping in 1D using the Markov chain formalism," *Eur. J. Phys.* 28, 1105–1124 (2007), doi:10.1088/0143-0807/28/6/008** [원문\*] — 유사. 후기 세대가 고정된 기본 모드로 수렴하고 세대마다 일정 비율로 감쇠.
- **Esposito, L. W.; House, L. L., "Radiative transfer calculated from a Markov chain formalism," *Astrophys. J.* 219, 1058 (1978), doi:10.1086/155870** [원문\*] — 유사. (천정각, 깊이) 상태의 흡수 연쇄, fundamental matrix로 해를 구함. 고유값 해석 없음.
- **Kutayiah, A. R. et al., "Markov chains for modeling complex luminescence, absorption, and scattering in nanophotonic systems," *Opt. Express* 29, 4249 (2021), doi:10.1364/OE.416670** [원문\*] — 보완. 정상 상태(고유값 1 고유벡터), 각도별 Kirchhoff·상세균형, LSC 적용.
- **Krieger, A.; Wolf, S., "The scattering order problem in Monte Carlo radiative transfer," *Astron. Astrophys.* 645, A143 (2021), doi:10.1051/0004-6361/202039133** [원문\*] — 보완. 비음수 분포가 지배 고유상태로 수렴(부록 A.3).

**고칠 곳 (A7):**
- 마르코프 연쇄로 광자 탈출을 다룬 첫 연구, 또는 후기 탈출과 스펙트럼 반경을 처음 연결한 연구라고 쓰면 안 됩니다.
- 차별점은 다음 네 가지로 쓰는 것이 안전합니다.
  - 수동 반환 연산자의 극각 상태공간
  - 상반성 대칭화로 얻는 실수 스펙트럼
  - λ₂ → N_mix
  - min p ≤ p∞ ≤ p

### A8. 용어 충돌

| 용어 | 기존 사용 | 위험 | 권고 |
|---|---|---|---|
| trapping ratio | 광학 용례 없음. OpenAlex 제목 8건은 모두 다른 분야: 퇴적물 trapping ratio(*Exp. Therm. Fluid Sci.* 2017, doi:10.1016/j.expthermflusci.2017.05.022), CO₂ 모세관 trapping ratio(*Energy Procedia* 2011, doi:10.1016/j.egypro.2011.02.375) [초록만]. 광학에는 "light trapping", 광섬유 "trapping efficiency"가 있음 | 중 (ξ가 클수록 포획이 약해져 뜻이 반대로 읽힘) | "late-pass escape ratio" 등으로 바꾸기를 권장 |
| angular mixing number | 정확한 구문 0건(OpenAlex, 웹). 축열조의 "MIX number"(Davidson et al., *J. Sol. Energy Eng.* 116, 94 (1994))만 인접 | 하 | 유지. N_mix = t_rel(Levin & Peres)임을 명시 |
| late-pass escape probability | 0건 | 하 | 유지. 준정상 분포(Darroch–Seneta)로 정의 |
| late-time escape probability | 0건 | 하 | 쓰지 않음(통과가 아니라 시간) |
| asymptotic escape probability | 선 복사 전달의 점근 탈출 확률(Puetter, *ApJ* 251, 446 (1981)), 확률 동역학(Qiao & Duan 2015) [초록만] | 중 | 쓰지 않음 |
| return matrix | 광학 용례 없음. 금융의 "asset return matrix" | 하 | 유지. 처음 나올 때 "(reflection–redistribution) matrix"를 병기. OPTOS는 "redistribution matrix"를 씀 |

---

## 2. 원고 문장별 근거 확인 (B1–B17)

### B1. 서론 1문단: 외부 구조 OLED EQE가 55–60%에서 정체

이전 조사(Nature 원고 Fig. 1a 원자료)에서 원문으로 확인한 기록입니다 [원문(Fig. 1a 조사)]. 모두 단일 소자, 매크로 렌즈 없음, 적분구(IS)·각도 적분(AR)·LCF로 측정한 값입니다.

| 연도 | EQE | 구조 | 측정 | 서지 |
|---|---|---|---|---|
| 2015 | 41.5% | 광추출 시트 (TADF) | IS | Kaji, H. et al., *Nat. Commun.* 6, 8476 (2015), doi:10.1038/ncomms9476 |
| 2015 | 44.7% | 상부 전극 위 고굴절 유기 MLA (TE) | LCF | Kim, J.-B. et al., *Org. Electron.* 17, 115–120 (2015), doi:10.1016/j.orgel.2014.12.006 |
| 2018 | 56.3% | SiO₂ 입자 산란 필름 | AR | Song, J. et al., *Nat. Commun.* 9, 3207 (2018) (원고 [6]) |
| 2019 | 61.7% | AgNW 전극 위·기판 뒷면 ZnO 나노패턴 (내부+외부, 백색 유연) | IS | Li, W. et al., *Adv. Opt. Mater.* 7, 1900985 (2019), doi:10.1002/adom.201900985 |
| 2020 | 44.3% | 유기 나노텍스처 미세구조 필름 (청색) | AR+IS | Han, Y. et al., *Org. Electron.* 87, 105892 (2020), doi:10.1016/j.orgel.2020.105892 |
| 2021 | 48% | MLA 패턴 사파이어 기판 (고굴절 기판) | AR | Lee, J. et al., *Adv. Opt. Mater.* 9, 2002182 (2021), doi:10.1002/adom.202002182 |
| 2023 | 58% | 초박형 기판에 내장한 역MLA | AR | Kim, J. et al., *ACS Photonics* 10, 1775 (2023) (원고 [35]) |
| 2025 | 57.2% | 유한 크기 준평면 구조 | AR | Kim, M. et al., *Nat. Commun.* 16, 11606 (2025) (원고 [34]) |

참고(내부 구조): Jeon, S. et al., *Sci. Rep.* 5, 8685 (2015), doi:10.1038/srep08685 — 50.65%; Jeon, S. et al., *Adv. Opt. Mater.* 6, 1701349 (2018), doi:10.1002/adom.201701349 — 백색 44.3%.

- 판단: 외부 구조만으로 55% 이상인 기록 세 편(Song 2018, Kim 2023, Kim 2025)은 모두 Yoo 그룹입니다. 다른 그룹의 외부 구조 기록은 41–48% 수준이고, 내부+외부 결합인 Li 2019가 61.7%입니다.
- 위치·권장 문장 (서론 1문단): "With external and combined structures, the highest reported EQEs of single-unit OLEDs lie near 55–62% [6,34,35; Li et al. 2019], a ceiling recently traced to the loss that returned light suffers on each round trip through the device [15]."

### B2. 서론 2문단: "randomization fraction" [20] Deckman 1983

- **Deckman, H. W.; Roxlo, C. B.; Yablonovitch, E., *Opt. Lett.* 8, 491–493 (1983), doi:10.1364/OL.8.000491** [원문, 저자 공개본]
  - 확인: 텍스처가 빛을 충분히 산란하지 못하면 입사광의 일부 β만 무작위화된다고 보고, "randomization fraction β"를 F_exp = βF_en + (1−β)F_flat로 정의함(식 4).
  - 판단: 원고 인용은 **정확**합니다. 수정 불필요.

### B3. 서론 2문단: 내부 분포가 Lambertian처럼 보이기까지의 통과 수 [21] Manzoor 2020

- **Manzoor, S. et al., *J. Appl. Phys.* 127, 063104 (2020), doi:10.1063/1.5131173** [원문]
  - 확인: 통과마다 남은 빛으로 정규화한 각도 분포와 Lambertian의 RMSD(식 5)로 평가. 실제 랜덤 피라미드는 2 pass, 이상적 랜덤 피라미드는 8 pass에 Lambertian 분포가 됨.
  - 판단: 원고 표현은 **정확**합니다.

### B4. 서론 2문단: 박막 태양전지의 부분 무작위화는 정량화되지 않음 [22] Battaglia 2012

- **Battaglia, C. et al., *J. Appl. Phys.* 112, 094504 (2012), doi:10.1063/1.4761988** [원문]
  - 확인: 산란체가 Lambertian에 얼마나 가까운지를 재는 "Lambertionality factor" a와 탈출 분율 b를 정의함. 매 통과 같은 분포인 경우와 첫 통과만 비Lambertian인 경우를 계산함. 연속 산란이 Lambertian에 다가가는 과정의 현실적 기술은 "yet to be developed"라고 씀.
  - 판단: "or not at all"은 **공정하지 않습니다.** 단일 산란은 정량화했고, 통과별 수렴만 미해결로 남겼습니다.
  - 권장 수정 (서론 2문단): "…, or by a single-event Lambertianity factor, with the approach to a Lambertian distribution over successive passes left open [22]."

### B5. 2.3절: Lambertian 재분배(q ∝ w)가 흔한 선택이라는 근거

- **Yablonovitch 1982 (원고 [16]); Schnitzer et al. 1993 (원고 [17])** — 완전 무작위화(photon gas) 가정의 원전.
- **Green 2002** [원문] (A6) — 이상적 Lambertian 다중 패스 해석해, LED 포함. 관계: **유사**.
- **David 2013** (원고 [19]) [원문, 1차 조사] — 거친 표면 산란을 정반사 성분 R_s(θ)와 확산 성분 R_d로 나누고, 실용 광선 추적에서는 확산 성분을 Lambertian으로 근사해도 된다고 씀(III-C절, 각주 7 [원문\*]). 관계: **유사**.
- **Neyts & Ullan Nieto 2006** [원문] (A3) — 산란 후 등방 재분배를 가정한 OLED 기하 모델. 관계: **유사**.
- 단서: 복사 전달 모델(Shiang & Duggal 2004, Bathelt et al. 2007 [초록만])은 등방 대신 Henyey–Greenstein 비대칭 인자 g를 씁니다.
- 위치·권장 문장 (2.3절): "Reduced models usually redistribute returned light isotropically, as in the photon-gas picture [16,17; Green 2002] and in geometrical OLED models with isotropic scattering [Neyts & Ullan Nieto 2006]; applied to the combined return of a structure whose escape depends on angle, this choice violates reciprocity."

### B6. 서론 2문단: 각도 분해 다중 패스 모델 계열 (누락 확인)

- **Greiner, H., "Light extraction from organic light emitting diode substrates: simulation and experiment," *Jpn. J. Appl. Phys.* 46, 4125–4137 (2007), doi:10.1143/JJAP.46.4125** (원고 [2]) [원문]
  - 확인: 구조를 특성화하는 도표 두 가지를 계산함. ① 출발 방출각 θ와 OLED 스택 반사율(각도 무관 가정)에 대한 다중 패스 총 탈출 확률, ② 기판 방출각별 공기 중 각도 분포. "광원 분포와 무관하게 구조의 효과를 평가하는 표"라고 설명함(3.3절). 결론: 반사광을 충분히 재지향하는 한 구조의 세부 형태는 부차적이며, 달성 가능한 전달률은 대략 스택 평균 반사율로 정해짐.
  - 관계: **유사**, 중요. 원고 서론 "What is missing is a common way to compare the structures themselves"와 직접 닿습니다.
  - 권장 수정 (서론 2문단): "Greiner characterized structures by tabulated escape probabilities versus emission angle and stack reflectance [2], but such tables are numerical and do not separate escape from redistribution."
- **Shiang, J. J.; Duggal, A. R., "Application of radiative transport theory to light extraction from organic light emitting diodes," *J. Appl. Phys.* 95, 2880–2888 (2004), doi:10.1063/1.1644037** [초록만; 방법은 GE 최종 보고서 Duggal et al., "OLEDs for General Lighting," DOE (2004), doi:10.2172/836831로 원문\*]
  - 내용: 극각을 25–60개 채널로 나눈 이산 방향 복사 전달과 경계 반사 행렬. 음극은 각도 무관 정반사 R_OLED 하나.
  - 관계: **유사.** 가장 중요한 누락 문헌.
  - 위치·권장 문장 (서론 2문단): "Angle-discretized radiative-transfer models with specular boundary matrices were applied to OLED substrates early on [Shiang & Duggal 2004]."
- **Altazin, S. et al., "Simulations, measurements, and optimization of OLEDs with scattering layer," SID Symp. Dig. Tech. Pap. 46, 564–567 (2015), doi:10.1002/sdtp.10478** [원문\*]
  - 내용: SETFOS 4.1이 쌍극자 방출을 두꺼운 층과 결정론적 net-radiation 알고리즘으로 결합하고, Mie 계산으로 얻은 각도 분해 BSDF를 씀. 쌍극자+광선 추적과 비교 검증.
  - 관계: **유사.** OPTOS 옆에 인용 권장.
- **Kovačič, M. et al., "Analysis and optimization of light outcoupling in OLEDs with external hierarchical textures," *Opt. Express* 29, 23701 (2021), doi:10.1364/OE.428021** [원문\*]
  - 내용: 통과별 분석(Fig. 7a). 첫 통과 투과는 텍스처와 무관하게 약 45%, 이득은 약 10회 통과 후 포화하고 그 수준은 스택 반사율이 정함.
  - 관계: **유사.** 서론 또는 통과 수 논의에 인용. 같은 저자의 결합 모델: *ACS Photonics* 5, 422–430 (2018), doi:10.1021/acsphotonics.7b00874 [초록만].
- **Mladenovski, S. et al., "Integrated optical model for organic light-emitting devices," *J. Appl. Phys.* 109, 083114 (2011), doi:10.1063/1.3576114** [원문\*] — 보완. 각도·편광 분해 비간섭 다중 반사, 평면 기판만.
- **Vandersteegen, P. et al., "Employing a 2D surface grating to improve light out coupling of a substrate emitting organic LED," Proc. SPIE 6486, 64860H (2007), doi:10.1117/12.701344** [원문\*] — 유사. 기판 각도 분포 × 각도 의존 격자 투과(다중 반사 포함, 비간섭 기판). 스택이 완전 흡수이면 격자 이득이 거의 없음.
- 그 밖: Bathelt et al., *Org. Electron.* 8, 293–299 (2007), doi:10.1016/j.orgel.2006.11.003 [초록만]; Chen, K.-Y. et al., *Opt. Express* 18, 3238 (2010), doi:10.1364/OE.18.003238 [초록만] — 보완.

### B7. 3.1절: CPS 모델의 복굴절 확장 원전

- **Wasey, J. A. E.; Safonov, A.; Samuel, I. D. W.; Barnes, W. L., "Effects of dipole orientation and birefringence on the optical emission from thin films," *Opt. Commun.* 183, 109–121 (2000), doi:10.1016/S0030-4018(00)00806-3** [초록만]
  - 내용: 평면 다층 쌍극자 모델을 복굴절 층으로 확장한 원전. 광축이 계면에 수직인 경우만 다룸(OLED에 해당).
  - 관계: **같음.**
- **Penninck, L.; De Visschere, P.; Beeckman, J.; Neyts, K., "Dipole radiation within one-dimensional anisotropic microcavities: a simulation method," *Opt. Express* 19, 18558 (2011), doi:10.1364/OE.19.018558** [원문\*]
  - 내용: 단축 매질에서 쌍극자 장의 평면파 전개와 산란 행렬. 광축 방향 제한 없음.
  - 관계: **같음** (일반형).
- **Moon, C.-K.; Kim, S.-Y.; Lee, J.-H.; Kim, J.-J., "Luminescence from oriented emitting dipoles in a birefringent medium," *Opt. Express* 23, A279 (2015), doi:10.1364/OE.23.00A279** [초록만]
  - 내용: 도핑 박막과 OLED로 검증. OLED에서 가장 많이 쓰임.
  - 관계: **같음.**
- 선택: Callens, M. K.; Yokoyama, D.; Neyts, K., "Anisotropic materials in OLEDs for high outcoupling efficiency," *Opt. Express* 23, 21128 (2015), doi:10.1364/OE.23.021128 [원문\*] — 보완 (이방성 수송층에 적용).
- 위치·권장 문장 (3.1절): "Birefringent layers are treated as uniaxial with the optic axis along the normal, using the anisotropic extension of the CPS model [Wasey et al. 2000; Penninck et al. 2011] as implemented for OLEDs [Moon et al. 2015]."

### B8. 3.2절: 산란층과 출사면을 adding method로 결합

- **van de Hulst, H. C., *Multiple Light Scattering: Tables, Formulas, and Applications*, Vols. 1–2 (Academic Press, 1980)** [미확인, 도서] — 같음. adding–doubling의 표준 원전(1권 4장).
- **Plass, G. N.; Kattawar, G. W.; Catchings, F. E., "Matrix operator theory of radiative transfer. 1: Rayleigh scattering," *Appl. Opt.* 12, 314 (1973), doi:10.1364/AO.12.000314** [초록만] — 같음. 두께가 다른 층을 결합하는 행렬 연산자 이론.
- **Prahl, S. A.; van Gemert, M. J. C.; Welch, A. J., "Determining the optical properties of turbid media by using the adding–doubling method," *Appl. Opt.* 32, 559 (1993), doi:10.1364/AO.32.000559** [초록만] — 같음. 주변과 굴절률이 다른 슬랩과 유리 경계.
- **Prahl, S. A., "The adding-doubling method," in *Optical-Thermal Response of Laser-Irradiated Tissue* (Plenum, 1995), pp. 101–129, doi:10.1007/978-1-4757-6092-7_5** [원문\*]
  - 내용: 경계 조건을 "필요한 경계 조건을 부과하는 층을 하나 더 더하는" 방식으로 포함함(5.3.4절). 계면이 자체 반사·투과 행렬(예: Fresnel)을 가짐.
  - 관계: **같음.** 원고의 "산란층 + 평면 출사면" 결합 방식과 정확히 일치.
- 보완: Hansen, J. E.; Travis, L. D., "Light scattering in planetary atmospheres," *Space Sci. Rev.* 16, 527–610 (1974), doi:10.1007/BF00168069 [초록만].
- 위치·권장 문장 (3.2절): "Scattering slabs and the exit face are combined with the adding method [van de Hulst 1980; Plass et al. 1973], the interface entering as a layer with its own Fresnel reflection and transmission operators [Prahl 1995]."

### B9. 3.2절: 높이장 몬테카를로 광선 추적의 표준 인용

- **Greiner 2007** (원고 [2]) [원문] — 같음 (OLED). OLED 내부 각도 분포를 입력으로 미세 굴절·산란 구조를 몬테카를로 광선 추적(3.1–3.2절).
- **McIntosh, K. R.; Baker-Finch, S. C., "OPAL 2: rapid optical simulation of silicon solar cells," Proc. 38th IEEE PVSC, 265–271 (2012), doi:10.1109/PVSC.2012.6317616** [초록만] — 유사. 텍스처 광선 추적과 계면 계수를 분리해 계산.
- **Baker-Finch, S. C.; McIntosh, K. R., "Reflection of normally incident light from silicon solar cells with pyramidal texture," *Prog. Photovolt.* 19, 406–416 (2011), doi:10.1002/pip.1050** [초록만] — 유사.
- **Holst, H. et al., "Application of a new ray tracing framework to the analysis of extended regions in Si solar cell modules," *Energy Procedia* 38, 86–93 (2013), doi:10.1016/j.egypro.2013.07.253** [초록만] — 유사 (Daidalos).
- 참고: Brendel, R., "SUNRAYS: a versatile ray tracing program for the photovoltaic community," Proc. 12th EU PVSEC (1994) — DOI 없음, 쪽 미확인.
- 위치·권장 문장 (3.2절): "Microlens arrays and textures are ray-traced on their height fields with the interface coefficients applied at each intersection, as in OLED substrate-extraction [Greiner 2007] and silicon-texture [McIntosh & Baker-Finch 2012] ray tracing."

### B10. 3.2절: 가우스 상관 무작위 표면, rms 기울기만 중요 (광선광학 유효 조건)

- **Ogilvy, J. A., "Wave scattering from rough surfaces," *Rep. Prog. Phys.* 50, 1553–1608 (1987), doi:10.1088/0034-4885/50/12/001** [원문\*]
  - 내용: 접평면(Kirchhoff) 근사의 조건은 곡률 반경 a에 대해 k a cos³θ ≫ 1이고, 높이나 기울기 크기에는 명시적 제한이 없음. 기하광학 그림자 효과는 k → ∞에서 성립하고 입사 기울기와 rms 기울기의 비에만 의존.
  - 관계: **같음.**
- **Thorsos, E. I., "The validity of the Kirchhoff approximation for rough surface scattering using a Gaussian roughness spectrum," *J. Acoust. Soc. Am.* 83, 78–92 (1988), doi:10.1121/1.396188** [초록만]
  - 내용: 정확한 적분방정식 해와 비교. 상관 길이가 가장 중요한 유효성 변수. rms 기울기각 약 20° 이하에서 정확도를 정량화.
  - 관계: **같음.** 다만 "kℓ ≳ 6" 같은 구체 수치는 원문 확인 전 인용하지 말 것.
- **Tang, K.; Dimenna, R. A.; Buckius, R. O., "Regions of validity of the geometric optics approximation for angular scattering from very rough surfaces," *Int. J. Heat Mass Transf.* 40, 49–59 (1997), doi:10.1016/S0017-9310(96)00073-7** [초록만] — 같음. 무작위 프로파일에서 기하광학을 정확한 전자기 해와 직접 비교.
- 도서: Beckmann, P.; Spizzichino, A., *The Scattering of Electromagnetic Waves from Rough Surfaces* (Pergamon, 1963; Artech House, 1987) [미확인].
- 척도 불변성 표현 주의:
  - 원고 "only the slope distribution matters"는 한 번만 상호작용하는 미세면 모델에서만 엄밀합니다.
  - 더 안전한 표현: "Ray optics has no length scale, so the results are invariant under uniform scaling of the height field; for a fixed correlation shape they depend on the height and correlation length only through the rms slope."
  - 근거: Barrick, D. E., "Rough surface scattering based on the specular point theory," *IEEE Trans. Antennas Propag.* 16, 449–454 (1968), doi:10.1109/TAP.1968.1139220 [초록만] (정반사점 단면적이 기울기 통계로 정해짐).
- 위치·권장 문장 (3.2절 또는 범위): "The ray treatment requires radii of curvature and correlation lengths well above the wavelength [Ogilvy 1987; Thorsos 1988; Tang et al. 1997]."

### B11. 범위: 마이크로미터 MLA에 광선광학 적용

- **Sun, Y.; Forrest, S. R., "Organic light emitting devices with enhanced outcoupling via microlenses fabricated by imprint lithography," *J. Appl. Phys.* 100, 073106 (2006), doi:10.1063/1.2356904** (원고 [9]) [초록만]
  - 내용: MLA를 몬테카를로와 FDTD로 함께 계산함(초록에서 확인).
  - 미확인 수치: 렌즈 지름이 λ에 가까워지면 추출 효율이 떨어지고, 지름 약 5 µm(약 10λ) 이상에서 3D 몬테카를로 광선 추적 값과 일치한다는 내용(검색 결과 요약에서만 확인). 원문 Fig. 1b를 확인한 뒤 수치를 인용해야 함.
  - 관계: **같음** (원문 확인 시).
- **Park, J. S.; Lee, J. W.; Park, D. J.; Ko, J.-H., "Correlation between the size of microlenses and the outcoupling efficiency of organic light emitting diodes studied by finite-difference time-domain simulation," *J. Korean Phys. Soc.* 74, 885–891 (2019), doi:10.3938/jkps.74.885** [초록만]
  - 내용: 2D FDTD만 사용. 렌즈 크기가 λ보다 작아지면 평면 소자 수준으로 떨어짐. 광선 추적과의 비교는 없음.
  - 관계: **유사.** 요청서 후보 중 실재 확인.
- **Slootsky, M.; Forrest, S. R., "Full-wave simulation of enhanced outcoupling of organic light-emitting devices with an embedded low-index grid," *Appl. Phys. Lett.* 94, 163302 (2009), doi:10.1063/1.3116644** [초록만]
  - 내용: 내장 저굴절 격자의 전파동 해석. MLA나 광선–파동 비교 없음.
  - 관계: **무관.** 요청서 후보이지만 B11 근거로 쓸 수 없음.
- 보완: Park, J. S.; Ko, J.-H.; Lee, J. W., "Finite-difference time-domain simulation of fractal-like microlens arrays for high outcoupling efficiency of organic light-emitting diodes," *J. Nanosci. Nanotechnol.* 19, 1695–1698 (2019), doi:10.1166/jnn.2019.16175 [초록만] (2D FDTD, 지름 0.5–0.2 µm에서 평면 이하로 떨어짐). Lenef et al. 2018 [원문\*]: 지름 0.5–1.5 µm 밀집 반구의 FDTD 투과율이 광선 한계에 근접.
- 3D 광선 추적과 전파동을 직접 비교한 다른 연구는 찾지 못했습니다.
- 위치·권장 문장 (범위): "Full-wave studies show the microlens gain falling as the lens size approaches the wavelength [Sun & Forrest 2006; Park et al. 2019]; the lenslets considered here are many wavelengths across."

### B12. 3.3절: 대리모델 전역 탐색과 패턴 탐색

- **Gutmann, H.-M., "A radial basis function method for global optimization," *J. Glob. Optim.* 19, 201–227 (2001), doi:10.1023/A:1011255519438** [초록만] — 같음.
- **Regis, R. G.; Shoemaker, C. A., "A stochastic radial basis function method for the global optimization of expensive functions," *INFORMS J. Comput.* 19, 497–509 (2007), doi:10.1287/ijoc.1060.0182** [초록만] — 같음. MATLAB surrogateopt 문서가 이 방법(삼차 RBF, 가중치 순환 merit 함수)을 인용.
- **Torczon, V., "On the convergence of pattern search algorithms," *SIAM J. Optim.* 7, 1–25 (1997), doi:10.1137/S1052623493250780** [초록만] — 같음.
- **Audet, C.; Dennis, J. E., "Analysis of generalized pattern searches," *SIAM J. Optim.* 13, 889–903 (2002), doi:10.1137/S1052623400378742** [초록만] — 같음. MATLAB patternsearch 문서가 인용(문서에는 2003년으로 표기).
- MATLAB Global Optimization Toolbox의 사용 버전을 소프트웨어 인용으로 추가하기를 권장.
- 위치·권장 문장 (3.3절): "Each campaign used a radial-basis-function surrogate (surrogateopt [Gutmann 2001; Regis & Shoemaker 2007]) followed by generalized pattern search (patternsearch [Torczon 1997; Audet & Dennis 2002])."

### B13. 5절 1문단: 같은 넓이의 굴절면은 내부 휘도의 밝은 부분을 원하는 대역으로 보낼 수 없다

- **판정: 쓰인 그대로는 근거가 없고, 일반 명제로는 선행 문헌과 어긋납니다.**
- **Lenef et al. 2018** [원문]
  - 거시적으로 균일한 넓이 A의 추출 표면, 즉 광원과 같은 넓이의 경우를 다룸.
  - 비균일 radiance에 대한 표면 추출 상한은 가장 밝은 입사 방향만 골라 투과시키는 이진 각도 투과에서 도달된다고 보임(식 12–15). 기하·회절·산란 표면 모두에 적용되는 열역학 논증.
  - 예: 높은 각의 광선을 받아들이는 표면이면 0.546까지 가능(평면 0.163).
- **Bortz, Shatz & Ries 1997** (원고 [31]) [원문] — 최대 전달은 가장 큰 radiance가 표적 위상공간을 우선 채울 때 달성됨. 상한만 다루고, 같은 넓이 굴절면이 이를 달성할 수 있는지는 말하지 않음.
- **Ries, H.; Shatz, N.; Bortz, J.; Spirkl, W., "Performance limitations of rotationally symmetric nonimaging devices," *J. Opt. Soc. Am. A* 14, 2855 (1997), doi:10.1364/JOSAA.14.002855** [초록만]
  - 동반 논문: Ries, H. et al., "Consequences of skewness conservation for rotationally symmetric nonimaging devices," Proc. SPIE 3139, 47–58 (1997), doi:10.1117/12.290221 [원문\*].
  - 내용: 회전대칭 장치에서는 skewness별로 전달 가능한 étendue가 정해짐. 넓이가 같은 두 원판은 skewness 분포가 정확히 맞으므로, 원판 대 원판 전달에서는 skewness가 étendue 이상의 제약을 더하지 않음.
  - 참고: 요청서의 제목 "Consequences of skewness conservation…"은 JOSA A 논문이 아니라 SPIE 동반 논문의 제목입니다.
- **Benítez, P.; Miñano, J. C.; Blen, J., "Squeezing the étendue," in *Illumination Engineering* (Wiley/IEEE, 2013), pp. 71–99, doi:10.1002/9781118462539.ch3** [초록만] — 같은 넓이의 무초점 렌즈렛 배열은 한 방향의 콜리메이션을 높이면 다른 방향을 같은 비율로 낮춤. 위상공간 모양은 바꿀 수 있어도 étendue는 줄일 수 없음.
- **Winston, R.; Miñano, J. C.; Benítez, P., *Nonimaging Optics* (Elsevier, 2005)** [미확인, 도서] — étendue와 skew 불변량(부록 A).
- 결론: étendue와 skewness 보존은 같은 넓이의 표면이 대역으로 보낼 수 있는 위상공간의 양을 정하지만, 어떤 내부 방향이 그 위상공간을 채우는지는 정하지 않습니다. "균일한 것처럼 행동한다"는 결론은 무작위화 극한(Yablonovitch 1982)이나 원고가 계산한 구조 계열의 관찰로만 뒷받침됩니다.
- 권장 수정 (5절 1문단): "Being coextensive with the source, a refractive film cannot enlarge the étendue available to a band; it can only select, through its local slopes, which internal directions are refracted into it [Lenef et al. 2018]. The radiance-ordered bound [31–33] requires the band to be filled from the brightest internal directions only. The rotationally symmetric microlenses and textures studied here fill each band from a broad range of internal angles, so in our simulations they extract close to what a uniform internal radiance would give."

### B14. 5절 2문단: 산란층은 두꺼울수록 p를 잃음

- **Altazin et al. 2015** [원문\*]
  - 내용: 외부 PET 필름의 입자 농도에 따라 출력이 0.25%(탁도 67%)에서 최대가 되고, 그 뒤로는 입자 흡수와 필름 반사율 증가로 줄어듦(4절, Fig. 10).
  - 관계: **같음.**
  - 위치·권장 문장 (5절 2문단): "For particle-loaded foils an optimum near 67% haze is predicted, beyond which back-reflection and absorption reduce extraction [Altazin et al. 2015]."
- **Shiang, J. J.; Faircloth, T. J.; Duggal, A. R., "Experimental demonstration of increased organic light emitting device output via volumetric light scattering," *J. Appl. Phys.* 95, 2889–2895 (2004), doi:10.1063/1.1644038** [초록만; 데이터는 GE 보고서로 원문\*]
  - 내용: 향상 배율이 scatterance S에 대해 최대를 가짐(ZrO₂ S ≈ 2, 형광체 필름 S ≈ 3, 최대 1.41배). S가 크면 빛이 손실 있는 음극으로 되돌아감.
  - 관계: **같음.**
- **Ang, P. Y. et al., "Inside or outside: evaluation of the efficiency enhancement of OLEDs with applied external scattering layers," *Sci. Rep.* 9, 18601 (2019), doi:10.1038/s41598-019-54640-x** [원문\*]
  - 내용: 적분구 EQE는 탁도와 함께 증가(최대 +24%)하지만, 탁도가 가장 큰 필름은 측면 확산 때문에 정면 휘도가 15–27% 줄어듦.
  - 관계: **보완.** 작은 픽셀의 적분구 측정이 감소를 가릴 수 있다는 단서.

### B15. 5절 3문단: 무작위·불규칙 MLA 권장

- **Lim, J. et al., "Enhanced out-coupling factor of microcavity organic light-emitting devices with irregular microlens array," *Opt. Express* 14, 6564–6571 (2006), doi:10.1364/OE.14.006564** [원문\*]
  - 내용: 회절을 피하려고 5 µm × 2.5 µm 렌즈를 무작위로 배치. 마이크로캐비티 색 변화가 줄고 정면 EQE가 약 1.8배(캐비티 효과 포함).
  - 관계: **같음.**
- **Park, J.; Hyun, W.; Park, K.; Im, S.; Park, O. O., "Enhanced light outcoupling efficiency in organic light-emitting devices using irregular microlenses fabricated with 3D colloidal arrays," *Sci. Adv. Mater.* 6, 2370–2377 (2014), doi:10.1166/sam.2014.2185** [초록만] — 같음. 불규칙 MLA로 약 +40%, 간섭 억제로 원래 스펙트럼 유지.
- **Galeotti, F.; Mróz, W.; Scavia, G.; Botta, C., "Microlens arrays for light extraction enhancement in organic light-emitting diodes: a facile approach," *Org. Electron.* 14, 212–218 (2013), doi:10.1016/j.orgel.2012.10.034** [초록만]; **Han, J. et al., "Enhanced outcoupling in down-conversion white organic light-emitting diodes using imprinted microlens array films with breath figure patterns," *Sci. Technol. Adv. Mater.* 20, 35–41 (2019), doi:10.1080/14686996.2018.1551040** [원문\*] — 유사. breath-figure로 만든 준규칙 MLA.
- 위치·권장 문장 (5절 3문단): "Randomly placed or irregular microlenses have been used on OLEDs to avoid diffraction and angular colour shift [Lim et al. 2006; Park et al. 2014]."

### B16. 범위: 소자 반사가 각도를 보존(정반사)한다는 가정

- **Shiang & Duggal 2004 / GE 최종 보고서** [원문\*]
  - 내용: 스택 반사를 각도 무관 정반사로 가정하고, 계산된 Al/유리 반사율 0.84–0.87(0–80°)로 정당화. 비정반사 항은 텍스처가 있는 음극에서만 필요하다고 봄. R_OLED = 0.79를 정반사 측정 부속장치로 측정.
  - 관계: **같음.**
- **Greiner 2007** [원문] — 스택 반사율을 입사각과 무관하게 두었고, 박막 이론으로 계산한 반사율을 써도 결과가 크게 달라지지 않았다고 씀. 관계: **같음**.
- **Levell et al. 2013** (원고 [7]) [원문\*] — 소자 반사율을 스펙트럼 평균 스칼라(41–86%)로 다룸. 정반사·확산을 나눌 수 있는 장치로 쟀지만 분리 값은 보고하지 않음. 관계: **보완**.
- **Riedel, B. et al., "Improving the outcoupling efficiency of indium-tin-oxide-free organic light-emitting diodes via rough internal interfaces," *Opt. Express* 18, A631–A639 (2010), doi:10.1364/OE.18.00A631** [원문\*] — 보완. 거친 계면 OLED의 정반사·확산 반사를 측정했고 차이가 작다고만 보고.
- 평면 금속 음극 OLED의 확산 반사 비율을 측정해 보고한 문헌은 **찾지 못했습니다.**
- 위치·권장 문장 (범위): "As in earlier substrate models, the device is treated as an angle-preserving reflector [Shiang & Duggal 2004; Greiner 2007]; diffuse reflection from planar metal cathodes is small, and its inclusion would add a non-diagonal factor to the return operator."

### B17. 범위: 두꺼운 기판의 비간섭 처리

- **Katsidis, C. C.; Siapkas, D. I., "General transfer-matrix method for optical multilayer systems with coherent, partially coherent, and incoherent interference," *Appl. Opt.* 41, 3978 (2002), doi:10.1364/AO.41.003978** [초록만] — 같음.
- **Centurioni, E., "Generalized matrix method for calculation of internal light energy flux in mixed coherent and incoherent multilayers," *Appl. Opt.* 44, 7532 (2005), doi:10.1364/AO.44.007532** [초록만] — 같음 (사입사 포함).
- 대안: Troparevsky, M. C. et al., *Opt. Express* 18, 24715 (2010), doi:10.1364/OE.18.024715 [초록만] — 유사 (무작위 위상 평균).
- OLED 전용:
  - **Mladenovski et al. 2011** (B6) [원문\*] — 기판이 파장보다 훨씬 두꺼워 기판→공기 추출을 비간섭으로 다뤄야 한다고 명시. η = η_OLED→sub · η_sub→air로 쓰고 기판 내 반사를 세기로 합산. 관계: **같음**.
  - **Bulović, V. et al., "Weak microcavity effects in organic light-emitting devices," *Phys. Rev. B* 58, 3730–3740 (1998), doi:10.1103/PhysRevB.58.3730** [초록만] — 유사. OLED 층과 기판의 빛을 파동광학과 광선광학으로 함께 기술.
- 위치·권장 문장 (3.1절): "The substrate is treated incoherently [Katsidis & Siapkas 2002; Centurioni 2005], as is standard in OLED outcoupling models [Bulović et al. 1998; Mladenovski et al. 2011]."

---

## 3. 직접 비교 연구 (C1–C3)

### C1. 같은 소자에서 주기 MLA, 무작위 MLA, 산란층 비교

- 세 구조를 한 소자에서 비교한 실험은 **찾지 못했습니다.** 두 구조 비교는 있습니다.
- **Han et al. 2020** (B1 표) [원문\*] — 같은 청색 인광 소자에서 무작위 산란 미세구조(탁도 53%)와 상용 반구 MLA 필름을 비교. EQE 32.8% → 44.3%(무작위) 대 40.0%(MLA). 기판 빛 추출 계산값은 48% 대 30%. 관계: **유사**.
- **Vandersteegen et al. 2007** (B6) [원문\*] — 같은 9 mm² OLED에서 2D 격자와 MLA가 적분구 측정으로 비슷한 성능(약 +30%). 관계: **유사**.
- **Rahadian, F.; Imai, K.; Fujieda, I., *Opt. Eng.* 46, 124001 (2007), doi:10.1117/1.2821869** [초록만] — 프리즘 필름 4종과 확산 필름 3종을 두 OLED에 비교, 최대 1.60배. 반사율이 높은 OLED가 항상 더 큰 이득. 관계: **유사**.

### C2. 각도 보존 반사와 확산 반사를 실험으로 구분한 연구

- **Ang et al. 2019** [원문\*] — 필름의 직접·확산 투과·반사 스펙트럼과 405 nm 각도별 산란 측정. 관계: **유사**.
- **Altazin et al. 2015** [원문\*] — 입자층의 각도별 투과를 측정하고 Mie 계산과 비교해 BSDF 입력으로 씀. 관계: **유사**.
- **GE 보고서, Levell 2013, Riedel 2010** [원문\*] — 소자 정반사율이나 정반사·확산 반사를 측정했지만 분리 결과는 제한적. 관계: **보완**.
- OLED용 MLA 필름의 고니오미터 BSDF 측정은 **찾지 못했습니다.**

### C3. 2024–2026 OLED 광추출 리뷰

- **Kang, K. et al., "Nanostructures in organic light-emitting diodes: principles and recent advances in the light extraction strategy," *Laser Photonics Rev.* 18, 2400547 (2024), doi:10.1002/lpor.202400547** [초록만, 오픈 액세스] — 내부 나노구조 중심, 산란층 포함.
- 외부 구조 비중이 큰 대안: **Shinar, R.; Shinar, J., "Light extraction from organic light emitting diodes (OLEDs)," *J. Phys. Photonics* 4, 032002 (2022), doi:10.1088/2515-7647/ac6ea4** [초록만]; **Salehi, A. et al., "Recent advances in OLED optical design," *Adv. Funct. Mater.* 29, 1808803 (2019), doi:10.1002/adfm.201808803** [초록만].
- 위치·권장 문장 (서론 1문단): "Recent reviews summarize internal and external extraction strategies [Kang et al. 2024; Shinar & Shinar 2022]."

---

## 4. 원문을 확인하지 못한 논문 (PDF 링크)

| 항목 | 논문 | 링크 | 비고 |
|---|---|---|---|
| A4 | Barchiesi, *JOSA A* (2026) | https://doi.org/10.1364/JOSAA.600340 | 조기 게재본, Optica 로그인 필요 |
| B6·B14 | Shiang & Duggal, *JAP* 95, 2880 (2004) | https://pubs.aip.org/aip/jap/article-pdf/95/5/2880/18707068/2880_1_online.pdf | AIP는 기관 접근 가능하나 PDF가 브라우저 탭에서 열리지 않음 |
| B14 | Shiang, Faircloth & Duggal, *JAP* 95, 2889 (2004) | https://doi.org/10.1063/1.1644038 | 같은 이유 |
| A3 | Rau, Paetzold & Kirchartz, *PRB* 90, 035211 (2014) | https://doi.org/10.1103/PhysRevB.90.035211 | 판정 영향 낮음 |
| A1 | Moon, *JOSA* 30, 195 (1940) | https://doi.org/10.1364/JOSA.30.000195 | 판정 영향 낮음 |
| A5 | Jacquez & Kuppenheim (1955); Goebel (1967); Tardy (1991); Hanssen (1996); Prokhorov et al. (2003) | https://doi.org/10.1364/JOSA.45.000460 ; https://doi.org/10.1364/AO.6.000125 ; https://doi.org/10.1364/JOSAA.8.001411 ; https://doi.org/10.1364/AO.35.003597 ; https://doi.org/10.1364/AO.42.003832 | 판정 영향 낮음 |
| B5·B6 | Bathelt et al., *Org. Electron.* 8, 293 (2007) | https://doi.org/10.1016/j.orgel.2006.11.003 | |
| B6 | Kovačič et al., *ACS Photonics* 5, 422 (2018) | https://doi.org/10.1021/acsphotonics.7b00874 | |
| B15 | Park et al. (2014); Galeotti et al. (2013) | https://doi.org/10.1166/sam.2014.2185 ; https://doi.org/10.1016/j.orgel.2012.10.034 | |
| C1 | Rahadian et al., *Opt. Eng.* 46, 124001 (2007) | https://doi.org/10.1117/1.2821869 | |
| C3 | Kang et al. (2024); Shinar & Shinar (2022); Salehi et al. (2019) | https://doi.org/10.1002/lpor.202400547 ; https://doi.org/10.1088/2515-7647/ac6ea4 ; https://doi.org/10.1002/adfm.201808803 | Kang 2024는 오픈 액세스 |
| B11 | Sun & Forrest, *JAP* 100, 073106 (2006) | https://doi.org/10.1063/1.2356904 | **우선 확인.** "지름 약 5 µm 이상에서 광선 추적과 일치" 수치(Fig. 1b). Chrome에서도 초록 페이지로 돌아감 |
| B13 | Ries, Shatz, Bortz & Spirkl, *JOSA A* 14, 2855 (1997) | https://doi.org/10.1364/JOSAA.14.002855 | skewness 제약의 원문 확인 |
| B7 | Wasey et al. (2000); Moon et al. (2015) | https://doi.org/10.1016/S0030-4018(00)00806-3 ; https://doi.org/10.1364/OE.23.00A279 | Moon 2015는 오픈 액세스이나 Optica 보안 확인 |
| B8 | Prahl et al., *Appl. Opt.* 32, 559 (1993) | https://doi.org/10.1364/AO.32.000559 | |
| B10 | Thorsos, *JASA* 83, 78 (1988); Tang et al., *IJHMT* 40, 49 (1997) | https://doi.org/10.1121/1.396188 ; https://doi.org/10.1016/S0017-9310(96)00073-7 | 유효 조건 수치 인용 전 확인 |
| B11 | Park et al., *JKPS* 74, 885 (2019) | https://doi.org/10.3938/jkps.74.885 | |
| B17 | Katsidis & Siapkas (2002); Centurioni (2005) | https://doi.org/10.1364/AO.41.003978 ; https://doi.org/10.1364/AO.44.007532 | 표준 방법 인용이라 판정 영향 낮음 |

사용자가 제공한 PDF로 원문을 확인한 논문: Koenderink & van Doorn 1983, Green 2002, Greiner 2007, Neyts & Ullan Nieto 2006.

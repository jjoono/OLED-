# References — OLED microlens-array limits paper

검증 상태: URL은 웹검색으로 실재 확인. 서지 상세(권/페이지)는 확인된 것만 표기, 미확인은 표시.
그룹 자기인용(⚑)은 self-overlap 관리 대상.

---

## Claim 1 — "freeform/역설계가 '더'를 약속하며 부상" (✅ 검증됨, 유지)

- Kim et al., "Inverse design of organic light-emitting diode structure based on deep
  neural networks," *Nanophotonics* 10 (2021). — OLED 구조 역설계(DNN).
  https://onlinelibrary.wiley.com/doi/10.1515/nanoph-2021-0434
- "Dual-Task Optimization Method for Inverse Design of RGB Micro-LED Light Collimator"
  (FDTD 역설계, ±20° 내 결합효율 ~30→60%).
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11820347/
- "Inverse design for scalable photonic systems," *Nature Reviews Materials* (2026). — 역설계 부상 리뷰.
  https://www.nature.com/articles/s41578-026-00915-5
- "Design of freeform microlens arrays with prescribed luminance distributions for
  MicroLED optical packaging," *Appl. Opt.* 64, 7875 (2025). — 프리폼 MLA로 목표 배광 설계.
  https://opg.optica.org/ao/abstract.cfm?uri=ao-64-27-7875

## Claim 2 — étendue 트레이드오프는 조각조각 기존 (⚠️ 절대표현 금지, 부분 인용)

- "Design Methodology for High Brightness Projectors," *J. Display Technol.* 4, 1 (2008). —
  étendue로 size vs efficiency 트레이드오프. http://www.optimesa.com/docs/Design%20Methodology%20for%20High%20Brightness%20Projectors.pdf
- "Directed display architecture" (US patents 10288884 / 11163165 / 10761327). —
  MLA 디스플레이 밝기 vs 면적/eyebox étendue 트레이드오프.
  https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10288884
- → 이들은 "개별 트레이드오프는 이미 알려짐"의 근거. 통합 achievable-region/possibility 지도는 미발견
  (그래서 "아무도 없음" 대신 "통합·실험검증된 바 없음"으로 서술).

## Claim 2 심화 — 트레이드오프/한계 선행 (차별화 문단 근거)

이미 발표된 "조각들" (우리는 이것들을 통합·손실분리·실험검증으로 넘어섬):
- MLA는 escape cone **밖** 광을 추출하되 cone **안** 추출은 줄고, back-surface **재활용(다중반사)**으로 보상
  — "Analysis of light out-coupling from microlens array," *Opt. Commun.* (2011).
  https://www.sciencedirect.com/science/article/abs/pii/S0030401811006511
- MLA-OLED 휘도효율 평가 기법(기초): *Opt. Express* 12, 5777 (2004).
  https://opg.optica.org/oe/fulltext.cfm?uri=oe-12-23-5777&id=81744
- 효율 ↔ **이미지 blur(공간)** 트레이드오프 (각도 조향 아님): MLA-film 연구.
  https://www.tandfonline.com/doi/full/10.1080/15980316.2018.1531073
- 방향성 ↔ 총효율 트레이드오프, 충전물질에 따라 opening-angle 최적 이동
  — micro-Horn collimator µLED, arXiv:2412.14027. https://arxiv.org/html/2412.14027
  (단, collimator/능동집속 구조 = 주기 굴절렌즈 클래스 밖)
- 회절/공진 기반 방향성 OLED (다른 클래스): *Adv. Photonics Res.* (2023).
  https://advanced.onlinelibrary.wiley.com/doi/10.1002/adpr.202200143
- 메타표면 phosphor µLED 방향성 (다른 클래스): *ACS Nano* (2024).
  https://pubs.acs.org/doi/10.1021/acsnano.4c13472

→ 차별화: 위는 (a) 효율, (b) 효율↔공간blur, (c) collimator/메타표면의 방향성↔효율을
   *각각* 다룸. **주기 굴절렌즈 어레이 × 자발광 확장광원**에 대해 효율×방향성을 하나의
   **닫힌 achievable region**으로 통합하고, **기하/손실을 분리**하고, **제작으로 봉투를 검증**한
   선행은 없음.

## 각도선택 필터 + 광재활용 (⚠ 선행 확정 — "제안" 금지, "검증"으로만 사용)

우리 Fig3(b)(c) 계산의 대상 구조는 **이미 알려진 개념**이다. 새 플랫폼 제안으로
쓰면 즉시 선행에 걸린다. 반드시 "우리 이론이 이 알려진 거동을 예측/설명한다"는
**검증 재료**로만 배치할 것.

- Buhl et al., "Resonance-Based Directional Light Emission from OLEDs,"
  *Adv. Photonics Res.* (2023). — **핵심**: DBR을 입사각 필터로 쓰면 저각 방출을
  억제하지만 "light interacts more often with absorbing materials inside the device"
  로 손실이 커진다고 명시. 우리가 계산한 트레이드오프가 이미 서술돼 있음.
  https://advanced.onlinelibrary.wiley.com/doi/10.1002/adpr.202200143
- "Using angle-selective optical film to enhance the light extraction of a thin-film
  encapsulated 3D reflective pixel for OLED displays" (2022). — OLED에 각도선택 필름
  직접 적용. https://pubmed.ncbi.nlm.nih.gov/36558597/
- "Broadband Angular Selectivity of Light at the Nanoscale," arXiv:1512.02761. —
  각도선택성 + photon recycling 리뷰. https://arxiv.org/pdf/1512.02761
- "Methods and apparatus for broadband angular selectivity" (US 10073191). —
  각도선택 필터 = TiO2/SiO2 다층막 정의까지 동일.
  https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10073191
- "High efficient OLED displays ... air-gapped bridges on quantum dot patterns for
  optical recycling," *Sci. Rep.* 7, 43063 (2017). — OLED 광재활용 실측 이득(적색 58.2%).
  https://www.nature.com/articles/srep43063
- LCD 백라이트 BEF/DBEF: off-axis 광을 저손실 캐비티로 되돌려 재활용 — 동일 원리의
  상용화 선례(정면휘도 ~1.6x).

## Reciprocity / Generality (기하 절반의 보편성 근거)

- U. Rau, "Reciprocity relation between photovoltaic quantum efficiency and
  electroluminescent emission of solar cells," *Phys. Rev. B* 76, 085303 (2007).
  — 준평형에서 EL 방출 ↔ 광전 흡수(EQE)의 상반성. 우리 "방향성 봉투가 방출=흡수 동형"의
  일반화 Kirchhoff 근거. https://doi.org/10.1103/PhysRevB.76.085303
- E. Yablonovitch, "Statistical ray optics," *J. Opt. Soc. Am.* 72, 899 (1982).
  — 방출↔흡수를 detailed balance로 연결, 4n²/2n², 손실이 유한 천장을 만든다.
  https://opg.optica.org/josa/abstract.cfm?uri=josa-72-7-899
- 점검 결과(요지): **동형** = 각도 재분배/étendue/zero-mean-slope 봉투(방출=흡수 acceptance).
  **비동형** = 총효율 천장(useful 채널이 탈출 vs 활성층으로 다름), 스펙트럼 적분 최적, 동작점.

## Foundational (정독)

- Yablonovitch 1982 (위) — 효율천장=손실 소관, 재활용, 4n².
- "Nonimaging optics: a tutorial," *Adv. Opt. Photon.* 10, 484 (2018). — étendue 보존,
  brightness theorem, sine-law 집광한계 (방향성 상한의 정석 언어).
  https://opg.optica.org/aop/fulltext.cfm?uri=aop-10-2-484&id=389885

## 메타표면은 실제로 조향됨 (스코핑 가드레일 근거)

- "Enhanced and directional electroluminescence from MicroLEDs using metallic or dielectric
  metasurfaces," *Communications Engineering* (2025). — 메타표면이 위상경사로 방향성 부여
  → 우리 "못 넘는다" 주장을 refractive height-continuous array로 한정해야 하는 근거.
  https://www.nature.com/articles/s44172-025-00401-w

---

## 핵심 선행/자기중복 (risk report와 연동)

- ⚑ Kim, Kim, Park, Song, Kim, Moon, Yoo, "Toward Near-Foldable Surface Light Sources...
  Ultrathin Substrates Embedded with Micron-Scale Inverted Lens Arrays," *ACS Photonics*
  10, 1775 (2023). — **자기 그룹**. trans-scale(CPS+BSDF+LightTools) MLA, EQE 58%.
  → 방법·venue·효율숫자 선점. [PDF: Kim2023_ACSPhotonics_IMLA_foldable.pdf]
  https://pubs.acs.org/doi/10.1021/acsphotonics.3c00017
- ⚑ "Near-planar light outcoupling structures with finite lateral dimensions for
  ultra-efficient and optical crosstalk-free OLED displays," *Nature Communications*
  s41467-025-66538-6 (2025). — **자기 그룹**. 스택+구조 공동설계 48% EQE.
  https://www.nature.com/articles/s41467-025-66538-6
- 반구(H/D→0.5) 최적: *Opt. Commun.* (microlens parameters).
  https://www.sciencedirect.com/science/article/abs/pii/S1566119916303585
- 주기 MLA 조향 불가(회절 프레임): "High-resolution beam steering using microlens arrays,"
  *Opt. Lett.* 31, 2861 (2006). https://opg.optica.org/ol/abstract.cfm?uri=ol-31-19-2861
- 확장광원 freeform + étendue blur: *Optica* 3, 840 (2016).
  https://opg.optica.org/optica/fulltext.cfm?uri=optica-3-8-840&id=348107
- cross-scale 3D OLED 픽셀 EM 시뮬: *Org. Electron.* (2022).
  https://www.sciencedirect.com/science/article/abs/pii/S1566119922003068
- 상용 trans-scale 툴: Fluxim Setfos. https://www.fluxim.com/emission-module

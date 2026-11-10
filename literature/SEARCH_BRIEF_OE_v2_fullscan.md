# 문헌 조사 요청서 (2차): OE 원고 v2 전체 점검

## 0. 배경

### 원고
- 대상: Optics Express 투고 예정 원고 `manuscript_draft.md` (영문), `manuscript_draft_ko.md` (한국어).
- 제목: *Reciprocity Constrains Angular Mixing in Passive Light-Extraction Structures: Two Spectral Metrics for OLED Outcoupling*

### 1차 조사 이후 바뀐 핵심 주장
1. **상반성 → 상세균형.** 반환 행렬이 에텐듀 가중치에 대해 상세균형 $w_i B_R(i\to j) = w_j B_R(j\to i)$을 만족하고, 그래서 스펙트럼이 실수다.
2. **상반적 축약 모델.** 축약 모델의 재분배는 상반성 때문에 반드시 $q \propto w(1-p)$여야 한다. 흔히 쓰는 램버시안 재분배 $q \propto w$는 상반성을 위반한다.
   - 342개 경우에서 오차가 평균 0.81/최대 10.5 %p에서 0.30/3.0 %p로 줄어든다.
3. **두 번째 지표.** $p_\infty = 1-\rho(B_R)$는 준정상 분포의 탈출 확률이다.
   - 상반성으로부터 $\min p(\theta) \le p_\infty \le p$가 성립한다. 레일리 몫으로 증명한다.
   - 트래핑 비 $\xi = p_\infty/p$를 정의한다.
4. **혼합의 한계.** 상반성 때문에 혼합은 구조가 투과시키지 않는 각도를 열 수 없다. 예를 들어 평면 계면은 완전히 혼합해도 η가 0.38에서 0.49로만 오른다.

### 이미 원고에 반영한 것 (다시 조사하지 않아도 됨)
- **1차 조사 결과:** survey 파일의 문헌들
- **이번에 추가한 문헌 (Crossref 확인 완료):**
  - 준정상 분포: Darroch–Seneta 1965, Collet–Martínez–San Martín 2013
  - 누설계: Altmann–Portela–Tél RMP 2013, PRL 2013
  - 유전체 캐비티와 부분 탈출: Ryu et al. PRE 2006, Clauß et al. PRE 2019
  - MLA 형상: Sun & Forrest JAP 2006, Peng et al. JDT 2005, Hwang et al. pssa 2014
  - 무질서 텍스처와 주름: Zhang–Biswas–Shinar JOSA B 2021, Koo et al. Nat. Photonics 2010
  - 방법론: Clarke & Parry 1985 (Helmholtz 상반성), Chance–Prock–Silbey 1978, Henyey–Greenstein 1941, Holm 1979
- **알려진 사실:** OPTOS(Eisenlohr OE 2015)가 재분배 행렬의 상반성 대칭을 이미 언급했다. 원고는 이를 인정하고, 새로움을 그 결과로 나오는 것들(실수 스펙트럼, q 형태, p∞ 범위, 보고 지표)로 한정했다.

## 1. 최우선: 새로움 위협 확인 (원문 확인 필수)

| # | 확인할 것 | 위협 시나리오 | 검색어 예시 |
|---|---|---|---|
| A1 | **상호반사(interreflection) 연산자의 고유모드.** Koenderink & van Doorn, "Geometrical modes as a general method to treat diffuse interreflections in radiometry", *JOSA* 73, 843 (1983) 및 후속 | 확산 상호반사를 고유모드로 분해한 고전 연구입니다. λ₂·ρ 해석이 이미 있으면 반드시 인용해야 합니다. | interreflection eigenmodes radiometry; light transport operator eigenvalue spectral radius |
| A2 | **컴퓨터 그래픽스의 광수송 연산자 스펙트럼.** Arvo, Torrance, Smits 1994 "A framework for the analysis of error in global illumination", Neumann 급수 수렴률 = 스펙트럼 반경, BRDF 상반성 + 고유값 | p∞ = 1 − ρ와 같은 개념이 그래픽스에 이미 있는지 확인합니다. | light transport operator spectral radius reciprocity; Neumann series convergence global illumination |
| A3 | **축약 모델 재분배 q ∝ w(1−p)의 선행 사례.** "부분 램버시안 + 상반성", "reciprocal Lambertian fraction" 등 | 원고의 핵심 새로움 중 하나입니다. 같은 형태가 태양전지·LED 모델에 있었는지 확인합니다. | partially Lambertian redistribution reciprocity; Lambertian fraction texture model detailed balance |
| A4 | **JOSA A 2026 "Height-map S-matrix formalism for light trapping in textured multilayer stacks"** (doi:10.1364/JOSAA.600340로 추정) | 1차 웹 검색에서 원문에 접근하지 못했습니다. 재분배 행렬, 고유값, 부분 무작위화를 다루는지 확인합니다. | 제목 그대로 |
| A5 | **적분구(integrating sphere) 이론.** sphere multiplier = 1/[1−ρ(1−f)], 비-램버시안 코팅의 각도 기억 | 무기억 공식과 형식이 같습니다. 각도 기억을 고유값으로 다룬 적분구 문헌이 있는지 확인합니다. | integrating sphere non-Lambertian coating multiplier angular |
| A6 | **Green, *Prog. Photovolt.* 10, 235 (2002)** "Lambertian light trapping in textured solar cells and light-emitting diodes: analytical solutions" (doi:10.1002/pip.404) | LED까지 다룬 해석적 다중 패스 해입니다. 부분 무작위화나 각도 의존 탈출을 다루는지 확인합니다. | 원문 직접 확인 |
| A7 | **OLED·LED 광자 재순환(photon recycling) 모델에서 고유값·스펙트럼 반경을 쓴 사례** | "OLED 분야 최초"라는 표현을 쓸 수 있는지 판단합니다. | photon recycling LED matrix eigenvalue; OLED substrate mode recycling eigenvalue |
| A8 | **용어 충돌.** "trapping ratio", "angular mixing number", "late-pass escape probability" | 다른 분야에서 다른 뜻으로 쓰이면 이름을 바꿔야 합니다. | 각 용어를 따옴표로 검색 |

## 2. 원고 문장별 근거 확인 (인용 정확성)

| # | 원고 위치 | 원고 문장 (요지) | 현재 인용 | 요청 |
|---|---|---|---|---|
| B1 | 서론 1문단 | 외부 광추출 구조를 써도 OLED EQE가 대체로 55–60%에서 정체 | [15] (Nature 원고) | 외부 구조로 보고된 **최고 EQE 기록 목록** (2015–2026, 녹색·백색 포함, 구조 종류와 EQE). 독립 근거 2–3편 추가용 |
| B2 | 서론 2문단 | 무작위화를 "randomization fraction"으로 정량화 | [20] Deckman 1983 | Deckman 1983이 실제로 이 개념을 정의했는지 확인. 아니면 맞는 원전 (예: "Lambertian fraction"을 쓴 태양전지 논문) |
| B3 | 서론 2문단 | 내부 분포가 램버시안처럼 보이기까지의 패스 수로 정량화 | [21] Manzoor 2020 | 이 표현이 정확한지 원문 확인 |
| B4 | 서론 2문단 | 박막 태양전지 부분 무작위화는 정량화하지 않은 채로 남음 | [22] Battaglia 2012 | 이 요약이 공정한지 확인 |
| B5 | 서론 2문단 | 램버시안 재분배($q \propto w$)를 쓰는 기존 모델 | 없음 | 각도 의존 탈출 + 램버시안 재분배를 쓴 **LED·OLED 모델 2–3편** (예: David 2013, Greiner 2007의 해당 여부). 2.3절 "흔한 선택"의 근거 |
| B6 | 서론 2문단 | 각도 분해 다중 패스 모델 계열 | OPTOS, David 2013, [15] | OLED 기판 모드 재순환을 각도 분해로 다룬 다른 모델 (SETFOS/Fluxim 광선추적, Greiner, Mladenovski, Riedel 등). 누락 여부 |
| B7 | 3.1절 | CPS 모델의 복굴절 확장 | [37] CPS + [15,34] | 복굴절·이방성 층 CPS 확장의 **원전** 1–2편 |
| B8 | 3.2절 | 산란층을 출사면과 adding method로 결합 | 없음 | adding-doubling 원전 (van de Hulst 1980, Prahl 1993 등) |
| B9 | 3.2절 | 높이장 몬테카를로 광선추적 | 없음 | 텍스처·MLA용 높이장 광선추적의 표준 인용 (예: Brendel "Sunrays", Baker-Finch & McIntosh, 미세면 광선추적) |
| B10 | 3.2절 | 가우스 상관 무작위 표면, rms 기울기만 중요 (광선광학) | Koo 2010 | **광선광학(기하광학) 근사의 유효 조건** (상관 길이 ≫ λ, 기울기 범위). Beckmann–Spizzichino, Kirchhoff 근사 유효성 문헌 |
| B11 | 3.2절 / 범위 | 마이크로미터 MLA에 광선광학 적용 | 없음 | MLA OLED에서 **광선 추적 대 FDTD 비교** (렌즈 크기에 따른 유효성). 후보: Park et al. JKPS 74, 885 (2019), Slootsky & Forrest APL 2009 |
| B12 | 3.3절 | 대리모델 전역 탐색 + 패턴 탐색 | 없음 | 대리모델 RBF 최적화(Gutmann 2001 등)와 패턴 탐색(Audet & Dennis 2006 등) 표준 인용 |
| B13 | 5절 1문단 | **확장 광원과 같은 넓이를 갖는 굴절면은 내부 휘도의 밝은 부분을 원하는 대역으로 사상할 수 없다** | 없음 (2.5절에 [31–33]) | **가장 강한 주장이라 근거가 꼭 필요합니다.** 에텐듀·skewness 보존, 확장 광원 + 동일 넓이 광학계의 한계를 엄밀히 다룬 문헌 (Ries et al. skewness 1997, Winston–Miñano–Benítez *Nonimaging Optics* 2005, "étendue squeezing" 등) |
| B14 | 5절 2문단 | 산란층은 두께가 늘수록 p를 잃음 | 없음 | 산란층 두께·탁도에 따른 추출 효율 최적점을 보고한 OLED 실험·시뮬레이션 1–2편 |
| B15 | 5절 3문단 | 무작위 MLA·물결면 권장 | Koo 2010 | **무작위·불규칙 MLA를 OLED에 적용한 실험** (breath figure, 자기조립, 불규칙 MLA로 색 변화 억제 등) 2–3편 |
| B16 | 범위 | 소자 반사가 각도를 보존(정반사)한다는 가정 | 없음 | 금속 전극 OLED의 확산 반사 성분이 작다는 근거, 또는 Levell [7] 같은 소자 반사 연구에서의 취급 |
| B17 | 범위 | 기판이 두꺼워 비간섭(incoherent) 처리 가능 | 없음 | 두꺼운 기판의 비간섭 처리 표준 인용 (예: Centurioni 2005 generalized TMM, 혹은 OLED 기판 모드 문헌) |

## 3. 참고: 직접 비교 연구 (있으면 Discussion에 인용)
- **C1.** 같은 기판 위에서 주기 MLA, 무작위 MLA, 산란층을 같은 조건으로 측정한 OLED 실험 비교
- **C2.** 다중 패스 광추출에서 "각도 보존 반사"와 "확산 반사"를 실험으로 구분한 연구 (예: 각도 분해 반사율, 외부 구조의 BSDF 측정)
- **C3.** 2024–2026년 OLED 외부 광추출 리뷰 (최신 리뷰 1편은 서론 인용 권장)

## 4. 결과 형식 (1차와 동일)
- 항목 번호(A1–A8, B1–B17, C1–C3)별로 정리합니다.
- 문헌마다 다음을 적습니다:
  - 서지 (저자, 제목, 저널, 권, 쪽, 연도, DOI)
  - **원문 확인 여부**: [원문] / [초록만]
  - 원고와의 관계: 같음 / 유사 / 보완 / 무관
  - 원고에 넣을 위치와 권장 문장 (영문 1문장)
- **A 항목은 위협도를 상, 중, 하로 판정합니다.** "상"이면 원고의 어느 주장을 약하게 고쳐야 하는지 적습니다.
- 파일명: `literature/OE_v2_fullscan_results.md`

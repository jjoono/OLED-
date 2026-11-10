# 문헌 조사 지침: 수동 광추출 구조의 각도 혼합(angular mixing) 지표와 일반화 radiance 상한

## 배경 (조사 목적)
Optics Express 투고 예정 시뮬레이션 논문의 신규성 확인용 조사입니다. 논문의 핵심 주장은 두 가지입니다.

1. **일반화 radiance(étendue) 상한**
   - 수동(선형, 이득 없음, 파장 변환 없음) 외부 광추출 구조가 한 번 통과에서 각 polar band로 보낼 수 있는 출력의 상한입니다.
   - 실제 OLED 기판 안의 비등방 radiance 분포에서, 밝은 위상공간 영역부터 채우는 "radiance 정렬" 방식으로 정의합니다.
   - 기판 안의 빛이 등방이면 Yablonovitch 상한(1/n², Lambertian 분할)으로 돌아갑니다.
   - 현재 결과: MLA·산란층·평면 계면 모두 이 상한의 35–80%에 머뭅니다.
2. **각도 혼합 지표**
   - 구조의 반사 재분배 행렬(통과마다 각도 분포를 넘기는 마르코프 행렬)의 두 번째 고유값 λ₂로 정의합니다.
   - N_mix = 1/(1−λ₂)는 "각도를 잊는 데 필요한 통과 수"입니다.
   - 구조를 (탈출확률 p, 각도 혼합 N_mix)로, 소자를 왕복 손실 A′로 기술하면 총 추출 효율과 배광 구성이 정해진다는 틀입니다.
   - η = p/[p+(1−p)A′] (Kim, Lee, Yoo, Gather, 왕복 손실 논문)은 λ₂ = 0인 특수한 경우입니다.
   - 현재 결과: 산란층 N_mix ≈ 1–2, 반구 MLA ≈ 8, 일부 freeform MLA ≈ 17, 평면은 ∞입니다. 또 비선택적 수동 구조 전체에서 40–60° 선택도가 0.33–0.36 창에 고정됩니다.

**알고 싶은 것**: 위 1, 2와 같거나 매우 비슷한 정의·결과가 이미 발표됐는지, 그리고 반드시 인용해야 할 가장 가까운 선행 연구가 무엇인지.

## 조사 범위

### A. 반드시 확인 (가장 중요)
1. 각도 재분배 행렬 또는 마르코프 연쇄의 **고유값, spectral gap, 혼합 시간(mixing time)**을 광포획·광추출 구조의 특성 지표로 쓴 연구
   - 태양전지 광포획: OPTOS(Tucher, Eisenlohr, Fraunhofer ISE), angular matrix framework, RayFlare, GenPro
   - LED / OLED / QD-LED / 페로브스카이트 LED 광추출
   - 발광 태양광 집광기(LSC)의 Markov chain 모델
2. **"degree of randomization", "Lambertianity", "angular randomization factor", "ergodicity"**를 정량 지표로 정의한 광포획 연구
   - 카오스·에르고딕 공진기를 이용한 광포획(chaotic cavity light trapping)
   - Lambertian 산란체와의 차이를 수치화한 연구
3. **Yablonovitch / étendue / brightness theorem의 일반화**
   - 비등방 내부 radiance, 각도 제한(angular restriction)에서의 상한
   - band(각도 구간)별 상한
   - 태양전지 쪽의 angle-restricted light trapping, LED·조명 쪽의 étendue-limited extraction, nonimaging optics(Winston, Miñano 등)의 phase-space 상한
4. OLED 외부 광추출에서 **"MLA가 되돌린 빛의 각도를 보존한다", "통과마다 탈출확률이 감소한다"**를 정량화한 선행 연구
   - Nature 원고(Kim et al.)가 인용한 Greiner 2007, Levell 2013(SPIE) 등의 device reflectivity 연구 포함

### B. 확인 권장
5. 수동 구조(MLA, 산란층, 텍스처, 격자) 전반에 대해 **공통 상한이나 통일된 분류 틀**을 제시한 리뷰나 이론 논문 (OLED, LED, 태양전지)
6. 광추출 구조의 출광 각도 분포가 **구조와 무관하게 Lambertian 근처로 수렴한다**는 관찰이나 이론
7. ML(DQN, BO 등) 기반 광추출 구조 최적화 논문 중 **물리적 상한 대비 위치**를 함께 제시한 것

### C. 용어 충돌 확인
8. "angular memory", "memory effect"는 산란 광학의 optical memory effect(Feng, Kane, Lee, Stone 1988; Freund 1988; 스페클 기울기 상관)와 충돌합니다. "angular mixing number", "angular persistence", "angular relaxation time"이 다른 분야에서 이미 다른 뜻으로 쓰이는지 확인해 주세요.

## 검색어 (영문, 조합해서 사용)
- "redistribution matrix" eigenvalue light trapping
- "angular redistribution" Markov chain solar cell / LED / OLED
- "mixing time" OR "spectral gap" light trapping OR outcoupling
- "degree of randomization" light trapping
- ergodic light trapping chaotic cavity
- Lambertian scatterer "partial randomization" light trapping limit
- generalized Yablonovitch limit anisotropic / angle-restricted
- étendue limit light extraction LED "radiance" bound phase space
- OLED substrate mode recycling microlens "angle" preserved reflection
- outcoupling structure classification microlens scattering layer unified model

## 조사 방법
- Google Scholar, Web of Science 또는 Scopus, arXiv, Optica Publishing, ACS, Wiley, Nature 계열을 사용합니다.
- 가장 가까운 논문 5–10편은 **원문(본문, 부록)을 직접 확인**합니다. 초록만 보고 판단하지 않습니다.
- 가까운 논문의 인용 목록(backward)과 피인용 목록(forward)을 한 단계씩 추적합니다.
- 2026년 11월 현재까지의 최신 논문을 포함합니다.

## 결과물 형식
`literature/angular_mixing_survey.md`로 저장하고 커밋해 주세요.

1. **결론 요약 (5줄 이내)**
   - 일반화 상한과 각도 혼합 지표 각각에 대해 "선행 있음 / 부분 있음 / 없음" 판정과 근거
2. **가장 가까운 선행 연구 표**

| 순위 | 서지 (저자, 연도, 저널, DOI) | 무엇을 했나 | 우리 정의와 같은 점 | 다른 점 | 인용 필요성 (필수/권장/참고) |
|---|---|---|---|---|---|

3. **용어 충돌 확인 결과**: 후보 이름별로 다른 분야에서의 기존 사용 여부
4. **차별화 문장 초안 (영문 2–3문장)**: 논문 서론에 쓸 "기존 연구와의 차이"
5. **원문을 확인하지 못한 논문 목록**: 접근 실패 사유 포함

## 주의
- 수치나 주장은 원문에서 확인한 것만 표에 적고, 초록이나 2차 자료 기반이면 표시해 주세요.
- 판단이 애매한 논문은 빼지 말고 "애매함"으로 표시해 포함해 주세요.

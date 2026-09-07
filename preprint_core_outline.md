# 코어 프리프린트 아웃라인 (arXiv v1, Letter 분량)

작성일: 2026-09-07 · 목적: UDC OLEDX 대비 **과학적 우선권 확보** — 원리 + 첫 정량 시연만. 나머지는 같은 arXiv 엔트리의 v2(본지 투고 원고)로.

## 0. 한 줄 claim (이 문장에 필요한 것만 넣는다)

> 아웃커플링 구조에 **(근)무손실 후면 반사체**를 결합하면 광자 재활용이 추출효율을 1로 몰아간다. 이를 닫힌형 법칙으로 유도하고, 그 순간 **내부 기생흡수가 지배 손실**이 됨을 보이며, 실험으로 **DBR 반사체 EQE 77%**(기록)와 **실측 추출효율 91.6%**(Ag 반사체, 모델과 0.1%p 일치)를 시연한다.

분량 목표: 본문 1,500–2,000 단어, 그림 3개(최소 2개), Methods 1쪽, SI 5쪽 이내. 준비 2–3주.

## 1. 제목 후보 (원리 우선, 소자 숫자는 제목에 안 씀)

1. *Lossless rear reflectors enable near-unity light extraction in organic light-emitting diodes*
2. *A photon-recycling law for near-unity light extraction in OLEDs*
3. *Near-unity photon extraction in OLEDs by lossless-reflector photon recycling*

## 2. 그림 구성 (기존 슬라이드 → 패널 매핑)

### Fig. 1 — 원리와 법칙 (전부 기존 자료)
- (a) 개념도: 아웃커플링 구조 + 후면 반사체에서의 다중 반사 광 경로, Al(손실) vs DBR/Ag(근무손실). [MSCA deck 개념도 재작도]
- (b) **마스터 커브** η_ext vs R_LED(1−A′_para): "w/o outcoupling" / "w/ microcavity or horizontal emitter" / "+Al" / "+lossless reflector" 영역 표시. [슬라이드 8/10]
  - 본문 Eq.(1) 재귀합 → Eq.(2) 닫힌형 η_ext ≈ P̄·BSDF_T / (1−(1−P̄·BSDF_T)(1−A′_para)), η_out = η_ext·η_sub^(N_R=0).
- (c) **EQE vs n_sub**, Al 전극 vs IZO+DBR: Al은 ~0.60에서 포화 후 감소, DBR은 ~0.82까지 단조 상승, "previous works" 표시. [슬라이드 11] → "왜 기존 연구가 ~60%에서 멈췄는가"를 한 패널로.
- (SI로) 전극별 각도–파장 반사율 맵 [슬라이드 9], 손실 회계 스택 [슬라이드 10].

### Fig. 2 — 실험 (2세트, 각 1패널)
- (a) **DBR 보조 반사체 (Exp #1)**: 소자 단면 모식도 + EQE–휘도 곡선 Ref(~23%) / MLA+Al(~61%) / MLA+DBR(**~77%**). 본문에 배면 누설 18.6% → 0.9% 수치만 언급(표는 SI). [슬라이드 12, 38]
- (b) **Ag 보조 반사체 (Exp #2)**: 소자 단면 모식도 + EQE–휘도 곡선 Ref / MLA, 실측 **55.8%** vs 시뮬 **55.9%**, 실측 추출효율 **91.6%**(vs 기존 75.4%, Song 2018) 표기. [슬라이드 13]
- 반드시 명기: 두 실험 모두 반사체는 **전극 뒤·바깥을 덮는 보조(supplementary) 반사체**이며 전극 자체가 아님; Exp #1 발광체 Ir(dmppy-ph)₂tmd + n=1.77 MLA 기판(HNPS-01), Exp #2 발광체 Ir(ppy)₂acac + 유리/n≈1.5 MLA 필름; DBR은 열증착 ZnS/LiF 4.5스택; η_ext의 분모는 시뮬 값(분자는 실측).

### Fig. 3 — 설계 규칙의 재정의 (시뮬, 기존 자료; 분량 부담 시 (c)만 Fig.1(d)로 흡수하고 나머지 삭제)
- (a) EQE(d_ETL, d_HTL) 맵: w/o MLA(FP 간섭무늬) vs w/ MLA(d_ETL≳100 nm에서 평탄, 최대 87.6%) — "두께 최적화 불필요". [슬라이드 16]
- (b) EQE vs Θ: 제안 구조 87→89.5%(평탄) vs 기준 25→49%(급경사) — "수평 배향 불필요". [슬라이드 15]
- (c) Max EQE vs k_ITO(96→55%) / vs n_Ag(95→46%) — "아웃커플링이 붙으면 기생흡수가 지배 손실; 설계 규칙 재정의". [슬라이드 19]

## 3. 본문 구조 (문단 단위)

1. **Intro ¶1**: IQE≈100% vs EQE 20–30%; 손실 채널; 아웃커플링+고굴절로 ~60%까지 오고 정체.
2. **Intro ¶2**: 정체 원인 = 다중 통과 시 후면 반사체/TCO 흡수; 본 연구 thesis(무손실 반사체 + 재활용 법칙). *선택*: "최근 산업계 발표(UDC OLEDX, 2026)가 반사체 기반 재활용의 중요성을 보여주나 정량 설계이론·무손실 반사체 시연은 보고된 바 없다" 1문장 — 저자 판단(§6 참고).
3. **Theory ¶**: Eq.(1)–(2) 유도 요지, 마스터 커브 해석(R_LED→1 근방 비선형성), η_out = η_ext·η_sub 분해.
4. **Simulation ¶**: EQE vs n_sub(트레이드오프 해소), 기생흡수 지배 전환.
5. **Exp ¶1 (DBR)**: 구조, 77%, Al→DBR +16%p, 누설 억제 수치.
6. **Exp ¶2 (Ag)**: 구조, 55.8/55.9 일치, η_ext 91.6%.
7. **Design rules ¶**: Fig. 3 요약 — 두께·HDR 둔감, 기생흡수 민감 → 설계 규칙 재정의.
8. **Outlook ¶ (짧게)**: 모델 예측 — 전극 통합형 Ag n=1.77에서 86.6%, PLQY≈1 발광체로 >90%(텍스트만, 그림 없음); 보조 반사체의 공간 비효율 한계 1문장; 전극 통합·TEOLED는 후속.

## 4. Methods (압축) + SI

- Methods: 소자 제작(2세트 각각 스택 명시), DBR 제작(ZnS/LiF 열증착; 두께 최적화 알고리즘 — GA/PSO 중 실제 사용한 것으로 통일), EQE 측정(적분구, 측정 휘도, 에지 발광 처리, 소자 수 N), trans-scale 시뮬(2-step, BSDF, PSO), lossy-EML(k>0) 포함 명시.
- SI: Pout 행렬 전개(전체식), point-source 측정기하 무관성(EQE 74.1/74.0/74.0%), BSDF 행렬, 반사율 맵, 손실 회계, 누설 표(0/2.5/4.5 stacks), 광학상수(ITO/IZO n,k).
- Data/Code availability: GitHub 레포(TMF + BSDF 코드) 인용 — 오픈 모델 자체가 OLEDX 대비 차별점.

## 5. 넣지 않는 것 (v2/본지 원고 또는 별도 논문으로)

| 항목 | 제외 이유 |
|---|---|
| 폴라리톤 OLED 전부 | UDC PEP 로드맵과 겹침 — 선점 빌미 주지 말 것 |
| TEOLED 픽셀화, EHD 프린팅, 유연/KK 염료 | 데이터 미성숙; 별도 논문 계획 |
| seed layer 스크리닝, IZO 공정·trap filling·EDS | 공정 디테일 — 경쟁자에게 레시피만 주는 셈 |
| 얇은 Ag 두께 시리즈(5–25 nm), 전극 조합 7행 표 | 코어 claim에 불필요 |
| MLA 형상 허용도·산란층 S–g 맵(design rule #3) | 경계선 — 본문 1문장으로 "구조 무관성"만 언급, 그림은 v2 |
| 95.1% 사다리 그림 | 텍스트 1문장으로 충분; 비현실 가정이 v1에 그림으로 남는 것 회피 |
| DBR 2.5스택 곡선, 소자 사진, J-V-L, 롤오프 분석 | SI 또는 v2 |

## 6. 게시 전 체크리스트

1. **특허 결정** (EPO 유예 없음) → 출원 의사 있으면 arXiv 전에 우선권 출원.
2. **PI·공저자 동의** (Gather 랩 DBR 데이터, KAIST 협업 데이터).
3. **오류 수정 확인**: HDR 등방=0.67 표기, 86.6% PSO 여부, Cs₂CO₃ 표기, DBR 재료(ZnS/LiF) 및 최적화 알고리즘, EML 두께 통일.
4. **실측/시뮬 라벨링**: 모든 수치에 (exp.)/(sim.) 명시; η_ext 정의 명시.
5. EQE는 최대값 + 지정 휘도(예: 100, 1000 cd/m²)에서의 값 병기; 소자 수 N 명시.
6. OLEDX 언급 여부 결정 — 언급 시 중립적 1문장 + URL 인용.
7. arXiv 카테고리: physics.optics (cross-list physics.app-ph). 게시일 = 저널 투고일과 맞출 필요 없음(코어는 먼저).

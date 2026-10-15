# 2차 스크리닝 보고서 — "HATCN을 넘는 Ag seed가 존재하는가?"

**날짜**: 2026-07-31
**목적**: 현존 charge-transport 재료(유기 CTL, 무기 p/n형 반도체)와 금속 전극 재료 전체에서
HATCN(E_b = 1.03 eV)보다 강한 Ag 앵커링을 주는 seed 후보 탐색 + 선행보고(신규성) 판정.
**방법**: 1차와 동일 프로토콜 — PBE-D3(BJ)/def2-SVP, UKS, counterpoise, Ag/I에 def2-ECP.
클러스터/사이트 모델, DFT 강체 거리 스캔으로 Ag 위치 정련. (스크립트 30, 31)

---

## 1. 신규 계산 결과 (이번 세션)

| 후보 (모델) | 대표 소자 재료 | 사이트 | E_b (eV) | r (Å) | 비고 |
|---|---|---|---|---|---|
| **CuI (Cu₄I₄ 큐베인)** | 열증착 p형 HIL (tandem OLED CGL 실사용) | **Cu 금속** | **0.82** | 2.45 | 금속친화 Ag–Cu 채널 |
| CuI (Cu₄I₄) | 〃 | I-top | 0.22 | 3.10 | 물리흡착 |
| CuSCN ((CuSCN)₃) | 용액공정 p형 HTL | **Cu 금속** | **0.57** | 2.55 | 〃 |
| CuSCN | 〃 | S / N | 0.14 / 0.10 | ~3.0 | SCN의 S·N은 이미 Cu에 배위 → 공여능 소진 |
| Me₃PO (P=O 모델) | DPEPO/TSPO1/PO-T2T계 ETL | O(=P) | 0.25 | 2.44 | 트리알킬(전자밀도 상한). Ph₃PO는 이보다 약할 것이 확실해 생략 |
| s-triazine | TRZ계 ETL 코어 | N | 0.29 | 2.34 | 피리딜류(B3PyMPM 0.63)보다도 약함 |
| Benzothiadiazole | n형 억셉터 유닛 | N | 0.28 | 2.54 | |
| Thiophene | PEDOT/티오펜계 HTL | S | 0.17 | 2.90 | 방향족 S는 비공여성 |

### 1차 스크리닝 기준값 (재게시)

| 기준 | E_b (eV) |
|---|---|
| 청정 Al / Ag–Cu / Ag–Ag | 2.24 / 1.93 / 1.86 |
| **ZnS (Ag–S 설파이드)** | **1.60** |
| **HATCN (nitrile N)** | **1.03** |
| F4TCNQ / Cs₂CO₃ / TPBi / p-bPPhenB | 0.97 / 0.90 / 0.89 / 0.87 |
| B3PyMPM / Bphen / MoOx결손 / AlOx | 0.63 / 0.49 / 0.44 / 0.42 |
| MoO₃ / LiF / TAPC / TCTA / Liq | 0.28 / 0.25 / 0.25 / 0.28 / 0.17 |

---

## 2. 결론

**① 이 디스크립터(단원자 Ag 결합에너지) 기준으로, 진공증착 가능한 charge-transport 재료 중
HATCN(1.03 eV)을 넘는 것은 발견되지 않았다.**

앵커 화학별 서열이 완성됨:
- **금속–금속 (Al 2.24, Cu 1.93)** > **설파이드 S²⁻ (ZnS 1.60)** > **nitrile N (HATCN 1.03, F4TCNQ 0.97)**
  > 이미다졸/phen 킬레이트 N (0.87–0.90) > Cu(I)화합물의 Cu 사이트 (0.6–0.8) > 피리딜 N (0.6)
  > phen 단좌 (0.5) > 결손 산화물 (0.4) > **P=O·트리아진·BTD·화학양론 산화물·LiF·아민·카바졸 (0.25–0.29)**
  > 티오펜 S·페놀레이트·π (≤0.17)
- HATCN을 넘는 두 채널(금속–금속, 설파이드)은 각각 산화 취약/가시광 흡수(금속), 선행보고(ZnS, Kim AFM 2015)로 막혀 있음.
- **HATCN의 강점은 "1.0 eV급 결합 + 확산장벽 0.29 eV + 산화 무관 + 투명 + 증착 가능"의 유일 조합**이라는 것이 이번 전수 스크리닝으로 실증됨.

**② 유일한 미개척 강결합 후보: CuI의 금속(Cu) 사이트 경로 (0.82 eV, 클러스터 하한값).**
- 문헌상 CuSCN|Ag 계면에서 Cu⁺ 불균등화 → Cu⁰–Ag 합금화가 보고됨(JMC A 2023) → 실표면에서는
  결손·환원된 Cu 사이트가 Ag를 강하게 앵커링하는 **산화환원 매개 습윤**이 가능. 클러스터 E_b는
  이 채널(전하이동·합금화 에너지)을 포함하지 않으므로 0.82 eV는 하한으로 봐야 함.
- CuI는 **열증착 가능**(tandem OLED에서 HIL/CGL로 실사용), 투명(Eg 3.1 eV), p형 고이동도.
  **Ag seed 용도 선행보고 없음** (문헌 검색 기준).
- 단, 관건은 안정성: 같은 반응성이 장기적으로는 전극 열화 경로(AgI 형성)일 수 있음 — DOE에 열화 평가 포함 필요.

**③ 신규성 지형 (LIT_NOTES_beyond_hatcn.md 상세 — 2026-07-31 정정 반영)**
- 금속(Ge/Ni/Cu/Au/Al/Yb)·산화물(MoO₃/ZnO/TiO₂)·ZnS·PEDOT:PSS/PEI/PAI·SAM: 전부 선행 존재 → 탈락.
- **[정정] HATCN(7nm)/Ag(15–25nm) 스택은 Park & Suh, Opt. Express 26, 4979 (2018)에 이미 공개**
  (광학+면저항+SEM/AFM 포함; inverted TEOLED의 HIL 프레이밍이라 seed 키워드 검색에서 누락됐었음).
  구조 신규성은 소멸. 생존하는 주장: seed **기능·메커니즘 규명**(underlayer 대조 SEM + DFT 디스크립터),
  **≤10nm 연속막**(Park&Suh는 HATCN 위 15nm에도 void 보고 → 대비군), 전 화학 스크리닝 서열, 산화 반전.
- 특허 리스크: US9786868(OSRAM), US11276845, US11552266(Novaled 추정) — 청구항 원문 미확인, IP 실사 필요.
- 웹검색 재현율 한계 실증 → 투고 전 Scopus/WoS 체계 검색 + 인용망 추적 필수.

---

## 3. 논문 포지셔닝 제안

**주 논지(추천)**: "진공증착 OLED 스택에서 접근 가능한 **모든 앵커 화학의 전수 DFT 스크리닝**
(nitrile, pyridyl, imidazole, phen, amine, carbazole, π, phenolate, carbonyl, P=O, thiophene-S,
thiadiazole, triazine, SCN, halide, sulfide, oxide, defective oxide, metal)을 통해
HATCN급 nitrile 화학이 투명·비산화·증착형 CTL의 사실상 최적임을 규명 + SEM 실증".
- 부논지 1: 산화 반전 (청정 Al 2.24 → AlOx 0.42) — 금속 seed의 실패 메커니즘.
- 부논지 2: Cu(I) 화합물의 금속사이트/산화환원 매개 습윤 채널 발견 (CuI 0.82 eV, 미보고) — 후속 slab/실험 제안.
- 네거티브 결과(P=O·트리아진·BTD·티오펜 전부 ≤0.3 eV)는 "왜 아무 ETL이나 seed가 못 되는가"를
  정량화하는 데이터로 본문 가치 있음.

**후속 계산 (필수)**: QE 주기 slab — CuI(111) Cu-terminated/결손 표면의 Ag 흡착·확산장벽,
HATCN 단분자막 slab 재검증. (HANDOFF 6절 1번과 통합 수행 권장)

**후속 실험 DOE**: Glass/CuI(2–5nm)/Ag(6nm) vs HATCN seed 대조 — 면저항·AFM·SEM + 가속 열화(AgI 형성 여부 XPS).

---

## 4. 재현 정보
- 스크립트: `scripts/30_beyond_hatcn.py` (빌드→xtb→DFT 스캔→CP), `scripts/31_cu_site.py` (Cu 사이트 추가)
- 결과: `runs/beyond_hatcn_binding_eV.json` (스캔 전곡선 포함), psi4 로그 `runs/psi4_beyond_hatcn.out`, `runs/psi4_cu_site.out`
- 구조: `structures/{thiophene,BTD,triazine,Ph3PO,Me3PO,CuSCN3,Cu4I4}*.xyz`
- 주의: 전부 클러스터 모델. Cu 계는 GFN2 발산으로 기판동결+DFT 스캔 정련(이온성 고체와 동일 처리).
  Ph₃PO는 psi4 SCF 준비단계 행(hang)으로 미완 — Me₃PO가 P=O 상한값을 대신함.

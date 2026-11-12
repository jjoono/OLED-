# 문헌 전수조사 지시서: 초박막 Ag seed layer / HATCN 기반 투명 상부전극

> 이 문서는 다른 세션(조사 담당)에 그대로 전달하는 작업 지시서입니다.
> 결과물은 아래 **5. 산출물** 형식을 정확히 따르십시오.

---

## 0. 배경 (조사자가 알아야 할 맥락)

투명/전면발광 OLED의 상부전극으로 열증착 초박막 Ag(5–10 nm)를 쓰려 한다.
Ag가 얇아도 연속막이 되게 하려면 아래에 seed layer가 필요하다. 연구팀의 계산
결과와 가설은 다음과 같다. 문헌으로 **지지되는지, 반박되는지, 아직 공백인지**
판정하는 것이 이번 조사의 목적이다.

1. Seed가 Ag를 잡는 세기(결합 에너지 E_b, 대략 0.5 eV 이상)가 핵 생성 밀도를
   정한다. HATCN(1.63 eV, 계산값)이 가장 강하고, ETL류(TPBi·Bphen·BCP 등)는
   0.4–0.8 eV, 정공수송 아민류(NPB·TAPC·TCTA)는 약 0.25 eV다.
2. Ag와 seed 사이 전하 이동의 **방향**이 중요하다. MoOx는 Ag를 강하게 잡지만,
   Ag 원자가 늘수록 Ag 쪽 양전하가 커져(Δq > 0) 연속막 형성에 불리하다는 가설이다.
3. HATCN의 약점은 결정화다. Park & Suh, *Opt. Express* **26**, 4979 (2018)은
   HATCN 7 nm 위의 Ag 15 nm에서 공극과 10.7 Ω/□을 보고했다. 연구팀은
   HATCN 5 nm에서 Ag 5–7 nm가 닫힌다. 두께와 결정화가 차이를 설명한다는 가설이다.
4. 결정화를 막으려고 HATCN에 다른 물질(TPBi 등)을 섞는 방안을 검토했다.
   하지만 실제 소자 n-ETL(Bphen:Cs₂CO₃) / HATCN / Ag에서 HATCN은 LUMO로
   전자를 나른다. 그래서 LUMO가 높은 물질을 섞으면 전기적으로 희석된다는 우려가 있다.
   대안은 "얇은 HATCN 단독(3–5 nm)" 또는 "HATCN / HATCN:혼합 이중층"이다.
5. 계산으로 HAT 유도체 701개를 스크리닝했다. Ag 결합을 키우는 축합 고리는
   흡수를 가시광 쪽으로 끌어와서, HATCN보다 투명한 유도체는 찾지 못했다.

---

## 1. 조사 범위와 핵심 질문

각 주제마다 아래 질문에 답할 수 있는 문헌을 **빠짐없이** 찾는다.
기간은 제한하지 않는다. 고전 문헌부터 2026년까지, 최신 문헌을 우선한다.

### A. 유기 seed layer 위 초박막 Ag 성장
- 유기물·무기물 seed(HATCN, MoO₃/MoOx, WO₃, ZnO, Cu, Ge, Al, Ca, PEI,
  티올, 질소 함유 고분자 등) 위에서 Ag의 퍼콜레이션 두께와 면저항, 투과율 데이터
- 증착 속도, 기판 온도, 표면에너지가 Ag 성장 모드(Volmer–Weber)에 미치는 영향
- 핵 생성 밀도, 섬 크기, 퍼콜레이션 두께를 정량 측정한 연구(SEM, AFM, in-situ 저항,
  GISAXS 등)

### B. HATCN 막 형태: 결정화, 두께 의존성, 열 안정성
- HATCN 박막의 결정화 시작 두께, 표면 거칠기, 열처리와 시간에 따른 변화
- HATCN 위 금속(Ag, Al, Au) 성장 연구. Park & Suh 2018을 인용한 후속 연구 전부 포함
- HATCN 혼합막(HATCN:NPB, HATCN:TAPC, HATCN:기타)의 형태와 결정화 억제 보고

### C. HATCN 기반 전하 생성층과 상부전극 구조
- n-도핑 ETL(Bphen:Cs₂CO₃, Bphen:Li, BCP:Cs, Liq 등) / HATCN 접합의 전하 생성
  메커니즘(계면 쌍극자, 준위 정렬, UPS/XPS)
- n-ETL / HATCN / Ag(또는 Al) 구조의 반투명 또는 상부 음극 소자 보고:
  구동전압, 효율, 투과율
- HATCN 두께에 따른 구동전압과 효율 변화

### D. 혼합·희석 억셉터층의 전기적 특성 (우려 4번 검증용)
- HATCN 또는 다른 억셉터(F4TCNQ, F6-TCNNQ, MoO₃)를 넓은 갭 물질로 희석했을 때
  전도도와 이동도 변화(희석률 의존성, 퍼콜레이션)
- 혼합 전하 생성층, 이중층·경사층 전하 생성층과 단일층 비교
- HATCN:TPBi 또는 유사한 "억셉터 + ETL 호스트" 조합의 직접 보고 (없으면 "없음"으로 명시)

### E. 금속/유기 계면 전하 이동 (가설 2번 검증용)
- Ag/HATCN, Ag/F4TCNQ, Ag/MoOx 계면의 XPS(Ag 3d, N 1s, Mo 3d) 화학 이동,
  일함수 변화, 금속 클러스터 크기에 따른 전하 상태
- 산화물 위 금속 클러스터의 전하 상태가 크기에 따라 달라지는 보고

### F. 박막 Ag 음극 투명·전면발광 OLED
- 유전체/금속/유전체 구조, 미세공진, 캡핑층을 포함한 박막 Ag 상부전극 소자 성능
- 결정화, 응집 등 seed 형태 문제로 생긴 소자 열화 보고

### G. 계산 방법 (방법론 논문용)
- 유기 분자 위 금속 원자·클러스터의 흡착 에너지 DFT 벤치마크(PBE0, D3, ωB97X 등)
- xTB/GFN2와 DFT 사이의 Δ-ML 보정, 능동 학습, 분자 계열별 보정
- 유기 기판 위 금속 박막 성장의 kMC 시뮬레이션
- 금속 seed나 접착층을 기계학습이나 하이스루풋으로 스크리닝한 선행 연구.
  **선행 사례가 있는지가 신규성 판단에 결정적이다.**

### H. HAT 유도체 합성과 소자 응용
- 헥사케토사이클로헥산 + 1,2-디아민 경로(HATCN, 벤조 축합 HAT, (CN)₂-벤조 HAT)
- 헥사아미노벤젠 + 1,2-디케톤 경로(아릴 치환 HAT: 페닐, C₆F₅, 4-CN-페닐 등)
- 위 유도체들의 EA/LUMO, 흡수 스펙트럼, 승화성, Tg 또는 결정성,
  OLED·OPV·전하 생성층 적용 보고

---

## 2. 검색 전략

- **데이터베이스:** Google Scholar, Web of Science 또는 Scopus(가능하면),
  출판사 사이트(ACS, Wiley, RSC, Optica, Elsevier, Nature, AIP, APS), arXiv/ChemRxiv
- **인용 추적:** 아래 앵커 논문들을 인용한 문헌(forward)과 이들이 인용한 문헌(backward)을
  모두 확인한다.
  - Park & Suh, *Opt. Express* 26, 4979 (2018)
  - 각 주제에서 처음 찾은 피인용 상위 3편
- **키워드 예시** (영어로 검색하고 조합해서 쓴다)
  - ultrathin silver / thin Ag electrode / percolation threshold / seed layer / wetting layer / nucleation layer
  - HAT-CN / HATCN / hexaazatriphenylene hexacarbonitrile + silver / crystallization / morphology / thickness
  - charge generation layer + HAT-CN + n-doped / Cs2CO3 / Bphen
  - transparent OLED / top-emitting OLED + thin Ag cathode
  - mixed / blended / diluted + HAT-CN / acceptor + conductivity / mobility
  - Ag cluster + charge transfer + XPS + organic acceptor / MoO3
  - adsorption energy + metal atom + organic molecule + DFT benchmark
  - delta machine learning + xTB / GFN2 + DFT
  - kinetic Monte Carlo + metal growth + organic substrate
  - hexaazatriphenylene derivatives + synthesis + electron affinity
- 한국어·중국어 학위논문과 특허(KIPRIS, Google Patents)도 HATCN seed와 상부전극
  관련은 확인한다.

---

## 3. 포함·제외 기준

- **포함:** 실험 수치(두께, 면저항, 투과율, 구동전압, 에너지 준위)가 있거나,
  메커니즘을 직접 다룬 1차 문헌. 핵심 리뷰는 별도 표기한다.
- **제외:** 수치나 메커니즘 없이 언급만 하는 문헌. 단, 앵커 논문을 인용한 문헌은
  제외하더라도 목록에 "제외 사유"와 함께 남긴다.
- 같은 그룹의 중복 결과는 대표 1편으로 묶고 나머지는 하위 항목으로 둔다.

---

## 4. 반드시 지킬 것

1. **모든 문헌은 DOI(또는 arXiv ID, 특허번호)로 실제 존재를 확인한다.**
   확인하지 못한 문헌은 표에 넣지 말고, 따로 "미확인" 목록에 둔다.
   **기억에 의존해 서지를 만들지 않는다.**
2. 수치는 원문 그림이나 표에서 읽은 값만 적고, 출처 위치를 표시한다(예: Fig. 3b, Table 1).
   단위를 반드시 적는다.
3. 우리 가설(0절 1–5번)과의 관계를 문헌마다 **지지 / 반박 / 무관 / 판단 불가** 중
   하나로 표시하고, 한 줄 근거를 단다.
4. 서로 모순되는 보고를 찾으면 별도 절에 모은다(예: HATCN 위 Ag가 잘 닫힌다는 보고 vs
   결정화 때문에 안 된다는 보고). 두께, 속도, 기판 같은 조건 차이를 함께 적는다.
5. 요약과 판정은 한국어로, 서지와 원문 인용은 영어 원문 그대로 쓴다.

---

## 5. 산출물

저장 위치: `dft_seedlayer_screen/notes/` (해당 저장소가 없으면 작업 디렉터리)

1. **`LIT_SURVEY_TABLE.csv`**: 문헌 1편당 1행. 열 구성:

   | 열 | 내용 |
   |---|---|
   | id | 일련번호 |
   | topic | A–H (복수 가능, 예: "A;B") |
   | authors_short | 제1저자 et al. |
   | year | |
   | journal | |
   | doi | 필수 |
   | seed_material | 예: HATCN, MoO3, HATCN:NPB |
   | seed_thickness_nm | |
   | metal | Ag/Al/Au 등 |
   | metal_thickness_nm | |
   | percolation_nm | 연속막이 되는 두께 |
   | sheet_res_ohm_sq | (두께와 함께) |
   | T_550_percent | 550 nm 투과율 |
   | dep_rate_A_s | 증착 속도 |
   | substrate | |
   | device_stack | 소자가 있으면 전체 구조 |
   | key_numbers | 기타 수치 (출처 위치 포함) |
   | key_finding | 한 줄 요약 |
   | hypothesis | 관련 가설 번호와 지지/반박/무관/판단 불가 |
   | note | 조건 차이, 의문점 |

2. **`LIT_SURVEY_SUMMARY.md`**: 아래 순서로 작성한다.
   1. 주제 A–H별 핵심 결론. 3–6줄로 쓰고, 근거가 되는 문헌 id를 단다.
   2. 0절 가설 1–5번 각각의 판정과 근거. 선행 연구가 이미 같은 결론을 냈는지도 적는다.
   3. 모순되는 보고 목록과 조건 차이
   4. **신규성 판단:** 다음 각각이 이미 보고되었는가?
      - "seed–Ag 전하 이동 방향(Δq)을 연속막 형성의 판별 기준으로 쓴 연구"
      - "계열별 DFT 앵커로 보정한 seed 스크리닝 도구"
      - "HATCN 결정화와 두께 창"
      - "HATCN 이중층 또는 혼합 seed"
   5. 논문에서 **반드시 인용해야 할 문헌** 15–25편. 우선순위와 인용 위치
      (서론, 방법, 결과 비교 등)를 함께 적는다.
   6. 실험 설계에 바로 쓸 수 있는 수치 목록.
      예: HATCN 결정화 시작 두께 보고값 범위, n-ETL/HATCN 계면 쌍극자 크기,
      희석 억셉터층의 전도도 감소율
   7. 문헌 공백. 아무도 측정하지 않았거나 보고하지 않은 것

3. **`LIT_SURVEY_UNVERIFIED.md`**: 존재를 확인하지 못한 후보 문헌과,
   제외한 문헌 및 그 사유

---

## 6. 진행 방식

- 먼저 주제별로 후보를 넓게 모은 다음(목표: 전체 150편 이상 검토, 표에는 관련 있는 것만)
  인용 추적으로 보강한다.
- 주제 하나를 끝낼 때마다 CSV에 추가하고 저장해서, 중간에 끊겨도 결과가 남게 한다.
- 끝나면 "검토한 문헌 수 / 표에 넣은 수 / 미확인 수"를 요약 맨 위에 적는다.

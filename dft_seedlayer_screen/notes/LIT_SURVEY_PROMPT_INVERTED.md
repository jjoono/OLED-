# 문헌 조사 보완 지시서: 역구조(상부 양극) 기준

> 1차 조사(`notes/lit_survey/`, 353편)를 **보완**하는 작업이다. 처음부터 다시 하지 않는다.
> 1차 조사 결과 파일 세 개를 먼저 읽고, 이미 있는 문헌(DOI 기준)은 다시 넣지 않는다.

## 0. 왜 보완하는가

1차 조사는 "n-ETL / HATCN / Ag 상부 음극"(정구조, CGL형)을 전제로 설계됐다.
그런데 연구실의 주력 구조는 **역구조 전면발광 OLED**다.

```
ITO(음극) / EIL / ETL / EML / HTL(NPB, TAPC 등) / HIL = seed (HATCN 또는 MoOx) / Ag 5–10 nm(상부 양극) / 캡핑층
```

이 구조에서 HIL은 정공 주입층이면서 동시에 Ag seed다. 1차 조사의 주제 A, B, E, G, H는 소자 구조와 상관없어서 그대로 쓴다.
아래 I1–I6만 새로 조사한다.

연구팀의 현재 주장(검증 대상)
1. HATCN이 MoOx보다 얇은 상부 Ag의 seed로 낫다.
   - 연구실 결과: HATCN 5 nm 위 Ag 5–7 nm가 닫힌다.
   - 계산: HATCN E_b 1.63 eV(다중 사이트 PBE0)다. 아래 HTL에서 전자를 받아 환원돼도 1.77 eV로 Ag를 계속 잡는다.
2. Park & Suh(Opt. Express 26, 4979 (2018), 역구조, NPB 100 / HATCN 7 / Ag 15 nm, Ag 0.2 Å/s)의 공극을 두 가지로 설명할 수 있다.
   - HATCN 결정화
   - 느린 Ag 증착 속도
   - 이 둘을 구분할 근거를 찾는다.
3. 결정화를 피하는 구조로 다음 경사 구조를 검토 중이다.
   - HTL / HTL:HATCN(p-도핑 영역) / 얇은 neat HATCN(Ag 쪽) / Ag

## 1. 새 조사 주제

### I1. 역구조 OLED의 상부 양극
- 역구조 OLED, TEOLED, 투명 OLED에서 다음 상부 양극을 쓴 소자를 모두 찾는다.
  - HIL: HATCN, MoO3/MoOx, WO3, V2O5, ReO3, F4TCNQ, F6-TCNNQ 박막
  - 상부 금속: Ag, Au, Al, Ag 합금
- 수치: HIL 두께, Ag 두께와 속도, Rs, T, 구동전압, 효율, 각도 특성
- **HATCN/Ag와 MoOx/Ag를 같은 논문에서 비교한 연구는 최우선으로 찾는다.**
- Kippelen 그룹 문헌 세 건은 반드시 본문을 확보한다(1차 조사에서 'PDF needed').
  - Najafabadi et al., Org. Electron. 2013, 10.1016/j.orgel.2013.02.024
  - Knauer et al., Org. Electron. 2013, 10.1016/j.orgel.2013.06.004
  - Najafabadi 박사학위논문, Georgia Tech 2014, hdl:1853/52990, 4.3절

### I2. 역구조 OPV와 페로브스카이트의 상부 양극 (MoO3/Ag, HATCN/Ag)
- 이 분야에는 "HTL / MoO3 / 얇은 Ag" 상부 양극의 정량 데이터가 OLED보다 많다.
  - MoO3 두께에 따른 Ag 퍼콜레이션, Rs, T
  - 산화 상태
  - 열화(Ag 확산, MoO3 환원)
- HATCN을 MoO3 대신 쓴 역구조 OPV와 페로브스카이트 사례도 찾는다.

### I3. HTL/HATCN 계면과 혼합
- 아민 HTL(NPB, TAPC, TCTA, spiro 등)과 HATCN 사이의 전하 이동, 준위 정렬, 쌍극자(UPS/XPS)
- 아민 HTL 위 HATCN의 막 형태: 결정성, 거칠기, 두께 의존성. 1차 조사의 기판은 SiOx, ITO, Bphen이었다.
- HTL:HATCN 혼합층과 경사층(p-도핑 HTL)
  - 전도도, 투과율
  - CT 흡수 스펙트럼(파장, 흡광계수)
  - 혼합비 의존성
- "혼합층 위 neat HATCN 얇은 층" 같은 이중층이나 경사층 보고

### I4. 유기층 위 상부 금속 증착의 손상과 확산
- 상부에서 증착한 Ag, Au, Al 원자가 유기층(특히 HATCN, MoO3, HTL)으로 침투하는 깊이
- 계면 화학 반응, 열 손상
- Ag 침투가 정공 주입이나 누설에 미치는 영향

### I5. 정공 전용 소자
- ITO 또는 Ag/HTL/HIL/Ag 정공 전용 소자의 HIL 두께 의존성(HATCN, MoOx). 특히 HIL 1–10 nm 구간
- HIL이 얇아질 때 주입이 무너지는 두께(피복률 문제)

### I6. 상부 Ag 양극의 안정성
- 85 °C 또는 장시간 구동 뒤 얇은 Ag 상부 양극의 응집, Rs 변화, 캡핑 효과
- HIL 종류(HATCN 대 MoOx)에 따른 차이가 있다면 우선한다.

## 2. 규칙

1차 지시서(`notes/LIT_SURVEY_PROMPT.md`)의 다음 규칙을 그대로 따른다.
- 4절: DOI로 확인, 기억에 의존한 서지 금지, 수치 출처 위치, 가설 판정
- 3절: 포함·제외 기준

새 가설 번호(I-H1–I-H3)는 0절의 주장 1–3에 대응한다.

## 3. 산출물

저장 위치: `dft_seedlayer_screen/notes/lit_survey/`

1. `LIT_SURVEY_TABLE_INVERTED.csv`
   - 1차 표와 같은 열 구성에 두 열을 더한다.
     - `device_polarity`: inverted / conventional / n.a.
     - `hil_role`: HIL / CGL / seed only / n.a.
   - id는 I001부터 시작한다.
   - 1차 표 문헌 중 역구조에도 중요한 것은 다시 넣지 말고, `LIT_SURVEY_SUMMARY_INVERTED.md`에서 1차 id로 참조한다.
2. `LIT_SURVEY_SUMMARY_INVERTED.md`
   1. I1–I6 주제별 결론
   2. 0절 주장 1–3의 판정
   3. **HATCN/Ag 대 MoOx/Ag 직접 비교 데이터 표**(같은 논문 안 비교, 다른 논문끼리의 비교를 나눠서)
   4. 역구조 기준 신규성 판단
      - Ag ≤10 nm에서 HATCN 두께 창
      - 경사 seed
      - MoOx 대 HATCN 메커니즘
   5. 역구조 논문에서 반드시 인용할 문헌. 1차 목록에 추가할 것만 적는다.
   6. 실험 설계 수치. HTL 위 HATCN 두께, Ag 속도, CT 흡수, 정공 전용 소자 기준값
   7. 문헌 공백
3. `LIT_SURVEY_UNVERIFIED_INVERTED.md`: 존재를 확인하지 못한 후보와 제외한 문헌, 제외 사유

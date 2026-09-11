# 문헌 조사 노트 — "HATCN보다 나은 Ag seed" 신규성 판정 (2026-07-31)

## A. 이미 보고된 seed/wetting layer (신규성 없음 — 논문 후보 탈락)

| 재료 | 근거 | 비고 |
|---|---|---|
| 금속: Ge, Ni, Cu, Au, Cr, Ti, Al, NiCr | Chen 2010 (Ge), Formica (Cu), 다수 리뷰 | 가시광 흡수 문제. Al seed: Nanomaterials 12, 3540 (2022) |
| Yb (Yb:Ag 합금·계면층) | Light Sci Appl 2024 (ZnO/Yb:Ag), 산업계 표준 | 합금 접근, seed는 ZnO |
| 산화물: MoO₃, ZnO/AZO, TiO₂, SnO₂, ITO(IMI), NiO계 저방사 스택 | MoO₃/Ag (PMC6071051), 저방사 유리 산업 | ZnO는 산업 표준 Ag seed |
| ZnS | Kim AFM 2015 (ZnS/Ag/ZnS, Cs₂CO₃/Ag) | 우리 DFT Ag–S 1.60 eV와 부합. 선행 있음 |
| 유기 고분자: PEDOT:PSS, PEI, PAI | ACS AMI 2012 (PEDOT:PSS 위 초평활 Ag), PEI nucleation inducer | thiophene-S·amine 화학 기지 |
| SAM: MPTMS(thiol), APTMS(amine) | 고전 문헌 | charge transport 아님 |
| a-C:H, Cs₂CO₃ | 리뷰 언급, Kim 2015 | |
| 질소 표면처리(N deployment) | ACS Appl Nano Mater 2020 | 분자 아닌 공정 |

## B. 미보고(신규성 유력) — charge transport 가능 재료

| 후보 | 역할 | Ag-seed 선행보고 | 위험/비고 |
|---|---|---|---|
| **CuSCN** | p형 투명 HTL (용액공정) | **없음** (2회 검색) | Ag/CuSCN 계면반응 보고 있음(JMCA 2023, PSC 상부전극 열화 맥락 — Cu⁺ 불균등화·합금화). 반응성 계면=습윤에 유리할 수 있으나 안정성 검토 필요. 용액공정이라 하부(글라스측) seed로만 적합 |
| **CuI** | p형 투명 HIL/HTL, **열증착 가능**, tandem OLED CGL로 실사용 (J Inf Disp 2023) | **없음** | 공정 호환성 최고(진공 열증착). Ag–I 결합. AgI 형성 가능성(계면 반응) 확인 필요 |
| **포스핀옥사이드 ETL (DPEPO/TSPO1/PO-T2T)** | ETL/호스트, 열증착 | **없음** (Ag–O=P 배위착물 문헌만 존재) | E_b가 관건 (DFT로 판정) |
| **Benzothiadiazole(BTD) 유닛** | n형 억셉터 빌딩블록 | **없음** | 증착형 BTD 소분자 흔치 않음(주로 폴리머) |
| F6-TCNNQ | p-도판트 | 직접 보고 없음 | 그러나 nitrile 화학=HATCN과 동류(F4TCNQ 0.97<1.03), 성능 우위 없을 듯 |

## C. HATCN 자체의 Ag-seed 선행보고 여부 — ⚠ 2026-07-31 정정

**[정정] HATCN 위에 Ag를 증착한 구조는 이미 공개되어 있음.** 사용자 제공 논문으로 확인:
- **Park & Suh, Opt. Express 26(4), 4979 (2018)** (경희대): inverted TEOLED에서 **HATCN(7nm)/Ag(15·20·25nm)** 스택을 제작, 투과/반사/흡수 + 면저항(15nm: 10.7Ω/□) + **SEM/AFM 모폴로지**까지 측정. 15nm Ag에 void/crack 관찰 → 25nm 권장이 결론. HATCN의 역할은 HIL(p-n-p-n 접합)로 개념화. 관련 선행: 같은 그룹 Synth. Met. 162, 402 (2012) (HATCN 저전압 TEOLED), ZnS/Ag/MoO₃ "quasi-perfect Ag" Org. Electron. 14, 3437 (2013).
- 일반화: **inverted TEOLED/tandem에서 Ag 상부양극을 HATCN(HIL/CGL) 위에 올리는 것은 보편적 구조** — "seed" 프레이밍이 아니라 "주입층" 프레이밍이라 seed/wetting/nucleation 키워드 검색에 안 걸렸음 (본 조사의 방법론적 한계).

**따라서 신규성 재판정:**
- ✗ 죽음: "HATCN/Ag 스택" 구조 자체, "HATCN 위 Ag 광학·면저항 측정" 수준의 주장.
- ○ 생존(방어 가능): ① HATCN의 **seed 기능 규명** — 동일 조건에서 underlayer만 바꾼 대조실험(Al/AlOx vs HATCN)로 성장모드 차이 실증 + DFT 디스크립터(E_b·확산장벽)로 메커니즘 설명. ② **≤10nm 초박 영역**에서 연속막 실현 (Park&Suh는 HATCN 위에서도 15nm에 void 보고 → 오히려 대비군으로 인용 가치). ③ 전 앵커화학 전수 스크리닝 서열. ④ 산화 반전(Al 2.24→AlOx 0.42) 메커니즘.
- 논문 프레이밍은 "HATCN을 처음 깔았다"가 아니라 **"왜/얼마나 얇게 되는가 + 이보다 나은 것이 존재하는가"**로 잡아야 함. Park & Suh 2018 필수 인용.
- 특허 지뢰: **US9786868**(OSRAM/Siemens 계열, "metal growth layer"), **US11276845**("OLED with silver contacts"), **US11552266**(전극+유기반도체층, Novaled 추정), EP2750213. 웹 페치 차단으로 청구항 원문 미확인 — **IP 실사 필요**.
- ⚠ 본 웹검색(스니펫 기반)의 재현율 한계가 실증됨 — 투고 전 **Scopus/WoS/Scholar 체계 검색 + 인용망 추적**(특히 Suh/Kwon/Leo 그룹 TEOLED 문헌군) 필수.

## D. 계산 스크리닝(DFT) 신규성

- "유기 CTL 분자 위 Ag adatom 결합에너지 체계적 DFT 스크리닝 → seed 선별" 형식의 선행 논문 미확인. Ge–Ag 결합 논거(Ag–Ge>Ag–Ag) 수준의 단편적 언급만 존재.
- → 방법론(디스크립터 E_b×확산장벽 스크리닝) + 신규 재료(CuI/CuSCN/PO계) 조합이면 신규성 주장 가능.

## E. 검색에 걸린 2차 참고

- 인버티드 OLED에서 phen계 ETL에 Ag 도핑(nucleophilic ETL) — Org Electron 2020 (S1566119920304201): phen–Ag 상호작용 활용 선례(단, seed 아님).
- APL 118, 213301 (2021): ALD half-reaction으로 만든 nucleation induction layer(유연 OLED 상부 금속 양극) — 무기 ALD 접근.
- Nat Commun 11, 2874 (2020) 계열: Cu-doped Ag, "T>100%" 초박 Ag 전극 벤치마크.

## 종합 우선순위 (신규성 × 예상성능 × 공정성)

1. **CuI** — 열증착 p형, tandem OLED 실사용, Ag seed 미보고. DFT E_b 확인 중.
2. **CuSCN** — S 사이트 강결합 예상(ZnS 1.60 eV 유추), 미보고. 용액공정 한계.
3. **포스핀옥사이드 ETL** — 미보고, DFT 판정 대기.
4. BTD — 미보고이나 증착형 소분자 부재로 실용성 낮음.
5. (탈락) 금속류·산화물·ZnS·PEDOT/PEI/PAI — 전부 선행 존재.

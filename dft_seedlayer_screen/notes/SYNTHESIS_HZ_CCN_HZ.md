# Hz-C(CN)₂-Hz 합성 경로 검토 (2026-11-26)

목표: 2,2-bis(5,8-dicyano-s-heptazin-2-yl)malononitrile, C₁₉N₂₀, MW 508. 비교 대상: Hz-NCN-Hz = N,N-bis(5,8-dicyano-s-heptazin-2-yl)cyanamide.

## 1. 문헌 현황

- **heptazine에 탄소 친핵체를 붙인 보고는 찾지 못했다.** 알려진 heptazine 유도체는 거의 모두 cyameluric chloride(C₆N₇Cl₃)의 Cl을 N·O·S·P 친핵체로 치환해 만든다. C–C 결합은 Friedel–Crafts 아릴화 예가 하나 있을 뿐이고, 리뷰도 아릴 치환 heptazine이 드물다고 적는다.
- **같은 탄소에 전자 부족 아릴을 두 번 붙이는 SNAr은 어렵다.** 아릴아세토나이트릴에서 첫 생성물이 더 산성이어서 탈양성자화되면 두 번째 SNAr이 막히고, Pd 촉매가 필요했다(OBC 2011). malononitrile에서는 첫 생성물 [Ar–C(CN)₂]⁻가 tricyanomethanide처럼 전하가 퍼진 약한 친핵체가 된다.
- **이런 음이온은 반응 자리가 둘이다.** tricyanomethanide는 니트릴 질소로 공격해 케텐이민을 만드는 예가 보고돼 있고(PMC7318174), 실릴화가 C/N 선택성을 바꾼다(Angew. 2024).

## 2. 제안 경로 (추정, 실험 검증 필요)

| 단계 | 반응 | 위험 |
|---|---|---|
| 1 | Cl₃Hz + Na[CH(CN)₂] 1당량, MeCN 또는 THF, −30 → 0 °C → Na[Cl₂Hz–C(CN)₂] | 낮음–중간. 안정한 음이온 염으로 분리할 수 있을 것으로 본다. |
| 2 | Na[Cl₂Hz–C(CN)₂] + Cl₃Hz (가열, 극성 비양성자 용매, Ag⁺로 Cl 활성화 또는 실릴화 음이온 사용) → Cl₂Hz–C(CN)₂–HzCl₂ | **높음.** 약한 친핵체, C 공격 대 N 공격(케텐이민) 경쟁. 아래 계산 결과 참조. |
| 3 | Cl 4개 → CN (TMSCN/촉매, CuCN, 또는 Pd/Zn(CN)₂) | **높음.** CN⁻이 heptazine 고리 탄소를 공격해 [Hz–C(CN)₂]⁻를 이탈기로 내보내 다리가 끊어질 수 있다. |

대안
- **CN 먼저:** Cl₃Hz → Cl(NC)₂Hz(부분 시안화) → 다리 형성. 3단계의 다리 절단은 피하지만, heptazine 위 CN 자체가 탄소 음이온에 치환될 수 있다(tricyanotriazine의 CN은 이탈기다).
- **Pd 촉매 2차 아릴화:** 말로노나이트릴의 Pd α-아릴화는 알려져 있으나, 이미 아릴화된 [Ar–C(CN)₂]⁻의 2차 아릴화와 heptazine 염화물의 Pd 호환성은 미지수다.

## 3. 계산 점검 (scripts 183, 184)

GFN2-xTB/ALPB(MeCN) 구조, PBE0-D3(BJ)/def2-SVP + ddCOSMO(MeCN) 단일점. 반응 에너지만(장벽·엔트로피 제외), 배좌 1개.

| 비교 | PBE0 (eV) | GFN2 (eV) | 의미 |
|---|---|---|---|
| 케텐이민 Hz–C(CN)=C=N–Hz − 목표 Hz–C(CN)₂–Hz | **−0.17** | −0.62 | 질소로 붙은 이성질체가 더 안정 |
| 카보다이이미드 Hz–N=C=N–Hz − Hz–N(CN)–Hz | **−0.19** | −0.64 | Hz-NCN-Hz도 같은 문제 |
| Hz–C(CN)₂–Hz + CN⁻ → Hz(CN)₃ + [Hz–C(CN)₂]⁻ | −2.24 | −4.32 | 다리 절단이 크게 유리 |
| Hz–N(CN)–Hz + CN⁻ → Hz(CN)₃ + [Hz–NCN]⁻ | −1.48 | −3.60 | NCN 다리가 0.76 eV 덜 끊어짐 |

해석
- 두 목표 분자 모두 다리가 **열역학적으로 덜 안정한 이성질체**다. 차이(0.2 eV)가 작아 반응 조건(속도론)으로 목표 이성질체를 얻을 여지는 있다. 그러나 가열이 필요한 2단계에서는 질소 결합 이성질체가 섞이거나 주생성물이 될 가능성을 고려해야 한다.
- 다리 절단 에너지는 작은 CN⁻의 연속체 용매화가 부족해 과장됐을 수 있다. 다만 **C(CN)₂ 다리가 NCN 다리보다 0.76 eV 더 잘 끊긴다**는 상대 비교는 의미가 있다. [Hz–C(CN)₂]⁻가 매우 안정한 이탈기라는 뜻이고, 3단계 시안화와 완성 분자의 친핵체(물, 아민) 내성 모두에 불리하다.
- 염소 중간체 기준의 2·3단계 경쟁(C 대 N 공격, Cl→CN 치환 대 다리 절단)은 계산 중이다.

## 4. 판단

- Hz-C(CN)₂-Hz는 계산 성능은 1위지만 **합성 난이도는 후보 중 가장 높다.** 2단계(같은 탄소의 2차 SNAr)와 3단계(다리 절단) 모두 선례가 없고 계산상으로도 불리하다.
- 실제 반응에서 생길 가능성이 큰 **케텐이민·카보다이이미드 이성질체**도 seed 후보로 계산해 볼 가치가 있다. 두 분자는 다리가 선형·평면에 가까워 결정화 위험이 올라갈 수 있다.
- 실험은 값싼 cyanuric chloride로 1·2단계 모델 반응을 먼저 해 C/N 선택성과 다리 안정성을 확인하는 것이 순서다.

## 출처

- Frontiers in Chemistry 2021, heptazine chemistry: https://www.frontiersin.org/journals/chemistry/articles/10.3389/fchem.2021.717569/full
- Emergence of s-Heptazines (review): https://www.researchgate.net/publication/326868607
- Synthesis of some heptazine derivatives (Friedel–Crafts): https://www.researchgate.net/publication/338148080
- Renewing accessible heptazine chemistry: https://pmc.ncbi.nlm.nih.gov/articles/PMC6549597/
- Pd-catalysed electron-deficient aryl acetonitriles (OBC 2011): https://pubs.rsc.org/en/content/articlehtml/2011/ob/c0ob00903b
- Pd-catalysed dicyanomethyl polyphenyls with sodium malononitrile: https://www.thieme-connect.com/products/ejournals/pdf/10.1055/s-1995-3984.pdf
- Tricyanomethanide as N-nucleophile: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7318174/
- Tricyanomethane or dicyanoketenimine (Angew. 2024): https://onlinelibrary.wiley.com/doi/full/10.1002/anie.202413565
- EP1854797A1 (cyameluric chloride): https://patents.google.com/patent/EP1854797A1/en

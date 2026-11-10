# 문헌 조사 결과: 각도 혼합 지표와 일반화 radiance 상한

- **조사일:** 2026-11-10 (2026년 11월까지 발표된 문헌 포함)
- **지침:** `literature/SEARCH_BRIEF_angular_mixing.md`
- **검색 경로:**
  - 웹 검색
  - OpenAlex, Crossref, Semantic Scholar, Europe PMC, arXiv API
  - 출판사 사이트(기관 접근이 되는 곳)
  - Wayback 보관본(Optica 봇 확인 페이지 대신 사용)
- **인용 추적:** 가까운 논문은 인용 목록과 피인용 목록을 한 단계씩 추적했습니다. 대상은 Bortz 1997, Muschaweck 2013, Lenef 2018, David 2013, Greiner 2007, Levell 2013, Shiang 2004, Kovačič 2021·2022, OPTOS 2015·2016, Li 2015, Manzoor 2020, Battaglia 2012, Altmann 2013 등입니다.
- **근거 표기:**
  - **[원문]** 이번 조사에서 본문·부록을 직접 확인
  - **[원문\*]** 보조 조사에서 원문을 확인했고, 핵심 문장은 다시 확인하지 않음
  - **[초록]** 초록·서지 정보만 확인
  - **[2차]** 다른 문헌의 기술에 의존
- **수치 표기:** 표의 수치는 표시한 근거 수준에서 확인한 것입니다. [초록]·[2차] 수치는 원문에서 확인하지 못한 값입니다.

---

## 1. 결론 요약

1. **일반화 radiance 상한: 선행 있음(원리).** radiance가 큰 광선부터 출력 위상공간을 채우는 상한은 Bortz·Shatz·Ries(1997)가 일반형으로 제시했고, Muschaweck·Wiesmann(2013)과 Lenef 등(2014, 2018)이 비Lambertian 내부 radiance를 가진 LED·형광체의 단일 통과 추출에 적용했으며, Lenef(2018)는 제한된 출력 각도(étendue)까지 다룹니다.
2. **일반화 상한에서 남는 신규성:** OLED 기판 radiance를 polar band로 나눈 적용, band 사이의 결합 제약, 실제 MLA·산란층·평면 구조를 이 상한으로 평가한 것의 조합은 찾지 못했습니다.
3. **각도 혼합 지표: 지표로서는 없음, 개념은 부분 있음.** 통과별 각도 재분배 행렬은 OPTOS·AMF·GenPro4·RayFlare·David(2013)·Kim 등(2026, 자기 프리프린트)이 이미 쓰지만, λ₂나 1/(1−λ₂)를 구조 지표로 쓴 예는 태양전지, LED·OLED·백라이트, 광학 상한의 세 갈래 조사에서 모두 찾지 못했습니다.
4. **가장 가까운 개념:** Manzoor(2020)의 통과별 Lambertian 편차(2 pass 대 8 pass), Battaglia(2012)의 미해결 과제 명시, David(2013)의 "bounce마다 ergodicity 감소", 복사 전달의 1/(1−g)이며, 1/(1−λ₂)는 마르코프 연쇄의 relaxation time(Levin & Peres)에 해당하므로 이 관계를 본문에 밝혀야 합니다.
5. **이름:** "memory"는 optical memory effect(Feng 1988, Freund 1988)와 충돌하므로 피해야 하고, "angular mixing number"는 기존 사용례를 찾지 못했습니다.

---

## 2. 가장 가까운 선행 연구

### 2-1. 일반화 radiance 상한 (주장 1)

| 순위 | 서지 (저자, 연도, 저널, DOI) | 무엇을 했나 | 우리 정의와 같은 점 | 다른 점 | 인용 필요성 |
|---|---|---|---|---|---|
| 1 | Lenef, Piquette, Kelso, *ECS J. Solid State Sci. Technol.* 7, R3211 (2018), doi:10.1149/2.0271801jss **[원문]** | 2법칙 논증으로 표면 추출 상한을 유도함. Lambertian 입력이면 어떤 표면이든 투과는 (n₂/n₁)² 이하. 임의 radiance 입력이면 radiance가 큰 영역은 투과 1, 나머지는 0인 이진 투과함수가 최대(Eqs. 12–15). 출력이 제한된 각도 범위, 더 일반적으로 제한된 étendue인 경우로 확장. YAG→공기 예시: Lambertian 0.299, 등방 체적 발광은 평면 0.163, 최적 표면 0.546 | 단일 통과. 비등방 radiance를 밝은 영역부터 채움. 균일 radiance면 (n₂/n₁)²로 환원. 기하·회절·산란 표면 모두에 성립 | LED·형광체 대상. polar band 분할과 band 간 결합 제약 없음. 실제 구조를 이 상한 대비 평가하지 않음 | 필수 |
| 2 | Muschaweck, Wiesmann, *Adv. Opt. Technol.* 2, 291 (2013), doi:10.1515/aot-2013-0027 **[원문]** | "직접(단일 통과) 추출할 수 있는 최대 분율"을 두 그릇 비유로 유도함. 작은 출력 위상공간을 내부 radiance가 큰 부분부터 채움. 균일 radiance면 étendue 비. 투명 InGaN 접합(L ∝ 1/cosθ, n_i = 2.4, n_o = 1.41)에서 59%, Lambertian 접합은 35% | radiance 정렬 단일 통과 상한. 균일이면 étendue 비로 환원 | 반구 전체 출력만 다룸. band 분해 없음. LED 대상 | 필수 |
| 3 | Bortz, Shatz, Ries, *Proc. SPIE* 3139, 59 (1997), doi:10.1117/12.290222 **[원문]** | 비균일 광원·표적에 대해 "cumulative sorted source étendue"를 정의함. 가장 큰 radiance가 표적 위상공간의 가중치 큰 영역부터 채울 때 전달량이 최대임을 보임(Eqs. 2.2–2.6). 회전대칭계는 skewness별로 같은 구성 | radiance 정렬 상한의 일반형. 표적 영역이 임의이므로 band 상한도 이 틀에 들어감 | 조명 광학 맥락. 추출 문제나 Yablonovitch 한계와 연결하지 않음 | 필수 (확인한 범위에서 가장 이른 일반형 진술) |
| 4 | Zhang, Hsu, Miller, *Optica* 6, 1321 (2019), doi:10.1364/OPTICA.6.001321 **[원문\*]**; O. D. Miller, F. Monticone, arXiv:2605.24738v2 (2026), §4.4 **[원문]**; D. A. B. Miller, arXiv:2608.12578 (2026) **[원문\*]** | 파동 brightness theorem. 출력 채널 전력의 상위 N개 합은 입력 coherency 행렬 고유값의 상위 N개 합 이하(약한 majorization). Miller & Monticone은 수동 선형계가 채널 전력 스펙트럼을 "doubly-(sub)stochastic mixing"으로만 재분배한다고 정리 | band(채널 집합)별 상한과 band 간 결합 제약이 바로 따름. "수동 = 혼합" 해석은 주장 2와도 연결됨 | 모드·파동 형식. OLED 대상 아님. 구조 평가 없음 | 권장 |
| 5 | Heikkilä, Oksanen, Tulkki, *Appl. Phys. Lett.* 99, 161110 (2011), doi:10.1063/1.3655155 **[원문]** | 텍스처 LED의 복사 전달 모형. 태양전지의 ergodic 한계가 LED 방출에도 적용됨(Kirchhoff). 기존 2n² 한계는 등방 방출일 때만 회복되고, 방출이 비등방이면 기하 인자로 들어감. 텍스처와 방출 패턴 설계로 2–4배 개선 여지 | 내부 각도 분포의 비등방성이 추출 상한을 바꿈 | 총량 상한. band별 정렬 상한이 아님 | 권장 |
| 6 | Muschaweck, Ries, *Proc. SPIE* 5942, 594206 (2005), doi:10.1117/12.615874 **[초록]**; Muschaweck, US 9,405,861 B2 (우선일 2010) **[원문\*]**; Lenef, Kelso, Piquette, *Opt. Lett.* 39, 3058 (2014), doi:10.1364/OL.39.003058 **[초록]** | 2005 논문과 특허는 위상공간을 radiance 순으로 정렬해 광원의 flux–étendue 특성 곡선을 만듦. 2014 논문은 Muschaweck 정리를 회절 표면으로 확장 | radiance 정렬 | 광원 특성화, 형광체 | 참고 (2014는 권장) |
| 7 | Ries, Shatz, Bortz, Spirkl, *JOSA A* 14, 2855 (1997), doi:10.1364/JOSAA.14.002855 **[초록]**, SPIE 동반 논문 *Proc. SPIE* 3139, 47 (1997) **[원문\*]** | 회전대칭 장치에서는 skewness 보존 때문에 étendue만으로 정한 한계보다 엄격한 성능 한계가 생김 | 불변량별 회계가 band별 회계와 비슷함 | 균일 radiance 가정. 회전대칭 구조에만 해당 | 참고 (애매함: 회전대칭 렌즈의 band 상한이 더 조여질 수 있음) |
| 8 | Schnitzer, Yablonovitch, Caneau, Gmitter, Scherer, *Appl. Phys. Lett.* 63, 2174 (1993), doi:10.1063/1.110575 **[원문\*]** | 완전 각도 무작위화 photon-gas 모형. dome 공식 n_l²/4n²을 위상공간 논증에서 나오는 일반 상한으로 언급 | 등방(균일 radiance) 경우의 상한 | 비등방 radiance 없음 | 참고 |
| 9 | 각도 제한 광포획 묶음: Campbell, Green, *IEEE TED* 33, 234 (1986) **[초록]**; Ulbrich et al., *pss(a)* 205, 2831 (2008) **[초록]**; Yu, Fan, *APL* 98, 011106 (2011) **[원문\*]**; Kosten et al., *Light Sci. Appl.* 2, e45 (2013) **[원문\*]**; Green, Zhou, *Nat. Commun.* 16, 251 (2025) **[원문\*]** | 수광·방출 각도를 제한하면 경로 증강 한계가 4n²/sin²θ로 커짐. Yu & Fan은 각도 합 규칙으로 band 사이에 증강을 배분함 | 각도 제한, band 배분 | 흡수 방향 문제. 균일 내부 radiance 가정 | 참고 |
| 10 | Gardner, *Phys. Fluids* 6, 839 (1963), doi:10.1063/1.1706823 **[2차]** | Liouville 제약 아래 위상공간 분포를 감소 순으로 재배열("restacking")해 가용 에너지 상한을 구함 | 감소 재배열로 상한을 구하는 다른 분야의 사례 | 플라즈마 | 참고 |

참고로, 현 OE 원고는 radiance 상한 근거로 Winston, Jiang, Ricketts (*Adv. Opt. Photon.* 10, 484, 2018)만 인용합니다(원고 ref. 12). 위 1–3번을 추가해야 합니다.

### 2-2. 각도 혼합 지표 (주장 2)

| 순위 | 서지 (저자, 연도, 저널, DOI) | 무엇을 했나 | 우리 정의와 같은 점 | 다른 점 | 인용 필요성 |
|---|---|---|---|---|---|
| 1 | Kim, Lee, Yoo, Gather, arXiv:2610.12065 (2026) **[원문]** (자기 선행) | 기판 안 빛의 각도 분포를 B_T(θ), B_R(θ_i, θ_r), R(θ)로 통과마다 행렬 곱해 전파(Suppl. Note 1, Eqs. S5–S8). η = p/[p+(1−p)A′]는 매 통과 각도 무작위화를 가정한 해. "microlens arrays return light largely at its original angle instead of randomizing it"이라고 지적하고, 통과별 p_k와 가중 평균 p̄를 씀 | 같은 행렬 틀. 각도 보존 문제를 이미 제기함 | 고유값, λ₂, N_mix 없음. 구조 분류 지표 없음 | 필수 |
| 2 | David, *J. Display Technol.* 9, 301 (2013), doi:10.1109/JDT.2013.2240373 **[원문]** | 거친 GaN 표면의 각도 분해 산란 행렬 S로 bounce를 합산함: C_N = T·Σ(R_p·S)^(n−1)·I (Eq. 6). 탈출까지의 bounce 수가 큰 각에서 발산하므로 "일정한 bounce 수" 개념은 맞지 않다고 봄(Fig. 18). 등방 분포에서 출발해도 bounce가 거듭될수록 ergodicity가 줄고 큰 각에 강도가 몰림(Fig. 19). 평균 bounce 수 N ≈ 9.5(R_p^N = C_ex로 정의)는 부적합하다고 지적 | 통과별 각도 재분배 행렬. 각도 기억 때문에 photon-gas 가정이 깨짐을 보임 | λ₂ 등 스펙트럼 지표 없음. LED 칩 내부의 거친 표면. 구조 분류 없음 | 필수 |
| 3 | Manzoor, Filipič, Onno, Topič, Holman, *J. Appl. Phys.* 127, 063104 (2020), doi:10.1063/1.5131173 **[원문]** | 텍스처 Si 안의 각도 분포를 통과마다 남은 빛으로 정규화해 추적함. Lambertian 대비 RMSD(Eq. 5)로 수렴을 평가. 실제 랜덤 피라미드는 2 pass, 이상적 랜덤 피라미드는 8 pass에 Lambertian 도달 | "각도를 잊는 데 필요한 통과 수"를 정량화함 | 경험적 지표. 초기 조건에 따라 달라짐. 고유값 없음. 흡수(태양전지) 대상 | 필수 |
| 4 | Battaglia, Boccard, Haug, Ballif, *J. Appl. Phys.* 112, 094504 (2012), doi:10.1063/1.4761988 **[원문]** (a, b 정의는 **[원문\*]**) | Lambertianity 인자 a와 탈출 분율 b(Lambertian이면 1/n²)를 정의함. 매 통과 같은 분포인 경우와 첫 통과만 비Lambertian인 경우를 식으로 다룸. 연속 산란이 Lambertian에 접근하는 현실적 기술은 "yet to be developed"라고 명시 | 부분 무작위화와 통과별 수렴을 미해결 과제로 제시함 | 지표를 정의하지 않음 | 필수 (동기 부여 인용) |
| 5 | Eisenlohr et al. (OPTOS), *Opt. Express* 23, A502 (2015), doi:10.1364/OE.23.00A502 **[원문]**. 같은 계열: Tucher et al., *Opt. Express* 23, A1720 (2015)과 24, A1083 (2016) **[원문\*]**; Li et al. (AMF), *Opt. Express* 23, A1707 (2015) **[원문\*]**; Santbergen et al. (GenPro4), *IEEE J. Photovolt.* 7, 919 (2017) **[원문\*]**; Pearce (RayFlare), *J. Open Source Softw.* 6, 3460 (2021) **[원문\*]** | 각도 채널 전력 벡터를 재분배 행렬로 통과마다 전파함. OPTOS는 열합 노름 1 이하로 Neumann 급수의 수렴을 보장함. Lambertian 산란 행렬은 모든 입사각에서 같은 분포를 만들고, 평면은 대각 행렬 | 같은 행렬 대상. Lambertian 행렬은 λ₂ = 0, 평면은 λ₂ = 1에 해당 | 고유값 해석 없음. 노름과 스펙트럼 반경은 수렴 판정에만 씀. 구조 지표 없음 | 필수 (최소 Eisenlohr 2015) |
| 6 | Altmann, Portela, Tél, *Rev. Mod. Phys.* 85, 869 (2013), doi:10.1103/RevModPhys.85.869 **[원문]**. 함께 볼 문헌: Altmann, *Phys. Rev. A* 79, 013830 (2009) **[원문\*]** | 누설 카오스계의 탈출률은 conditionally invariant measure로 정해짐. 단순 추정(Sabine·Eyring형)은 "비정반사 충돌이 궤적을 즉시 무작위화"하는 등 비현실적 조건에서만 정확함. 유전체 공진기의 방출 패턴도 이 불변 집합으로 설명됨(2009) | "λ₂ = 0일 때만 단순식이 정확하다"와 같은 구조. 남은 빛의 분포가 고정 분포로 수렴함 | 선도 고유값(c-measure)만 다룸. λ₂ 지표 없음. 광추출 구조가 아님 | 권장 |
| 7 | Kovačič, *Opto-Electron. Rev.* 30, e141542 (2022), doi:10.24425/opelre.2022.141542 **[원문]** | OLED 외부 텍스처에서 반사광을 스택 반사율로 감쇠시켜 8번째 통과까지 추적함. 처음 6–8회 통과에서 99% 이상이 추출되거나 흡수됨. 피라미드 텍스처는 첫 반사 뒤 내부 분포가 쌍극자 방향과 무관하게 거의 같아짐. 텍스처마다 T(θ)는 크게 달라도 최종 효율 차이는 작음 | OLED에서 통과별 각도 수렴(빠른 혼합)을 보임 | 행렬과 지표 없음. MLA 없음 | 권장 |
| 8 | Deckman, Roxlo, Yablonovitch, *Opt. Lett.* 8, 491 (1983), doi:10.1364/OL.8.000491 **[초록]**; Deckman, Witzke, Wronski, Yablonovitch, US 4,554,727 (1985) **[원문\*]** | 불완전 무작위화 비율 β를 정의: 측정 증강 = β × 완전 무작위화 증강 + (1−β) × 평면. 측정값 β는 샌드블라스트 0.4–0.6, HF 식각 0.1–0.4 **[원문\*]** | 무작위화 정도를 스칼라로 정의함 | 통과 동역학과 고유값 없음 | 권장 |
| 9 | Heikkilä et al. 2011 (2-1의 5번) **[원문]** | 출사면 반사를 정반사와 확산의 혼합으로 보고 albedo χ로 매개변수화함 | 1-파라미터 부분 무작위화(χ = 1이 완전 무작위화) | λ₂나 통과 수 해석 없음 | 참고 |
| 10 | 유사 개념. 복사 전달: van Rossum, Nieuwenhuizen, *Rev. Mod. Phys.* 71, 313 (1999), doi:10.1103/RevModPhys.71.313 **[원문\*]**; Rotter, Gigan, *Rev. Mod. Phys.* 89, 015005 (2017) **[원문\*]**. 마르코프 연쇄: Levin, Peres, *Markov Chains and Mixing Times*, 2nd ed., AMS (2017), doi:10.1090/mbk/107, §12.2 **[원문]** | 복사 전달: 수송 평균 자유 경로 ℓ* = ℓ/(1−g), 즉 방향을 잊는 데 약 1/(1−g)번 산란. 마르코프 연쇄: absolute spectral gap γ* = 1 − max\|λ\| (λ ≠ 1), relaxation time t_rel = 1/γ* | 1/(1−λ₂)와 형식이 같음. g는 회전 불변 산란 연산자의 l = 1 성분 고유값 | 산란 사건당 양이거나 일반 연쇄의 양. 구조의 통과당 재분배가 아님 | 필수 (심사에서 지적될 가능성이 커 미리 인용) |
| 11 | 마르코프 연쇄·휴리스틱. Hébert, Hersch, Becker, *JOSA A* 24, 2628 (2007), doi:10.1364/JOSAA.24.002628 **[원문]**; Kutayiah et al., *Opt. Express* 29, 4249 (2021), doi:10.1364/OE.416670 **[원문\*]**; McTague, Ratnaweera, Sheldon, arXiv:2607.20783 (2026) **[원문\*]**; RayFlare GitHub issue #65 (2024) **[원문\*]** | Hébert 등은 다층의 다중 반사·투과를 마르코프 연쇄로 기술하며, 각 면의 각도 분포가 반복 사이에 일정하다고 가정함. LSC 마르코프 모형은 정상 고유벡터(고유값 1)와 스펙트럼 반경(수렴 조건)만 씀. RayFlare issue는 "N회 상호작용 뒤 Lambertian" 휴리스틱을 제안함 | 마르코프 연쇄 틀. 기억 없는 가정을 명시함 | λ₂ 지표 없음 | 참고 |
| 12 | Greiner, *Jpn. J. Appl. Phys.* 46, 4125 (2007), doi:10.1143/JJAP.46.4125 **[초록]**; Levell et al., *Proc. SPIE* 8829, 88291L (2013), doi:10.1117/12.2023021 **[원문\*]** | Greiner: 스택 반사율이 기판→공기 전달의 궁극 한계(평면 약 50%, 구조화 80% 이상) **[초록]**, 추출 구조 특성화 도표 제시 **[초록]**. Levell: 평균 반사율이 약 60%보다 낮으면 산란 필름 이득이 없음을 실험으로 보임 **[원문\*]** | 소자 반사(A′에 해당)가 지배 변수 | 각도 기억과 고유값 없음. Greiner는 원문 미확인 | 권장 (Nature 원고가 이미 인용. Greiner는 애매함) |
| 13 | Shiang, Duggal, *J. Appl. Phys.* 95, 2880 (2004), doi:10.1063/1.1644037 **[초록]**; US 7,245,074 B2 **[원문\*]** | OLED 산란층의 이산 방향 복사 전달. 각도 채널 결합 행렬을 고유분해해 해를 구함. "많은 추출 방식을 매개변수화할 수 있는 일반적 방법"이라고 기술 **[초록]** | 각도 채널 행렬과 고유분해 | 고유값이 깊이 방향 전파 연산자의 것. 통과별 혼합 지표가 아님 | 참고 |
| 14 | LCD 백라이트 특허. O'Neill, Cobb, US 5,917,664 A (출원 1996) **[원문\*]**; Freier, US 7,898,520 B2 (우선일 2005) **[원문\*]** | 프리즘 필름 재순환에는 확산 반사체의 방향 무작위화가 필요함. 확산판을 넣어 프리즘 필름이 같은 빛을 반복해 되돌리지 않게 함. 부품을 각도 bin BSDF 행렬로 쌓아 계산 | 각도 기억이 재순환 이득을 막는다는 정성적 인식. 행렬 계산 | 스펙트럼 지표 없음 | 참고 |

### 2-3. 확인 권장 항목 요약 (지침 B5–B7)

- **공통 상한, 통일 분류 (B5)**
  - 모든 구조에 같은 열역학 상한이 적용된다는 진술은 Lenef(2018: 기하·회절·산란 표면 모두)와 Bortz(1997)에 있습니다.
  - OLED 구조를 공통 틀로 매개변수화하려는 시도로는 Shiang & Duggal(2004, 복사 전달) [초록]과 Greiner(2007, 특성화 도표) [초록]가 가장 가깝습니다.
  - (p, N_mix)와 A′로 구조와 소자를 나누는 분류는 찾지 못했습니다.
- **출력이 구조와 무관하게 Lambertian 근처로 수렴 (B6)**
  - 관찰 사례:
    - Kovačič(2022) [원문]
    - Manzoor(2020) [원문]
    - David(2013) 각주 7: 입사각 전체로 적분한 확산 산란은 대개 Lambertian이지만, 입사각별로는 아님 [원문\*]
    - Song et al.(*Nat. Commun.* 9, 3207, 2018): scatterance가 높으면 Lambertian형 출력 [원문\*]
  - 이론적 근거로 쓸 만한 문헌(우리 해석): 균일(Lambertian) 분포가 수동 재분배의 고정점이라는 점을 함의합니다.
    - Miller & Monticone(2026): 수동계는 doubly-(sub)stochastic 혼합만 한다 [원문]
    - GenPro4: 열합 1, 대칭(가역) 행렬 [원문\*]
- **ML 최적화를 물리 상한과 비교한 연구 (B7):** 찾지 못했습니다.
  - Zhao et al.(*Adv. Funct. Mater.* 2026, e78939, doi:10.1002/adfm.78939, DQN) [원문]은 BO와 무작위 탐색과만 비교합니다.
  - Kim et al.(*Nanophotonics* 10, 4533, 2021, DNN) [초록]과 Huang et al.(*Photonics Res.* 10, 269, 2022, doi:10.1364/PRJ.441188, DDQN) [초록]도 상한을 제시하지 않습니다.

---

## 3. 용어 충돌 확인 결과

| 후보 이름 | 다른 분야의 기존 사용 | 충돌 위험 | 권고 |
|---|---|---|---|
| memory effect, angular memory (effect) | 산란 매질의 optical memory effect: Feng, Kane, Lee, Stone, *PRL* 61, 834 (1988), doi:10.1103/PhysRevLett.61.834; Freund, Rosenbluh, Feng, *PRL* 61, 2328 (1988), 제목 "Memory effects in propagation of optical waves through disordered media", doi:10.1103/PhysRevLett.61.2328. 산란 매질 통과 이미징: Bertolotti et al., *Nature* 491, 232 (2012); Katz et al., *Nat. Photonics* 8, 784 (2014); Osnabrugge et al., *Optica* 4, 886 (2017, "generalized optical memory effect") [초록] | 높음 | 쓰지 않음. 현 OE 원고 Sec. 4.1과 Note S9의 "angular memory"도 "angle retention" 등으로 교체 권장 |
| angular mixing | 조명 광학의 mixing rod·적분기에서 각도·색 균일화: Bonenberger et al., *Adv. Opt. Technol.* 5, 157 (2016), doi:10.1515/aot-2015-0056 [초록]. Miller & Monticone(2026)은 수동 재분배를 "doubly-(sub)stochastic mixing"이라 부름 [원문] | 중간 (뜻이 호환됨) | 사용 가능. 처음 나올 때 정의를 명시 |
| angular mixing number | 정확한 구문의 사용례를 찾지 못함(OpenAlex 제목·초록, 웹 검색). 입자물리의 "mixing angle"은 다른 용어 | 낮음 | 권장 후보 |
| mixing number, N_mix | 건조 공학 침투 모형의 교반 횟수 N_mix: Schlünder, Mollekopf, *Chem. Eng. Process.* 18, 93 (1984) [2차]. 축열조 성층의 "MIX number": Davidson, Adams, Miller, *J. Sol. Energy Eng.* 116, 94 (1994) [초록]. 마르코프 이론의 mixing time t_mix는 relaxation time t_rel과 다른 양 | 중간 | 기호를 유지하면 "t_mix가 아니라 t_rel에 해당한다"고 명시. N_rel 같은 기호도 고려 |
| angular persistence | 세포 운동·화학주성의 방향 지속성: Kirkegaard, Goldstein, *IMA J. Appl. Math.* 83, 700 (2018) [초록]. 고분자 persistence (persistence length·number) | 중간 | 피함 |
| angular relaxation (time) | 분자 회전 이완: Fixman, Rider, "Angular relaxation of the symmetrical top," *J. Chem. Phys.* 51, 2425 (1969) [초록]. 전이 핵 아래 ⟨P_n(cosθ)⟩의 감쇠를 다루므로 수학적으로 같은 대상 | 중간–높음 (시간 영역의 뜻이 강함) | 피함 |
| angular randomization number | 사용례 없음. 인접 용어로 Deckman의 "randomization fraction β" | 낮음 | 대안 후보 |
| angular decorrelation | 스펙클의 조명각 상관: Shirley, George, *JOSA A* 6, 765 (1989) [초록]. 입자물리 jet의 azimuthal decorrelation | 중간 | 피함 |
| number of passes to randomize (서술형) | Manzoor(2020) [원문]; 카드 섞기 횟수: Bayer, Diaconis, *Ann. Appl. Probab.* 2, 294 (1992) [초록] | 낮음 | 서술 표현으로는 안전 |

---

## 4. 차별화 문장 초안 (영문)

**안 A (3문장)**

> Radiance-ordered bounds, which fill the accessible output phase space starting from the brightest rays, are established in nonimaging optics [Bortz 1997] and have been applied to single-pass extraction from LEDs and phosphors with non-Lambertian internal radiance [Muschaweck 2013; Lenef 2014, 2018]. Here we resolve this bound into polar-angle bands for the anisotropic radiance that the OLED microcavity launches into the substrate, and use it to place microlens, scattering, and planar external films on a common scale. Multi-pass extraction has so far been modelled either with photon-gas expressions that assume complete angular randomization at every pass [Yablonovitch 1982; Schnitzer 1993; Kim 2026] or with angle-resolved redistribution matrices evaluated numerically [Eisenlohr 2015; David 2013; Kim 2026], and the degree of randomization has been quantified only empirically [Deckman 1983; Battaglia 2012; Manzoor 2020]; we instead characterize each structure by the spectral gap of its pass-to-pass redistribution matrix, N_mix = 1/(1 − λ₂), the per-pass counterpart of the 1/(1 − g) factor of radiative transfer [van Rossum 1999] and of the relaxation time of a Markov chain [Levin 2017], which recovers the memoryless result in the limit λ₂ → 0.

**안 B (2문장, 짧은 버전)**

> Building on radiance-ordered extraction bounds for LEDs and phosphors [Muschaweck 2013; Lenef 2018], we resolve the single-pass bound into polar bands for the anisotropic substrate radiance of an OLED and benchmark microlens, scattering, and planar films against it. Beyond the escape probability, we describe each structure by the spectral gap of its pass-to-pass angular redistribution matrix, N_mix = 1/(1 − λ₂), the number of passes needed to lose angular information, which earlier matrix models computed implicitly but did not extract [Eisenlohr 2015; David 2013; Kim 2026] and whose memoryless limit λ₂ = 0 corresponds to η = p/[p + (1 − p)A′].

---

## 5. 원문을 확인하지 못한 논문 목록

| 논문 | DOI 등 | 사유 | 판정에 미치는 영향 |
|---|---|---|---|
| Greiner, *JJAP* 46, 4125 (2007) | 10.1143/JJAP.46.4125 | IOP 유료, 기관 접근 없음 | "추출 구조 특성화 도표"가 구조 분류 틀일 수 있음. 애매함. 원문 확인 1순위 |
| Neyts, Nieto, *JOSA A* 23, 1201 (2006) | 10.1364/JOSAA.23.001201 | Optica 구독 필요(로그인·결제 화면) | 산란 확률과 거울 반사율의 해석 모형. 애매함. 확인 권장 |
| Shiang, Duggal, *JAP* 95, 2880 (2004); Shiang, Faircloth, Duggal, *JAP* 95, 2889 (2004) | 10.1063/1.1644037; 10.1063/1.1644038 | 출판사 페이지에 본문 HTML이 없고 PDF를 열지 못함. 같은 내용의 특허 US 7,245,074 B2로 보완 | 낮음 |
| Lenef, Kelso, Piquette, *Opt. Lett.* 39, 3058 (2014) | 10.1364/OL.39.003058 | Optica 봇 확인 | 낮음 (2018 논문이 포괄) |
| Muschaweck, Ries (2005); Muschaweck (2012); Bortz, Shatz (2006) | 10.1117/12.615874; 10.1117/12.930605; 10.1117/12.678600 | SPIE 비공개 | 낮음 |
| Winston, Miñano, Benítez, *Nonimaging Optics* (2005), Shatz & Bortz 장 "Consequences of symmetry" | ISBN 978-0-12-759751-5; 10.1016/B978-012759751-5/50010-5 | 유료 도서 | radiance 정렬 상한이 더 이른 시기에 진술됐을 수 있음 |
| Winston, Jiang, Ricketts, *Adv. Opt. Photon.* 10, 484 (2018) | 10.1364/AOP.10.000484 | Optica 봇 확인 | 낮음 (원고가 이미 인용) |
| Ries (1982); Smestad et al. (1990); Markvart (2007/2008); Ries et al. *JOSA A* (1997) | 10.1364/JOSA.72.000380; 10.1016/0165-1633(90)90047-5; 10.1088/1464-4258/10/01/015008; 10.1364/JOSAA.14.002855 | 유료. Ries 1997은 SPIE 동반 논문 원문으로 보완 | 낮음 |
| Campbell, Green (1986, 1987); Ulbrich et al. (2008); Peters et al. (2010); Höhn et al. (2014); Yu, Raman, Fan, *PRL* (2012); Callahan et al. (2012); Guo, Fan, *PRL* (2023) | 10.1109/T-ED.1986.22472; 10.1063/1.339189; 10.1002/pssa.200880457; 10.1016/j.solmat.2010.04.009; 10.1364/OE.22.00A715; 10.1103/PhysRevLett.109.173901; 10.1021/nl203351k; 10.1103/PhysRevLett.130.146202 | 유료 또는 Optica 봇 확인 | 낮음 (균일 radiance를 가정한 각도 제한 한계) |
| Green, *Prog. Photovolt.* 10, 235 (2002); Deckman et al., *Opt. Lett.* (1983); Schuster et al., *OE* (2014); Yablonovitch, Cody, *IEEE TED* 29, 300 (1982) | 10.1002/pip.404; 10.1364/OL.8.000491; 10.1364/OE.22.00A542; 10.1109/T-ED.1982.20700 | 유료 또는 차단. Deckman의 내용은 특허 원문으로 보완 | 낮음 |
| Windisch et al. (2001, 2002); Krummacher et al. (2006); Rahadian et al. (2007); Lin et al. (2008); Chen et al. (2010); Möller, Forrest (2002); Wei, Su (2004); Kumar et al. (2011) | 10.1063/1.1397758; 10.1109/2944.999177; 10.1063/1.2234550; 10.1117/1.2821869; 10.1364/OE.16.011044; 10.1364/OE.18.003238; 10.1063/1.1435422; 10.1364/OPEX.12.005777; 10.1016/j.optcom.2011.06.005 | 유료 또는 Optica 봇 확인 | 낮음 (초록상 통과별 각도 분석이나 스펙트럼 지표 없음) |
| Kovačič et al., *ACS Photonics* 5, 422 (2018); *Opt. Express* 29, 23701 (2021); *Adv. Photonics Res.* 6, 2400188 (2025) | 10.1021/acsphotonics.7b00874; 10.1364/OE.428021; 10.1002/adpr.202400188 | 유료 또는 Optica 봇 확인 | 낮음 (같은 저자의 2022 논문 원문 확인) |
| Şahin, Ilan, *JOSA A* 30, 813 (2013); Gloge, *Bell Syst. Tech. J.* 51, 1767 (1972) | 10.1364/JOSAA.30.000813; 10.1002/j.1538-7305.1972.tb02682.x | 유료 | 낮음 (유사 개념) |

---

## 부록. 심사 대응 메모 (조사 중 확인한 기술적 주의점)

1. **손실이 있는 행렬의 λ₂.** 투과와 손실이 있으면 재분배 행렬은 substochastic입니다.
   - 열을 정규화한 조건부 연쇄를 쓰거나, |λ₂/λ₁|을 쓴다고 정의해야 합니다.
   - 남은 빛이 수렴하는 분포는 quasi-stationary 분포입니다. 이는 Altmann(2013)의 c-measure, David(2013)의 큰 각 집중과 같은 대상입니다.
2. **음수·복소 고유값.** 정반사 성분이 각을 번갈아 바꾸면 절댓값이 큰 음수 고유값이 생길 수 있습니다.
   - Levin & Peres의 t_rel은 absolute spectral gap(1 − max|λ|) 기준입니다.
   - 따라서 N_mix = 1/(1 − |λ₂|)로 정의하거나, λ₂가 절댓값 최대임을 확인해 두는 것이 안전합니다.
3. **실수 스펙트럼과 Lambertian 고정점(우리 해석).** 가역(reciprocal) 구조의 반사 행렬은 étendue 가중(cosθ sinθ dθ)에 대해 detailed balance를 만족합니다.
   - 대각인 소자 반사율 R(θ)를 곱하거나 열을 정규화해도 대칭 행렬과 닮음이므로 고유값은 실수입니다.
   - 무손실이면 균일 radiance(Lambertian)가 정상 분포입니다. 이 점은 GenPro4의 "열합 1 + 대칭" [원문\*], Miller & Monticone의 doubly-stochastic 혼합 [원문]과 맞습니다.
   - "비선택 수동 구조의 출력 구성이 Lambertian 분할 근처에 고정된다"는 관찰의 근거로 쓸 수 있습니다.
4. **통과 수 정의의 선택.** 1/(1−λ)와 e-folding 통과 수 −1/ln|λ|는 다릅니다. 예를 들어 λ = 0.5에서 각각 2와 1.44입니다. 어느 쪽을 쓰는지 본문에 명시해야 합니다.
5. **"등방"이라는 말.** Muschaweck(2013)과 Lenef(2018)는 "Lambertian"을 균일 radiance로, "isotropic"을 체적 등방 발광(표면에서 L ∝ 1/cosθ)으로 씁니다. 주장 1의 "기판 안 빛이 등방이면"은 "uniform (Lambertian) radiance"로 쓰는 것이 안전합니다.
6. **band 상한의 동시 달성.** band별 상한은 각각 달성할 수 있지만 동시에 달성할 수는 없습니다.
   - 여러 band의 합집합에 대한 상한(majorization형 결합 제약)을 함께 제시하면 Lenef(2018)와의 차이가 분명해집니다.
   - Miller & Monticone(2026)의 파동 버전이 근거가 됩니다.
7. **인용 보강.** 현 OE 원고에는 Bortz 1997, Muschaweck 2013, Lenef 2018, David 2013, Manzoor 2020, Battaglia 2012, Eisenlohr 2015, van Rossum 1999(또는 Rotter 2017), Levin & Peres 2017이 없습니다. 2절의 "필수" 항목부터 추가해야 합니다.

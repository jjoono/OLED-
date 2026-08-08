# References — OLED microlens-array limits paper

동기화 기준: **paper_plan_v7.md** (2026-08-08 반영). 주 목록은 v7 계획서의 참고문헌 [1]–[15]와
번호·순서 동일. 검증 상태: URL/DOI는 웹검색으로 실재 확인된 것 유지. 서지 상세(권/페이지)는
확인된 것만 표기, 미확인은 "[서지정보 확인 필요]"로 표시 (임의로 채우지 말 것).
그룹 자기인용(⚑)은 self-overlap 관리 대상.

---

## 주 참고문헌 목록 (v7 plan [1]–[15], 번호 고정)

### [1] Brütting 2013 — OLED outcoupling review
- Brütting, W.; Frischeisen, J.; Schmidt, T. D.; Scholz, B. J.; Mayr, C.
  "Device Efficiency of Organic Light-Emitting Diodes: Progress by Improved Light
  Outcoupling." *Phys. Status Solidi A* **2013**, *210*, 44–65.
- 내용: 유기층/기판 도파 모드, 금속 전극 손실(SPP), TIR 등 OLED 광추출 손실 채널과
  개선 전략의 표준 리뷰.
- 사용처: §1 서론 첫 문단 — 평면 OLED의 손실 모드 배경 근거 ([1,2]로 병기).

### [2] Yablonovitch 1982 — statistical ray optics
- Yablonovitch, E. "Statistical Ray Optics." *J. Opt. Soc. Am.* **1982**, *72*, 899–907.
  https://opg.optica.org/josa/abstract.cfm?uri=josa-72-7-899
- 내용: 방출↔흡수 detailed balance, 4n²/2n² 한계, 손실이 유한 천장을 만든다는 고전 결과.
- 사용처: §1 서론 손실 모드 배경 ([1,2]). (구판에서 쓰던 "보편성 정리" 근거 역할은
  P0 리뷰 후 폐기 — 일반 배경 인용으로만 사용.)

### [3] Möller & Forrest 2002 — classic MLA
- Möller, S.; Forrest, S. R. "Improved Light Out-Coupling in Organic Light Emitting
  Diodes Employing Ordered Microlens Arrays." *J. Appl. Phys.* **2002**, *91*, 3324–3327.
- 내용: 외부 MLA로 OLED 기판 모드를 추출하는 고전(원조) 논문. 전기 구조 무손상
  외부 광학 접근의 출발점.
- 사용처: §1 서론 — 외부 MLA가 대면적/유연 OLED의 매력적 해결책이라는 배경 [3–5].

### [4] Wrzesniewski 2012 — molded polymer MLA (top-emitting)
- Wrzesniewski, E.; Eom, S.-H.; Cao, W.; Hammond, W. T.; Lee, S.; Douglas, E. P.; Xue, J.
  "Enhancing Light Extraction in Top-Emitting Organic Light-Emitting Devices Using
  Molded Transparent Polymer Microlens Arrays." *Small* **2012**, *8*, 2647–2651.
  https://doi.org/10.1002/smll.201102662
- 내용: 성형 폴리머 MLA로 top-emitting OLED 광추출 향상 — (준)반구형 MLA의 높은 효율 실증.
- 사용처: §1 서론 — 반구형/준반구형 MLA의 고효율 실증 사례 [4,5].

### [5] Qu 2018 — embedded hemispherical MLA, very high EQE
- Qu, Y.; et al. "Efficient, Nonintrusive Outcoupling in Organic Light Emitting Devices
  Using Embedded Microlens Arrays." *ACS Photonics* **2018**.
  https://doi.org/10.1021/acsphotonics.8b00255 [서지정보 확인 필요: 권/페이지]
- 내용: embedded hemispherical MLA로 매우 높은 EQE 달성 (nonintrusive 외부 광학).
- 사용처: §1 서론 — embedded hemispherical MLA의 very high EQE 보고 [4,5].
- ⚠ **번호 불일치 확인 필요**: plan §4.1은 "[5] 본 연구 그룹의 선행 MLA-OLED 연구
  (실험 배광·EQE 검증)"라고 쓰는데, 리스트 [5]는 Forrest 그룹의 Qu 2018임.
  §4.1의 의도는 아래 미채택 섹션의 ⚑ Kim 2023 *ACS Photonics* (자기 그룹)일 가능성이
  높음 — 원고 확정 전에 [5]를 분리·재번호하거나 §4.1 인용을 수정할 것.

### [6] Kim 2021 — DNN inverse design of OLED
- Kim, S.(?); et al. "Inverse Design of Organic Light-Emitting Diode Structure Based on
  Deep Neural Networks." *Nanophotonics* **2021**, *10*.
  https://onlinelibrary.wiley.com/doi/10.1515/nanoph-2021-0434
  [서지정보 확인 필요: 저자 명단·호/페이지]
- 내용: DNN 기반 OLED 구조 역설계 — inverse design이 OLED 광학에 적용되는 흐름의 예.
- 사용처: §1 서론 — freeform illumination optics/inverse design의 부상 [6,7].

### [7] Appl. Opt. 2025 — freeform MLA for microLED packaging
- "Design of Freeform Microlens Arrays with Prescribed Luminance Distributions for
  MicroLED Optical Packaging." *Appl. Opt.* **2025**, *64*, 7875.
  https://opg.optica.org/ao/abstract.cfm?uri=ao-64-27-7875 [서지정보 확인 필요: 저자]
- 내용: 목표 배광(luminance distribution)을 내는 freeform MLA 설계 — micro-LED 패키징.
- 사용처: §1 서론 — freeform/역설계가 micro-LED·OLED 광학 패키징에 적용됨 [6,7].

### [8] Buhl 2023 — resonance-based directional OLED (DBR angular filter)
- Buhl, M.; et al. "Resonance-Based Directional Light Emission from Organic
  Light-Emitting Diodes." *Adv. Photonics Res.* **2023**.
  https://doi.org/10.1002/adpr.202200143 [서지정보 확인 필요: 권/논문번호]
- 내용: 공진(DBR) 기반 방향성 OLED. DBR을 입사각 필터로 쓰면 저각 방출을 억제하지만
  "light interacts more often with absorbing materials inside the device"로 손실 증가를
  명시 — 우리의 recycling loss 트레이드오프 계산과 정합.
- 사용처: §2.4 — 회절/공진 요소는 방향성 광분배 가능(스코핑 가드레일) [8,9];
  §2.6 — angular-selective recycling이 이미 알려진 개념임의 근거 [8,14,15].
  ⚠ 선행 확정 — 새 구조 "제안"으로 쓰지 말고 known-route reference로만 배치.

### [9] Commun. Eng. 2025 — metasurface directional microLED
- "Enhanced and Directional Electroluminescence from MicroLEDs Using Metallic or
  Dielectric Metasurfaces." *Commun. Eng.* **2025**.
  https://doi.org/10.1038/s44172-025-00401-w [서지정보 확인 필요: 저자·논문번호]
- 내용: 메타표면(금속/유전체)이 위상경사로 µLED 방출에 방향성 부여 — 실제로 조향됨.
- 사용처: §2.4 — 포화 결론을 "coextensive refractive MLA" 클래스로 한정하는
  스코핑 가드레일 [8,9] (metasurface/aperture expansion은 조향 가능).

### [10] Winston 2018 — nonimaging optics tutorial
- Winston, R.; Jiang, L.; Ricketts, M. "Nonimaging Optics: A Tutorial."
  *Adv. Opt. Photon.* **2018**, *10*, 484–511.
  https://opg.optica.org/aop/fulltext.cfm?uri=aop-10-2-484&id=389885
- 내용: étendue 보존, brightness theorem, sine-law 집광 한계 — 방향성 상한의 정석 언어.
- 사용처: §2.5 — radiance/étendue envelope를 "numerical frontier 해석 기준"으로
  사용하는 근거 [10,11] (global impossibility theorem 아님 — v7 표현 준수).

### [11] Rau 2007 — reciprocity relation (solar cells)
- Rau, U. "Reciprocity Relation between Photovoltaic Quantum Efficiency and
  Electroluminescent Emission of Solar Cells." *Phys. Rev. B* **2007**, *76*, 085303.
  https://doi.org/10.1103/PhysRevB.76.085303
- 내용: 준평형에서 EL 방출 ↔ 광전 흡수(EQE)의 상반성 관계 (태양전지 맥락).
- 사용처: §2.5 — [10]과 병기되어 radiance/étendue envelope 문단에 등장 [10,11].
- ⚠ **적용 범위 주의 — 상반성/광자재활용 일반 맥락에만 인용, radiance envelope
  논증에는 인용 금지 (외부 리뷰 지적).** 구판에서 "방향성 봉투 = 방출=흡수 동형"의
  일반화 Kirchhoff 근거로 오용했음. Rau 2007은 준평형 diode의 흡수↔방출 상반성이지
  수동 외부 광학계의 angular power 상한을 주지 않음. §2.5 본문 인용 위치를
  재검토하고, envelope 논증 자체는 [10] (étendue/brightness theorem)에만 기대게 할 것.

### [12] Xiang 2013 — microcavity systematic study
- Xiang, C.; Koo, W.; So, F.; Sasabe, H.; Kido, J. "A Systematic Study on Efficiency
  Enhancements in Phosphorescent Green, Red and Blue Microcavity Organic Light-Emitting
  Devices." *Light: Sci. Appl.* **2013**, *2*, e74. https://doi.org/10.1038/lsa.2013.30
- 내용: RGB 인광 microcavity OLED의 효율 향상 계통 연구 — cavity가 substrate-side
  source distribution을 바꾸는 대표 사례.
- 사용처: §2.6 — source/cavity engineering이 OLED angular emission control의
  확립된 방법이라는 근거 [12,13].

### [13] Song 2018 — lens-free scattering, >50% EQE ★ 직접 경쟁 선행
- Song, J.; et al. "Lensfree OLEDs with over 50% External Quantum Efficiency via
  External Scattering and Horizontally Oriented Emitters." *Nat. Commun.* **2018**,
  *9*, 3207. https://doi.org/10.1038/s41467-018-05671-x
- 내용: 렌즈 없이 외부 산란층 + 수평 배향 emitter만으로 EQE >50% 달성.
- 사용처: §2.6 — source/cavity engineering 경로의 근거 [12,13].
- ★ **DIRECT COMPETING PRIOR — §2.x에서 명시적으로 논의 필수.** 렌즈 없이도 총 EQE가
  매우 높게 나온다는 사실은 본 논문의 포화 논지를 오히려 강화함: 총 추출은 lens-free
  경로로도 달성 가능하며, 본 논문의 질문은 "freeform MLA가 그 위에 angular control을
  추가로 제공하는가"임. 이 구분을 본문에서 직접 서술하지 않으면 리뷰어가
  "MLA 자체가 불필요" 반론으로 사용할 수 있음.

### [14] 2022 — angle-selective optical film for OLED displays
- "Using Angle-Selective Optical Film to Enhance the Light Extraction of a Thin-Film
  Encapsulated 3D Reflective Pixel for OLED Displays." **2022**.
  https://pubmed.ncbi.nlm.nih.gov/36558597/
  [서지정보 확인 필요: 저자·저널명(MDPI 계열 추정)·권/논문번호]
- 내용: OLED 디스플레이 픽셀에 각도선택 광학 필름 직접 적용.
- 사용처: §2.6 — angular-selective recycling이 known route임의 근거 [8,14,15].

### [15] Kim 2017 — optical recycling with air-gapped bridges
- Kim, H.-J.; et al. "High Efficient OLED Displays Prepared with the Air-Gapped Bridges
  on Quantum Dot Patterns for Optical Recycling." *Sci. Rep.* **2017**, *7*, 43063.
  https://doi.org/10.1038/srep43063
- 내용: OLED 디스플레이 photon recycling 실측 이득 (적색 58.2%).
- 사용처: §2.6 — angular-selective recycling/photon recycling known-route 근거 [8,14,15].

---

## 미채택/보류 문헌 (not in v7 list)

v7 계획서의 [1]–[15]에 포함되지 않은 구판 항목. 사유 1줄씩. "보류"는 원고 확장 시
재고 가능, "미채택"은 현 논지에서 제외.

### 자기 그룹 (⚑ self-overlap 관리 — v7 리스트 누락 여부 확인 필요)
- ⚑ Kim, Kim, Park, Song, Kim, Moon, Yoo. "Toward Near-Foldable Surface Light Sources...
  Ultrathin Substrates Embedded with Micron-Scale Inverted Lens Arrays." *ACS Photonics*
  **2023**, *10*, 1775. https://pubs.acs.org/doi/10.1021/acsphotonics.3c00017
  [PDF: Kim2023_ACSPhotonics_IMLA_foldable.pdf]
  — 사유: v7 번호 리스트에 없음. 단, plan §4.1의 "본 연구 그룹의 선행 MLA-OLED 연구 [5]"가
  이 논문을 가리키는 것으로 보임 (위 [5] 항목의 번호 불일치 메모 참조). **리스트 편입 유력 —
  원고 확정 시 결정.** trans-scale(CPS+BSDF+LightTools) 방법·EQE 58% 선점 → self-overlap 관리.
- ⚑ "Near-planar light outcoupling structures with finite lateral dimensions for
  ultra-efficient and optical crosstalk-free OLED displays." *Nat. Commun.* **2025**.
  https://www.nature.com/articles/s41467-025-66538-6
  — 사유: v7 리스트에 없음 (보류). 스택+구조 공동설계 48% EQE — §2.6 source/cavity
  경로 논의 확장 시 재고; self-overlap 관리 대상.

### Inverse design / freeform 관련 (미채택)
- "Dual-Task Optimization Method for Inverse Design of RGB Micro-LED Light Collimator"
  (PMC11820347). — 사유: collimator(능동집속) 클래스는 본 논문 범위 밖; 역설계 부상
  근거는 [6,7]로 충분 (미채택).
- "Inverse design for scalable photonic systems." *Nat. Rev. Mater.* (2026).
  — 사유: 배경 리뷰 중복, [6,7]로 대체 (미채택).

### étendue/트레이드오프 "조각" 선행 (P0 리뷰 피벗으로 불필요)
- "Design Methodology for High Brightness Projectors." *J. Display Technol.* 4, 1 (2008).
  — 사유: theorem-style achievable-region 주장을 P0 리뷰 후 철회하면서 étendue
  트레이드오프 콜라주 논증이 사라짐; envelope 언어는 [10]으로 일원화 (미채택).
- US patents 10288884 / 11163165 / 10761327 (directed display architecture).
  — 사유: 동일 — theorem-style 논증 철회로 불필요, 특허 인용 회피 (미채택).
- "Analysis of light out-coupling from microlens array." *Opt. Commun.* (2011).
  — 사유: escape-cone 안팎 재활용 분석 — 현 v7 본문 서술에 인용 슬롯 없음 (보류;
  §2.5 논의 확장 시 재고).
- MLA-OLED 휘도효율 평가 기법. *Opt. Express* 12, 5777 (2004).
  — 사유: 기초 방법론, [3–5]로 대체 (미채택).
- 효율 ↔ 이미지 blur 트레이드오프 (MLA-film, J. Inf. Disp. 2018).
  — 사유: 공간 blur 축은 out of scope — v7은 angular 축만 다룸 (미채택).
- micro-Horn collimator µLED, arXiv:2412.14027.
  — 사유: collimator = 주기 굴절렌즈 클래스 밖; §2.6 aperture expansion 경로에
  예시가 필요해지면 재고 (보류).
- 메타표면 phosphor µLED 방향성. *ACS Nano* (2024).
  — 사유: 메타표면 조향 가드레일 역할은 [9] *Commun. Eng.* 2025로 대체 (미채택).

### 각도선택/재활용 관련 (부분 대체)
- "Broadband Angular Selectivity of Light at the Nanoscale." arXiv:1512.02761.
  — 사유: 각도선택+recycling 리뷰 — [8,14,15]로 커버; DBR 계산 §4.4 서술 보강에
  필요해지면 재고 (보류).
- "Methods and apparatus for broadband angular selectivity" (US 10073191).
  — 사유: 특허 인용 회피; TiO2/SiO2 다층막 선행 존재 사실만 내부 리스크 메모로 유지 (미채택).
- LCD 백라이트 BEF/DBEF 상용 선례 (정면휘도 ~1.6x).
  — 사유: 인용 가능한 1차 서지 없음 — 본문 언급 시 별도 출처 발굴 필요 (보류).

### 기타 (미채택/보류)
- 반구(H/D→0.5) 최적 microlens parameters. *Opt. Commun.* / *Org. Electron.* 계열
  (S1566119916303585). — 사유: hemisphere near-optimum 주장은 본 논문 자체 결과로
  제시 — 선행 인용 불필요 (미채택).
- "High-resolution beam steering using microlens arrays." *Opt. Lett.* 31, 2861 (2006).
  — 사유: 회절 프레임의 "주기 MLA 조향 불가" 논거는 P0 피벗(정리형 주장 철회) 후
  미사용 (미채택).
- 확장광원 freeform + étendue blur. *Optica* 3, 840 (2016).
  — 사유: v7 리스트 미포함이나 §2.5 lateral mixing/extended-source 논의와 관련 —
  리뷰 대응 시 재고 (보류).
- cross-scale 3D OLED 픽셀 EM 시뮬레이션. *Org. Electron.* (2022).
  — 사유: EM-스케일 시뮬은 본 방법론(CPS+ray tracing)과 다른 축, out of scope (미채택).
- Fluxim Setfos (상용 trans-scale 툴). https://www.fluxim.com/emission-module
  — 사유: 도구 웹페이지 — 참고문헌 아닌 방법 각주 후보 (보류).

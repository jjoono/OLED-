# 미확인·제외 문헌 목록 (2026-11-12)

`LIT_SURVEY_SUMMARY.md`와 `LIT_SURVEY_TABLE.csv`의 짝 문서다. 표에 넣지 않은 문헌을 모두 여기에 모았다.

| 구분 | 수 | 위치 |
|---|---|---|
| 미확인 후보: 문헌 | 36건 | 1.1절 |
| 미확인 후보: 특허 번호(원문 미열람) | 18건 | 1.2절 |
| 문헌은 표·제외 목록에 있으나 특정 수치·서술을 확인하지 못함 | 7건 | 1.3절 |
| 찾지 못함(없을 가능성) | 4건 | 1.4절 |
| Park & Suh 2018 피인용 문헌 | 21편(표 5, 제외 16) | 2절 |
| PDF needed(존재는 확인, 본문·초록은 못 봄) | 33편 | 3절 |
| 제외 문헌 전체(사유 포함) | 261편 | 4절 |

- 요약 맨 위의 "미확인 54건"은 1.1절(36건)과 1.2절(18건)을 더한 값이다.
- 미확인 후보는 표에 넣지 않았다. 서지는 검색 결과나 저장소 레코드에 보인 그대로 적었다. 기억으로 서지를 만들지 않았다.
- 1.1절의 내용은 각 주제 에이전트의 기록(`_work/*_unverified.md`)을 옮긴 것이다. 여러 주제에 나온 후보는 하나로 합쳤다.

---

## 1. 미확인 후보 (표에 넣지 않음)

### 1.1 문헌 36건
식별자(DOI, arXiv ID)가 없거나, 존재를 확인하지 못했거나, 내용을 원문에서 확인하지 못한 후보다. "주제·가설" 열의 내용은 검색 요약이나 초록에서 본 것이므로 인용하지 않는다.

| # | 후보 (보인 그대로) | 확인된 것 | 확인하지 못한 것 | 주제·가설 관련 | 다음 행동 |
|---|---|---|---|---|---|
| A-1 | D. Hosseini, *Chemically Installed Organic Metal Nucleation Inducers on Molecular Materials for Ultrathin Metal Electrodes*, Ph.D. dissertation, Loyola University Chicago, 2026 (ecommons.luc.edu/luc_diss/4289). 주제 E에서도 같은 후보를 찾음 | 기관 저장소 페이지. 초록만 열람 가능 | DOI 없음. 본문은 2028-02-06까지 엠바고 | A, H1. 초록 내용: 펜타센 위 o-MNI로 Ag 퍼콜레이션 11 → 8 nm. TPBi 위에서는 효과가 작고, 저자는 TPBi 자체의 핵 생성 경쟁 때문으로 해석. MADN 처리면 위 Ag 9 nm는 72%, 8.7 Ω/sq | 학술지 논문이 나오면 DOI로 교체 |
| A-2 | S. Jeong, GIST 학위논문(고분자 MNI 초박 Ag 전극), scholar.gist.ac.kr/handle/local/32859 | 검색 결과 요약 | DOI, 본문 | A. 작용기별 퍼콜레이션 수치가 있을 수 있음 | 표 [46] Jeong 2017 AFM과 대조 |
| A-3 | "Novel ultra-thin silver electrodes using an amine-containing nucleation inducer for organic solar cells", AP-HOPV18 학회 초록(nanoge) | 학회 초록 페이지 | DOI | A. 아민 고분자 seed, 기준선 약 15 nm | 동료심사 논문 아님 |
| A-4 | Bellchambers et al.(Hatton 그룹), 혼합 분자 단분자층 + Al 0.8 nm 하이브리드 seed(Cu 막, 퍼콜레이션 <4 nm) | 검색 요약 | DOI | A, H4. '분자층 + 극박 금속' 이중 seed 개념 | DOI 확인 |
| A-5 | "Sergeant et al.": Ag 5.5–6 Å/s, 기판 −5 °C가 섬 성장 최소의 최적점(내부 보고서 `litsearch_ultrathin_Ag/…report.md` 6.7절에 Lee & Lee 2015 리뷰 재인용으로 적힘) | 2차 인용 | 원문, DOI | A. 속도·냉각 조건 | 원문 확인 전 사용 금지 |
| B-1 | E. Najafabadi, *Stacked Inverted Top-Emitting White Organic Light-Emitting Diodes*, Ph.D. thesis, Georgia Tech, 2014, hdl:1853/52990. 주제 F에서도 같은 후보를 찾음 | 핸들이 저장소로 302 리다이렉트됨(존재). 검색 요약에 4.3절 "Novel Ag/HATCN Anode for Green Inverted Top-Emitting OLEDs"(p. 94) | 본문(저장소가 403 응답) | B, F, H3. Ag/HATCN 양극의 상세 데이터가 있을 수 있음 | 브라우저로 직접 열람 |
| B-2 | K. A. Knauer, *High-performance single-unit and stacked inverted top-emitting electrophosphorescent organic light-emitting diodes*, Georgia Tech thesis. 주제 F에서도 같은 후보를 찾음 | 저장소 레코드 URL | 연도, 본문, 식별자 | B, F. 같은 Kippelen 그룹의 Ag/HAT-CN 양극 | 직접 열람 |
| CD-1 | C. Falkenberg, *Optimizing Organic Solar Cells: Transparent Electron Transport Materials for Improving the Device Performance*, Ph.D. thesis, TU Dresden, 2012, urn:nbn:de:bsz:14-qucosa-89214 | 존재 확인. 전문 PDF를 읽음. 쪽수는 PDF 기준이다. 12 nm HATCN/유리는 RMS 2.4 ± 0.4 nm(p. 60–61, Fig. 3.2). 측면 연결은 10 nm 초과, 닫힌 막은 17 nm 초과(p. 61). n-도핑 HATCN 전도도 7.6e-6 S/cm(AOB 5 wt%), 2.2e-4 S/cm(NDN1 7 wt%)(p. 61). IP 9.57 eV, EA 4.6–5.1 eV(4.8 eV 채택), 광학 갭 3.3 eV(p. 59–64) | DOI 없음(규칙상 표에서 뺌) | B, D. H3 지지(유리 위라 조건이 다름). H4 관련 | 저널판(표 [131])에 같은 수치가 있는지 확인 |
| CD-2 | Lee, J.-H. et al., "Molecular structure of Hexaazatriphenylene-Hexanitrile (HATCN) Layer and its application to High Performance Transparent Inverted OLEDs", 한국고분자학회 학술대회 초록, 2013(OpenAlex W2346147194) | OpenAlex 레코드 | DOI, 원문 | B. 표 [97] Lee 2012와 같은 내용으로 추정 | — |
| CD-3 | 곽태호, 주성후, "CGL의 재료에 따른 청색 형광 Tandem OLED의 발광 특성", 2014(OpenAlex W2156474020) | OpenAlex 레코드 | DOI | C. 표 [126](10.5695/jkise.2014.47.4.210)과 같은 논문일 가능성이 높음 | 중복 여부 확인 |
| CD-4 | Li, L.-l. et al., "The Performance of Organic Solar Cells with Doping HATCN into CuPc as Anode Buffer Layer", Jilin Normal Univ. J., 2014(OpenAlex W2351864256) | OpenAlex 레코드 | DOI, 원문 | D. HATCN:CuPc 혼합 버퍼 | — |
| CD-5 | Kim, J. et al., "Study of Self-assembled Organic Layer Formation at the HATCN/Au Interface", 한국진공학회 초록, 2013(OpenAlex W1905222285) | OpenAlex 레코드 | DOI | E | — |
| CD-6 | Kang, H. et al., "Electron spectroscopy study of interface electronic structures of HATCN on Cu(111) surface", 한국진공학회 초록, 2009(OpenAlex W2242864065) | OpenAlex 레코드 | DOI | E | — |
| CD-7 | SID Digest P-203(연도·저자 미확인): Yb 도핑 n-CGL의 갭 상태와 HAT-CN LUMO의 상호작용으로 전하 생성(sid.org/Portals/sid/Files/SID-Digest/docs/P-203.pdf) | PDF URL | DOI(10.1002/sdtp.*) | C | DOI 확인 |
| CD-8 | Sims, 학위논문, 2014(Augsburg OPUS 40003): HATCN을 TPD에 1–10% 공증착 | 검색 요약 | 원문, DOI | D(혼합) | 원문 확인 |
| E-1 | Rösch et al.: 역구조 OPV에서 Ag가 MoO3로 확산해 MoO3 일함수가 바뀜 | JVSTA 2019 관련 검색 스니펫 | 논문 특정, DOI | E, H2 | — |
| E-2 | JVSTA 2019(Pisa 그룹 추정): 금속/MoO3 계면의 Mo 환원이 산화물 두께, 금속 반응성, 일함수 차이에 따라 달라짐 | arpi.unipi.it의 PDF 링크 | 저자, DOI | E, H2 | — |
| E-3 | AVS 2013 SS1-ThA8: F16CuPc/PFP 공여–수용 혼합막을 금속(Au, Cu, Ag(111)) 위에 올렸을 때의 코어 준위 이동 | 학회 초록 | DOI 없음 | E | — |
| E-4 | AVS 2007 SS1-MoM9(Heyde, Nilius, Pacchioni, Freund): MgO/Ag(001) 위 Au 대전의 STM 관찰 | 학회 초록 | DOI 없음 | E. 정식 논문은 표 [206] Sterrer 2007로 대체함 | — |
| E-5 | 2024년 열증착 Ag–TCNQ 박막: TCNQ0과 TCNQ− 사이의 전하 이동도 −0.5(C≡N 진동 이동에서 추정) | 검색 요약 | 논문, DOI | E | — |
| E-6 | AgTCNQ에서 TCNQ 진동 1455/2227 cm−1이 각각 73/17 cm−1 적색 이동 | 검색 요약 | 원 출처 | E | — |
| E-8 | Lee J. et al.(Yi Y. 그룹 추정), PCCP 2016: BCP–Ag 착물의 XPS | 기억에서 나온 단서뿐 | 존재, DOI | E | 존재부터 확인 |
| E-11 | Ag_n/Ti(0001)(n = 2–4): 클러스터가 커질수록 Ag 원자당 음전하가 줄어듦 | 검색 요약 | 논문, DOI | E, H2(Δq(n) 추세의 사례) | — |
| E-12 | "Greiner & Lu, TMO/금속 계면의 Auger parameter 분석"(Adv. Mater. 2013 또는 Anal. Bioanal. Chem. 2012로 제시됨) | Gemini가 제안함 | Crossref에서 찾지 못함. 존재 미확인 | E | 인용 금지 |
| G-1 | kMC of Ag nucleation on Si(111)-(7×7), arXiv cond-mat/0405277 | URL만 확인 | arXiv 메타데이터 | G(재구성 반도체 표면이라 관련성 낮음) | — |
| G-2 | Kabiraj & Mahapatra: 2D 전극의 Li 흡착을 고처리량 DFT와 crystal-graph ML로 스크리닝 | 검색 스니펫 | DOI | G(신규성 ②의 다른 분야 사례) | — |
| G-3 | Han et al.: Cu 단층/TaN(111)에서 liner 유무에 따른 DFT 접착으로 liner 적합성 판정 | 검색 스니펫 | Crossref 검색 실패 | G(신규성 ②) | — |
| G-4 | Simka et al.: DFT 접착 순서 TiN < Si 도핑 TiN < TaN < Ta가 100 Å Cu 젖음 실험과 일치 | 검색 스니펫 | Crossref 검색 실패(MRS Proceedings로 추정) | G(신규성 ②) | — |
| G-5 | ZnO/Ag/ZnO에서 Al–O 초박 seed층이 계면에너지를 낮춰 2D 성장을 유도(제1원리) | 검색 스니펫 | 제목, DOI | A, G(Yun 그룹 연작일 가능성) | 주제 A와 교차 확인 |
| G-6 | Starikova et al. 2007: hexaazatrinaphthylene 유도체의 IP/EA(실측 EA −2.65 ~ −4.59 eV, DFT) | 검색 스니펫 | DOI | H | 표 [345] Barlow 2007과 대조 |
| H-1 | Komatsubara & Yoshida, JSAP 2019 춘계 학술대회 11p-M112-6: LEIPS에 193 nm 대역통과 필터를 써서 EA > 6 eV를 측정(TCNQ, HAT-CN 포함) | 학회 초록 페이지 | DOI 없음, 수치 | H. HATCN EA의 확정값(> 4.9 eV 상한 해소)이 있을 수 있음 | 후속 논문 확인 |
| H-2 | Universidad Rey Juan Carlos 학위논문 "Novel synthetic approaches towards hexaazatriphenylene (HAT) derivatives"(burjcdigital.urjc.es handle 10115/11436) | 저장소 레코드 | 접속 실패, DOI 없음 | H(Segura 그룹 추정) | — |
| H-4 | "Linear and Nonlinear Optical Properties of Ramified Hexaazatriphenylenes: Charge Transfer Contributions to the Octupolar Response" | 다른 논문 참고문헌 목록의 제목 | DOI | H. 공여 치환 HAT의 CT 흡수 | — |
| H-5 | "Recent Progress of Hexaazatriphenylene-based Electrode Materials for Rechargeable Batteries"(NJTech 저장소 레코드) | 제목 | DOI | H(배터리 리뷰, 관련성 낮음) | — |
| H-6 | 2022년 J. Mater. Chem. C의 HATN계 순적색 TADF 논문 | 검색 요약 | 서지 | H. 확인된 HATNA TADF는 표 [327] Zhou 2019 | — |
| H-7 | Nasielski et al. 1990(10.1016/0040-4039(90)80128-9): "헥사아미노벤젠과 벤질을 축합하면 혼합물이 생긴다" | 논문 존재는 Crossref로 확인 | 서술 내용(검색 요약에서만 봄) | H(합성 경로) | 표 Kohne & Praefcke 1985 행의 note에 "(검색 요약, 미검증)"으로 표시함 |

### 1.2 특허 번호 18건 (번호만 확인, 원문 미열람)
검색 결과에서 번호와 요약만 봤다. 청구 범위와 실시예는 읽지 않았다.
- **Konica Minolta 투명전극 특허군(9건)**
  - 번호: US10355236B2, US9947889B2, US9917263B2, EP2890221, US11107998, US9966557, US9882072, US20150311452, EP2799226
  - 검색 요약상 내용: 질소 함유층/Ag, [n/M] 또는 ΔEef 지표, Al 중간층, 고표면에너지층 등
  - 같은 계열에서 원문을 확인한 것은 표에 넣었다: US9923166B2 [41], US9876189B2 [42], US9287521B2 [49], US10535834B2 [30]
- **OTI Lumionics NIC(핵 생성 억제 코팅) 특허군(8건)**
  - 번호: US11751415B2, US12178064B2, US11581487, US10700304, US11706969, US12590244, US12069939, US12101954
  - 원문을 확인한 US11730012B2(Example 4 kMC)는 표 [15]에 넣었다. US12492313B2는 그 계속출원이다.
- **기타(1건): US10693102.** 투명전극 발광소자에서 Al 도핑 LiF와 스퍼터 IZO 사이에 HAT-CN 20 nm 버퍼를 넣었더니 효율이 낮아졌다는 실시예가 있다고 한다(검색 요약).

### 1.3 문헌은 있으나 특정 수치·서술을 원문에서 확인하지 못한 항목 (인용 금지)

| # | 문헌 | 상태 | 확인하지 못한 수치·서술 |
|---|---|---|---|
| 1 | Najafabadi et al. 2013, Org. Electron., 10.1016/j.orgel.2013.02.024 | 제외, PDF needed | "Ag/HAT-CN anode … current densities up to three orders of magnitude higher than devices without HAT-CN"(검색 요약 문구) |
| 2 | Saragi et al. 2012, Synth. Met., 10.1016/j.synthmet.2012.07.021 | 제외, PDF needed | HATCN FET 이동도 2.5e-4 cm²/Vs(표 [88] 서론에서 재인용), "domain-like feature in nanoscale size"(검색 요약) |
| 3 | Liu et al. 2024, Mater. Sci. Semicond. Process., 10.1016/j.mssp.2023.108045 | 제외, PDF needed | "Liq/MoO3 CGL exhibit both lower electron and hole injection ability than Liq/HAT-CN"(검색 요약) |
| 4 | Cho, Pyo & Suh 2012, Synth. Met., 10.1016/j.synthmet.2011.12.028 | 제외, PDF needed | DNTPD 사이에 HAT-CN 약 5 nm를 넣어 p-도핑 HTL과 비슷한 전압을 얻음(검색 요약) |
| 5 | Hollerer et al. 2017, ACS Nano, 10.1021/acsnano.7b02449 | 제외 | 펜타센/Ag(001)에 MgO 초박막을 넣으면 흡착에너지가 2.65 → 4.71 eV가 된다는 수치가 이 논문 것인지 확인하지 못함 |
| 6 | Juárez et al. 2011, Chem. Eur. J., 10.1002/chem.201101198, 표 [337] | 표 | "tri-HAT is 0.75 V better acceptor than HAT"(URJC 저장소 초록 스니펫). 표의 key_numbers에는 넣지 않음 |
| 7 | Falkenberg et al. 2011, Sol. Energy Mater. Sol. Cells, 10.1016/j.solmat.2010.11.024, 표 [131] | 표(초록 근거) | 학위논문(1.1절 CD-1)의 형태 수치(RMS, 측면 연결 두께)가 저널판에도 있는지 |

### 1.4 찾지 못한 것 (없을 가능성)
- **Ag(I)–HAT(CN)6 배위 고분자나 결정 구조.** 검색했지만 보고를 찾지 못했다. CCDC를 직접 검색해야 한다.
- **HATCN 박막의 n, k 또는 흡수 시작 파장을 직접 보고한 논문.** 찾지 못했다. 팀의 자체 nk(`nk_JH_total.mat`의 HATCN)가 지금으로서는 유일한 정량 자료다.
  - 확인 후보: Opt. Mater. X 2025의 OLED 박막 19종 광학 상수 논문(10.1016/j.omx.2025.100419, 3절)
  - 초록에서는 HATCN이 들어 있는지 알 수 없었다.
- **Google Scholar 기준 앵커 피인용 목록.** 봇 차단 때문에 직접 조회하지 못했다. Crossref, OpenAlex, OpenCitations COCI, Semantic Scholar의 합집합으로 대신했다(2절).
- **KIPRIS(한국 특허)와 CNKI(중국 학위논문).** 접근하지 못했다. 한국어·중국어 웹 검색에서는 HATCN seed/Ag 관련 학위논문이나 특허를 찾지 못했다.

---

## 2. Park & Suh 2018(Opt. Express 26, 4979) 피인용 문헌 처리 결과

- **피인용 출처**
  - 출처별 수: Crossref is-referenced-by 20편, OpenAlex 20편, OpenCitations COCI 20편, Semantic Scholar 16편
  - 합집합은 21편이다. COCI에만 있는 것이 10.1002/sdtp.18677이고, OpenAlex와 S2에만 있는 것이 10.1109/eltech.2019.8839348이다.
- **결론:** 21편 모두 앵커를 '얇은 Ag 상부 전극'이나 소자 구조의 일반 사례로 인용했다. **HATCN 결정화의 근거로 인용한 논문은 없다.**

### 2.1 표에 넣은 피인용 문헌 (5편)
| id | 문헌 | DOI | 표에 넣은 이유 |
|---|---|---|---|
| [77] | Kwon et al. 2021, J. Ind. Eng. Chem. | 10.1016/j.jiec.2020.09.029 | Ag 위 HATCN 캡핑층의 형태(Rq), 접촉각, XPS |
| [11] | Le et al. 2024, Org. Electron. | 10.1016/j.orgel.2024.107061 | EIL–Ag 결합에너지 계산값(Table 1) |
| [12] | Ok et al. 2024, Adv. Opt. Mater. | 10.1002/adom.202402832 | Al/Cu 버퍼, 85 °C 240 h 음극 안정성 |
| [7] | Ok et al. 2025, SID Symp. Dig. Tech. Pap. | 10.1002/sdtp.18677 | LiF/Al/Ag 결합에너지와 형태 |
| [230] | Riahi et al. 2023, Adv. Electron. Mater. | 10.1002/aelm.202300675 | Ag/WO3 다층 상부 전극(대면적 균일도) |

### 2.2 제외한 피인용 문헌 (16편)
사유는 4절의 B 표에도 같이 있다.

| # | DOI | 연도 | 제목 | 제외 사유 |
|---|---|---|---|---|
| 1 | 10.1016/j.nxnano.2025.100132 | 2025 | OLED-on-silicon (OLEDoS) microdisplays: Technology challenges, design considerations, and adaptation in eXtended Reality (XR) ecosystem – Review | OLEDoS 리뷰. 앵커를 일반 사례로 인용. HATCN/Ag 1차 데이터 없음(초록 수준) |
| 2 | 10.1021/acsami.5c11689 | 2025 | Thermal Cross-linked Electron Transport Polymers for Suppressing Efficiency Roll-off in Green Solution-Processed Inverted OLEDs | 용액공정 역구조 OLED용 ETL 고분자. HATCN 형태, Ag 성장 데이터 없음(초록) |
| 3 | 10.1016/j.mseb.2024.117271 | 2024 | Conjugated polymeric nanocomposite-based light-generating active materials for OLED applications: A review | 고분자 나노복합 발광체 리뷰. HATCN, 박막 금속 성장과 무관 |
| 4 | 10.1016/j.mseb.2023.117118 | 2024 | Solution processed PANI-rGO-SiO2 nanocomposite with inflated electron injection and transport capabilities as ETL for OLED applications | 용액공정 나노복합 ETL. HATCN, 박막 금속 성장과 무관 |
| 5 | 10.1002/adfm.202214530 | 2023 | Highly Efficient Top‐Emitting Infrared‐to‐Visible Up‐Conversion Device Enabled by Microcavity Effect | Mg:Ag 상부 음극을 쓴 전면발광 상향변환 소자. HATCN seed, 형태 데이터 없음(초록) |
| 6 | 10.3390/coatings12101583 | 2022 | Transparent Conductive Indium Zinc Oxide Films: Temperature and Oxygen Dependences of the Electrical and Optical Properties | 스퍼터 IZO 투명전도막. HATCN, 박막 Ag와 무관 |
| 7 | 10.1016/j.vacuum.2022.111481 | 2022 | Nitrogen-modified a-C:H layers applied for light emitting diodes | a-C:H 층. HATCN, 박막 Ag와 무관 |
| 8 | 10.1007/s11664-021-09338-0 | 2021 | Organic Devices: Fabrication, Applications, and Challenges | 유기 소자 전반 리뷰. HATCN, 박막 Ag 성장의 1차 데이터 없음 |
| 9 | 10.1016/j.mtener.2021.100704 | 2021 | Fabrication of an oxide/metal/oxide structured electrode integrated with antireflective film to enhance performance in flexible organic light-emitting diodes | 스퍼터 ZnO/AgOx/ZnO OMO(원자 산소 매개 Ag). HATCN 없음. 앵커는 일반 인용(로컬 전문 확인). 주제 A(산소 매개 Ag 성장)와는 관련 가능 |
| 10 | 10.1007/s11664-021-08776-0 | 2021 | Realizing 8 cd A−1 Current Efficiency for Solution-Processed Inverted Top-Emitting Polymer Light-Emitting Diodes | 용액공정 역구조 전면발광 PLED. 초록이 없고 HATCN 형태나 Ag 성장 데이터의 단서도 없음. 유료 |
| 11 | 10.1016/j.rser.2020.110043 | 2020 | Organic light emitting diode devices: An energy efficient solid state lighting for applications | OLED 일반 리뷰. HATCN 형태나 금속 성장 데이터 없음(S2/OpenAlex/COCI 피인용) |
| 12 | 10.1002/adom.201901421 | 2020 | Microcavity‐Enhanced Blue Organic Light‐Emitting Diode for High‐Quality Monochromatic Light Source with Nonquarterwave Structural Design | DBR 미세공진 청색 OLED. 앵커는 미세공진 맥락으로 인용. HATCN/Ag 성장 데이터 없음(초록) |
| 13 | 10.1016/j.orgel.2019.105611 | 2020 | Determination of emitting dipole orientation in organic light emitting diodes | 발광 쌍극자 배향 측정법. 무관 |
| 14 | 10.1109/eltech.2019.8839348 | 2019 | Light Outcoupling Integrated Photonic Crystals for OLEDs | 광결정 광추출 학회 논문. 무관(OpenAlex/S2 피인용, COCI에는 없음) |
| 15 | 10.1016/j.orgel.2018.08.051 | 2018 | Tailoring Mg:Ag functionalities for organic light-emitting diodes | Mg:Ag 음극 연구. 로컬 전문에 HATCN 언급 0회. 앵커는 Ag Rs·두께 맥락으로 인용. 주제 A/F 후보(Mg:Ag) |
| 16 | 10.1016/j.orgel.2018.09.009 | 2018 | Thermally stable benzo[f]quinoline based bipolar host materials for green phosphorescent OLEDs | 호스트 재료 합성(Suh 그룹). 앵커는 소자 구조 맥락으로 인용. HATCN/Ag 형태 없음 |

### 2.3 앵커의 역방향 참고문헌 중 HATCN/Ag 관련 항목
- [25] Jeon 2012, Org. Electron., 10.1016/j.orgel.2012.02.008: 제외. 전기 특성만 다룬다.
- [26] Cho, Pyo & Suh 2012, Synth. Met., 10.1016/j.synthmet.2011.12.028: 제외, PDF needed(3절).
- [27] Najafabadi 2013, Org. Electron., 10.1016/j.orgel.2013.02.024: 제외, PDF needed(3절).
- [28], [29]: ZnS/Ag/MoO3 전극과 OPV용 얇은 금속 전극. HATCN과는 무관하고 주제 A에서 다뤘다.
- 추가로 HATCN 관련 핵심 문헌 7편(Frank 2010, Kim 2013, Yamaguchi 2020, Lapkin 2026, Opitz 2026, Joo 2023, Bröker 2010)의 피인용도 추적했다. 이 추적으로 Lapkin 2026(표 [73])과 Bröker 학위논문(10.18452/16257, 제외·통합)을 찾았다.

---

## 3. PDF needed — 열면 표에 넣을 수 있는 문헌 (33편)

존재와 DOI는 확인했지만 본문이나 초록을 열지 못해 표에서 뺀 문헌이다. 대부분 Elsevier 본문이다. 우선순위는 가설 판정과 신규성 판단에 주는 영향으로 정했다. 링크는 DOI 주소다.

### 3.1 우선순위 상 (14편)
| # | 문헌 | 연도 | 제목 | 왜 필요한가 |
|---|---|---|---|---|
| 1 | [10.1016/j.orgel.2013.02.024](https://doi.org/10.1016/j.orgel.2013.02.024) | 2013 | High-performance inverted top-emitting green electrophosphorescent organic light-emitting diodes with a modified top Ag anode | Najafabadi 2013. HAT-CN/Ag 상부 양극(Kippelen 그룹), 앵커 역방향 [27]. Ag 두께·형태 데이터가 있으면 H3의 비교군 |
| 2 | [10.1016/j.synthmet.2011.12.028](https://doi.org/10.1016/j.synthmet.2011.12.028) | 2012 | Low voltage top-emitting organic light emitting devices by using 1,4,5,8,9,11-hexaazatriphenylene-hexacarbonitrile | Cho, Pyo & Suh 2012. Suh 그룹 HATCN TEOLED, 앵커 역방향 [26] |
| 3 | [10.1016/j.synthmet.2012.07.021](https://doi.org/10.1016/j.synthmet.2012.07.021) | 2012 | Electron mobility in hexaazatriphenylene hexacarbonitrile field-effect transistors | Saragi 2012. HATCN FET 이동도와 형태(요약 3절 모순 #6) |
| 4 | [10.1016/j.tsf.2009.05.011](https://doi.org/10.1016/j.tsf.2009.05.011) | 2009 | Efficient hole injection by doping of hexaazatriphenylene hexacarbonitrile in hole transport layer | HAT-CN 도핑 HTL. H4(혼합·희석) |
| 5 | [10.1016/j.orgel.2012.06.039](https://doi.org/10.1016/j.orgel.2012.06.039) | 2012 | Determination of the interface energy level alignment of a doped organic hetero-junction using capacitance–voltage measurements | Kim J.-J. 그룹. C–V로 도핑 이종접합(n-ETL/HATCN 추정)의 계면 정렬 결정. H4 |
| 6 | [10.1016/j.orgel.2021.106202](https://doi.org/10.1016/j.orgel.2021.106202) | 2021 | Efficient inverted organic light-emitting devices using a charge-generation unit as electron-injection layers | 역구조 OLED에서 CGU를 전자주입층으로 사용. 팀 구조와 같은 물리 |
| 7 | [10.1016/j.mssp.2023.108045](https://doi.org/10.1016/j.mssp.2023.108045) | 2024 | Effects of charge generation layers on multiple guest/host bilayer-based tandem OLEDs | Liu 2024. Liq/HAT-CN 대 Liq/MoO3 CGL 비교 |
| 8 | [10.1016/S0167-5729(96)00011-8](https://doi.org/10.1016/S0167-5729(96)00011-8) | 1997 | Ultrathin metal films and particles on oxide surfaces: structural, electronic and chemisorptive properties | Campbell 1997 총설. 금속/산화물 성장 모드와 전하 이동의 표준 인용처 |
| 9 | [10.1021/acsaelm.0c00156](https://doi.org/10.1021/acsaelm.0c00156) | 2020 | Ultrathin Ag Transparent Conducting Electrode Structure for Next-Generation Optoelectronic Applications | Kwon 그룹. 유기 HIL 사이 Ag 9 nm(Rs 17.6 Ω/sq, T 79.2%). HIL이 HAT-CN이면 HATCN 위 Ag ≤10 nm의 유일한 데이터 |
| 10 | [10.1016/j.orgel.2014.08.004](https://doi.org/10.1016/j.orgel.2014.08.004) | 2014 | Highly efficient inverted top emitting organic light emitting diodes using a horizontally oriented green phosphorescent emitter | Kim 2014. 역구조 TEOLED. 상부 전극이 HAT-CN/Ag인지 확인 |
| 11 | [10.1016/j.orgel.2013.06.004](https://doi.org/10.1016/j.orgel.2013.06.004) | 2013 | Stacked inverted top-emitting green electrophosphorescent organic light-emitting diodes on glass and flexible glass substrates | Knauer 2013. Kippelen 그룹 적층 역구조 TEOLED의 상부 Ag 양극 구조 |
| 12 | [10.1109/IITC.2018.8430425](https://doi.org/10.1109/IITC.2018.8430425) | 2018 | A First-Principles Density Functional Theory Based Framework for Barrier Material Screening | Cu 배선 배리어 DFT 스크리닝 프레임워크. 신규성 ②. OpenAlex 초록이 다른 논문과 섞여 있음 |
| 13 | [10.1016/j.ensm.2023.103035](https://doi.org/10.1016/j.ensm.2023.103035) | 2023 | A systematic study on the metallophilicity of ordered five-atomic-layer MXenes using high-throughput automated workflow and machine learning | 금속 음극용 고처리량 DFT+ML 금속친화도 스크리닝. 신규성 ②의 다른 분야 선례 |
| 14 | [10.1016/j.apsusc.2009.11.044](https://doi.org/10.1016/j.apsusc.2009.11.044) | 2010 | Interaction of bathocuproine with metals (Ca, Mg, Al, Ag, and Au) studied by density functional theory | BCP–금속(Ag 포함) DFT. 팀 ETL E_b와 직접 비교 대상 |

### 3.2 우선순위 중 (16편)
| # | 문헌 | 연도 | 제목 | 왜 필요한가 |
|---|---|---|---|---|
| 1 | [10.1007/s12274-024-6718-y](https://doi.org/10.1007/s12274-024-6718-y) | 2024 | High-performance ultrathin Ag electrodes by chemical bond anchoring Ag atoms for stretchable organic light-emitting devices | 화학결합 앵커링 초박막 Ag(Nano Research 2024). 주제 A |
| 2 | [10.1016/j.vacuum.2017.07.007](https://doi.org/10.1016/j.vacuum.2017.07.007) | 2017 | Study of injection and transport properties of metal/organic interface using HAT-CN molecules as hole injection layer | 금속/HAT-CN 주입 특성 |
| 3 | [10.1016/j.orgel.2009.05.007](https://doi.org/10.1016/j.orgel.2009.05.007) | 2009 | P-type doping of organic wide band gap materials by transition metal oxides: A case-study on Molybdenum trioxide | MoO3를 넓은 갭 호스트에 도핑(희석)한 고전 연구. H4 |
| 4 | [10.1016/j.orgel.2013.11.033](https://doi.org/10.1016/j.orgel.2013.11.033) | 2014 | Highly efficient p-dopants in amorphous hosts | 비정질 호스트 속 p-도펀트 효율. H4 |
| 5 | [10.1016/j.orgel.2021.106353](https://doi.org/10.1016/j.orgel.2021.106353) | 2022 | Lateral current suppression in tandem organic light-emitting diodes by adopting a buffer layer | 탠덤 CGL 측면 전류 억제 버퍼(HATCN 측면 전도와 관련 가능) |
| 6 | [10.1016/j.orgel.2014.10.015](https://doi.org/10.1016/j.orgel.2014.10.015) | 2014 | Role of hole injection layer in intermediate connector of tandem organic light-emitting devices | 탠덤 연결층 안 HIL(HATCN 등)의 역할 |
| 7 | [10.1016/j.orgel.2015.05.033](https://doi.org/10.1016/j.orgel.2015.05.033) | 2015 | The role of gap states on energy level alignment at an α-NPD/HAT(CN) 6 charge generation interface | α-NPD/HAT(CN)6 전하생성 계면의 갭 상태(p측) |
| 8 | [10.1016/j.orgel.2020.105979](https://doi.org/10.1016/j.orgel.2020.105979) | 2021 | Investigation on the mechanism of charge generation in organic heterojunctions: Analysis of I–V and C–V characteristics | 유기 이종접합 전하생성 메커니즘(I–V, C–V) |
| 9 | [10.1016/j.optmat.2022.112345](https://doi.org/10.1016/j.optmat.2022.112345) | 2022 | Interfacial charge transfer study of hexacarbonitrile-based intermediate connecter in blue tandem organic light-emitting diodes | 헥사카보니트릴 연결층의 계면 전하 이동(청색 탠덤) |
| 10 | [10.1016/j.orgel.2017.10.042](https://doi.org/10.1016/j.orgel.2017.10.042) | 2018 | Highly efficient tandem organic light-emitting devices adopting a nondoped charge-generation unit and ultrathin emitting layers | 비도핑 CGU 탠덤(LiF/Al/HAT-CN 추정) |
| 11 | [10.1016/j.orgel.2017.11.022](https://doi.org/10.1016/j.orgel.2017.11.022) | 2018 | Highly-efficient tandem organic light-emitting device employing bis-4,6-(3,5-di-3-pyridylphenyl)-2-methylpyrimi-dine doped with cesium azide in charge generation unit | B3PYMPM:CsN3 n-CGL |
| 12 | [10.1016/j.jiec.2025.06.044](https://doi.org/10.1016/j.jiec.2025.06.044) | 2026 | Achieving highly stable tandem organic light-emitting diodes using interfacial property control of charge generation layer | CGL 계면 제어로 탠덤 안정화(n-CGL/HAT-CN 계면일 가능성) |
| 13 | [10.1039/c4tc02082k](https://doi.org/10.1039/c4tc02082k) | 2014 | A novel intermediate connector with improved charge generation and separation for large-area tandem white organic lighting devices | Ding 2014. 전하 생성·분리를 개선한 연결층. RSC 본문 403 |
| 14 | [10.1016/j.surfrep.2007.07.001](https://doi.org/10.1016/j.surfrep.2007.07.001) | 2007 | Interaction of nanostructured metal overlayers with oxide surfaces | Fu & Wagner 2007 총설. 금속 나노구조/산화물 상호작용 |
| 15 | [10.1039/d6el00029k](https://doi.org/10.1039/d6el00029k) | 2026 | Energetics and chemistry at electron selective interfaces for p-i-n perovskite solar cells: an in situ investigation | 전자 선택 계면 총론. 2차 스니펫상 BCP 위 Ag 착물(N 1s 새 성분) |
| 16 | [10.1007/BF01112848](https://doi.org/10.1007/BF01112848) | 1991 | Energy-adjusted ab initio pseudopotentials for the second and third row transition elements: Molecular test for M2 (M=Ag, Au) and MH (M=Ru, Os) | Ag2 분자로 시험한 ECP. 팀의 Ag2 1.86 eV 검증과 관련 |

### 3.3 우선순위 하 (3편)
| # | 문헌 | 연도 | 제목 | 왜 필요한가 |
|---|---|---|---|---|
| 1 | [10.1016/j.ssc.2004.04.007](https://doi.org/10.1016/j.ssc.2004.04.007) | 2004 | The effect of heat treatment on bistable Ag-TCNQ thin films | Ag–TCNQ 박막 열처리. 스니펫상 N 1s 새 상태(Ag → TCNQ 전자 이동) |
| 2 | [10.1016/j.tetlet.2012.01.137](https://doi.org/10.1016/j.tetlet.2012.01.137) | 2012 | A hexaazatriphenylene-based organogel that responds to silver(I) with high selectivity under aqueous condition | Ag(I)에 선택 반응하는 HAT 유기겔 |
| 3 | [10.1016/j.tet.2004.01.083](https://doi.org/10.1016/j.tet.2004.01.083) | 2004 | Hexaazatriisothianaphthenes: new electron-transport mesogens? | 헥사아자트리이소티아나프텐(티에노 축합 HAT). 표 [345] Barlow 2007과 연결 |

### 3.4 확인 권장 (PDF needed는 아니지만 열어 볼 가치가 있음)
- [10.1016/j.omx.2025.100419](https://doi.org/10.1016/j.omx.2025.100419) (2025) "Comprehensive optical characterization of organic thin films for OLED …". OLED 막 19종의 광학 상수를 다룬다. HATCN이 포함됐다면 문헌 공백 11(HATCN n, k)을 메울 수 있다.
- [10.1021/acsaelm.0c00156](https://doi.org/10.1021/acsaelm.0c00156)은 3.1에 있다. HIL이 HAT-CN인지만 확인하면 바로 표에 넣을 수 있다.

---

## 4. 제외 문헌 전체 (261편, 사유 포함)

- 주제 에이전트별로 묶었다. 한 문헌을 여러 에이전트가 제외했으면 처음 기록한 에이전트 아래에 한 번만 적었다.
- 사유가 "대표 행 …에 통합"인 것은 같은 그룹의 중복 결과다. 수치는 대표 행의 `note`에 하위 항목으로 넣었다.
- "[Park&Suh 2018 피인용]"이 붙은 16편은 앵커 피인용 문헌이다(2.2절).
- 사유가 영어로만 적혀 있던 B 에이전트의 66건은 한국어로 옮겼다. 원문 사유는 `_work/B_excluded.jsonl`에 있다.

### 4.1 A. seed 위 초박막 Ag (41편)
| # | DOI | 연도 | 제목 | 제외 사유 |
|---|---|---|---|---|
| A-1 | 10.1039/d5mh02144h | 2026 | High-performance transparent metal mesh electrodes utilizing a metal-vapor-desorption layer for organic light-emitting diode applications | 금속 증기 탈착층으로 Ag 메시 패터닝(2026) — 메시는 범위 밖, 응축계수 개념 사례로만 노트에 언급 |
| A-2 | 10.1007/s40820-026-02083-1 | 2026 | Enhancing π-Delocalization and Suppressing Traps via Doping in Electron Transport Materials for Efficient Semitransparent Organic Photovoltaics | PDINN vs D-PDINN 위 Ag 12 nm(OSC) — 소자 T만, percolation 비교 없음(같은 계열 10.1002/smll.75633 행 참조) |
| A-3 | 10.1007/s40820-026-02330-5 | 2026 | Polar-Twisted Electron-Transport Layer Simultaneously Breaks Efficiency-Flexibility-Cost Limits in Organic Solar Cells | OSC ETL 위 Ag 20 nm — 얇은 Ag 성장 정보 없음 |
| A-4 | 10.1016/j.vacuum.2025.114958 | 2026 | Oxygen-modulated wettability and morphology of Ag thin films on ZnO: A molecular dynamics study | ZnO 위 산소 조절 Ag 젖음 MD 연구(2026) — topic G(시뮬레이션)로 이관 권장 |
| A-5 | 10.1002/advs.202410931 | 2025 | α‐ZnO Manipulated Growth of Ag Gird on AgNWs Enables High Conductive Flexible Electrode for Large‐Area Monolithic Organic Photovoltaics | α-ZnO/AgNW 하이브리드 그리드 — 균일 초박막 아님 |
| A-6 | 10.1016/j.vacuum.2024.113949 | 2025 | Study of the optimal deposition conditions of Ag(O2) and Ag(N2) seed layers on ultrathin silver transparent conductive films | Ag(O2)/Ag(N2) 스퍼터 seed(percolation 3.9 nm) — 스퍼터 전용; Ag(O) 효과는 10.1016/j.apsusc.2025.165195, 10.1016/j.apsusc.2022.154421로 대표 |
| A-7 | 10.1016/j.vacuum.2025.114630 | 2025 | ZnS-enabled ultrathin Ag layers with high transparency and conductivity for flexible optoelectronic devices | ZnS 1 nm 스퍼터 seed(7.5 nm) — ZnS 근거는 10.1002/adfm.201502542로 대표 |
| A-8 | 10.1016/j.apsusc.2025.162342 | 2025 | Exploring SiO  as an effective adhesion promoter for Ag on glass and polymer substrates | SiO 부착 촉진층(2025) — 접착 중심, 초록 미확인 |
| A-9 | 10.1016/j.cej.2024.150839 | 2024 | Optimizing Ag films towards efficient flexible quantum-dot light-emitting diodes | PEI seed(QLED, Ag 16 nm) — PEI 근거 중복, 얇은 Ag 아님 |
| A-10 | 10.1016/j.cej.2024.156475 | 2024 | Mitigating interfacial recombination enabling efficient semitransparent organic photovoltaics | 11-MUA SAM(OPV, Ag 12 nm) — 소자 T만, percolation 정보 없음 |
| A-11 | 10.1007/s12274-024-6718-y | 2024 | High-performance ultrathin Ag electrodes by chemical bond anchoring Ag atoms for stretchable organic light-emitting devices | 화학결합 앵커링 초박막 Ag(Nano Research 2024) — 초록·본문 접근 불가(로그인 리다이렉트). PDF needed |
| A-12 | 10.37188/cjl.20230090 | 2023 | 基于超薄透明Ag2O/Ag阳极的高效有机电致发光器件 | 같은 그룹 Ag2O seed 중복 결과 — 대표 행 10.1016/j.apsusc.2022.154421 note에 통합 |
| A-13 | 10.1016/j.mssp.2022.107188 | 2023 | a-IGZO thin-film transistors with transparent ultrathin Al/Ag bilayer source and drain for active neural interfaces | Al/Ag TFT 전극 — Al seed 근거 중복 |
| A-14 | 10.1016/j.solener.2022.12.047 | 2023 | Transparent-conductive and infrared-shielding WO3/Ag/WO3 multilayer heterostructures | WO3/Ag/WO3, Ag ≥16 nm만 — ≤10 nm 정보 없음(같은 두께에서 WO3가 유리보다 Rs 높음, DB) |
| A-15 | 10.1016/j.solener.2021.09.065 | 2021 | Bifacial semitransparent perovskite solar cells with MoOx/Cu/Ag/MoOx multilayer transparent electrode | MoOx/Cu/Ag PSC 전극 — Cu seed 근거 중복 |
| A-16 | 10.1116/6.0000244 | 2020 | Manipulation of thin silver film growth on weakly interacting silicon dioxide substrates using oxygen as a surfactant | 같은 그룹(Sarakinos) — 대표 행 10.1021/acsanm.0c00736 note에 통합(O2 계면활성제) |
| A-17 | 10.1016/j.apsusc.2020.146989 | 2020 | Minimizing optical loss in ultrathin Ag films based on Ge wetting layer: Insights on Ge-mediated Ag growth | 같은 그룹(Yun) Ge seed — 대표 행 10.1016/j.surfin.2024.105208 note에 통합 |
| A-18 | 10.1016/j.cap.2020.06.016 | 2020 | Investigation of ultra-thin and flexible Au–Ag–Au transparent conducting electrode | Au 2/Ag 4/Au 2 삼중층을 100 °C 기판에서 증착 — Ag 단독 percolation 분리 불가 |
| A-19 | 10.1016/j.solener.2020.06.013 | 2020 | Semitransparent perovskite solar cells with ultrathin silver electrodes for tandem solar cells | Cu 1 nm seed(PSC) — Cu seed 근거(10.1126/sciadv.abg0355, 10.1021/am303147w)와 같은 방향, 새 정보 적음 |
| A-20 | 10.1016/j.surfcoat.2020.125752 | 2020 | Influence of aluminum interlayer on optical properties of very thin silver thin film | Al 1 nm/Ag 6 nm(16.5 Ω/sq) — Al seed 근거 중복(10.1016/j.orgel.2016.03.002 등) |
| A-21 | 10.1002/adfm.202005959 | 2020 | Transparent Fused Nanowire Electrodes by Condensation Coefficient Modulation | 응축계수 조절 Ag 나노와이어 망 — 메시/NW는 범위 밖(개념만 노트에 인용) |
| A-22 | 10.1021/acsami.9b08594 | 2019 | Correlating Nanostructure, Optical and Electronic Properties of Nanogranular Silver Layers during Polymer-Template-Assisted Sputter Deposition | 같은 그룹(Roth/Schwartzkopf) — 초록만, 수치 없음; 10.1021/acsami.5b02901 note에 통합 |
| A-23 | 10.1016/j.mtcomm.2019.100592 | 2019 | Ultrathin continuous silver film for efficient infrared emission | 플라스마 처리 Ge seed, 적외선 방사 목적 — OLED 상부전극 관련성 낮음 |
| A-24 | 10.1016/j.apsusc.2019.04.043 | 2019 | Degradation mechanism of protected ultrathin silver films and the effect of the seed layer | MPTMS/ZnO/SiO2 seed, 스퍼터 Ag — 열화 메커니즘 중심 |
| A-25 | 10.1021/acsaem.9b00878 | 2019 | Stabilizing Silver Window Electrodes for Organic Photovoltaics Using a Mercaptosilane Monolayer | MPTMS 단분자층으로 Ag 13 nm 안정화(Hatton 2019) — 성장/percolation이 아닌 안정성 중심 |
| A-26 | 10.1002/admi.201900608 | 2019 | An Oxide‐Dispersed Preparation Strategy for Silver Ultrathin Films with Low Percolation Threshold Thickness, Subnanometer Smoothness, and Prominent Durability | 산화물 분산 Ag 초박막(2019) — 공정 기반(seed 아님), 후속 확인 후보 |
| A-27 | 10.1515/nanoph-2020-0214 | 2018 | Highly transparent and conductive metal oxide/metal/polymer composite electrodes for high‐efficiency flexible organic light‐emitting devices | PAI 고분자 seed(용액) — 고분자 seed 근거 중복 |
| A-28 | 10.1063/1.4963262 | 2016 | Continuous ultrathin silver films deposited on SiO2 and SiNx using a self-assembled monolayer | SiO2/SiNx 위 SAM으로 연속 초박 Ag(APL 2016) — 초록 미확보, 스퍼터 추정; 후속 확인 후보 |
| A-29 | 10.1021/acs.langmuir.5b01244 | 2015 | Phenomenological Model of the Growth of Ultrasmooth Silver Thin Films Deposited with a Germanium Nucleation Layer | Ge seed 성장 모델 — 대표 행 10.1021/nl8027476 note에 통합(초록만) |
| A-30 | 10.1021/acsami.5b01471 | 2015 | Ge Wetting Layer Increases Ohmic Plasmon Losses in Ag Film Due to Segregation | Ge 젖음층의 Ag 내 편석이 손실 증가 — Ge 부적합 근거로 노트에만 언급 |
| A-31 | 10.1016/j.orgel.2014.06.021 | 2014 | Flexible organic light-emitting diodes with ZnS/Ag/ZnO/Ag/WO3 multilayer electrode as a transparent anode | ZnS/Ag/ZnO/Ag/WO3 양극 — ALD ZnO 위 18 nm 등 두꺼운 조건, 관련성 낮음 |
| A-32 | 10.1002/adma.201306212 | 2014 | Interfacial Engineering of Ultrathin Metal Film Transparent Electrode for Flexible Organic Photovoltaic Cells | 초박 금속막 계면공학(Adv Mater 2014, 피인용 ~200) — 초록 미확보, 후속 확인 후보 |
| A-33 | 10.1016/j.orgel.2013.09.014 | 2013 | ITO-free flexible organic light-emitting diode using ZnS/Ag/MoO3 anode incorporating a quasi-perfect Ag thin film | 같은 그룹 중복(KAIST Choi) — 대표 행 10.1002/adfm.201502542의 note에 수치 통합 |
| A-34 | 10.1002/adfm.201301359 | 2013 | Transparent Ultrathin Oxygen‐Doped Silver Electrodes for Flexible Organic Solar Cells | 산소 도핑 Ag(공정 기반, seed 아님) — Ag(O) 효과는 Yun 그룹 행으로 대표 |
| A-35 | 10.1002/adfm.201002021 | 2011 | Ultrathin Transparent Au Electrodes for Organic Photovoltaics Fabricated Using a Mixed Mono‐Molecular Nucleation Layer | MPTMS/APTMS 혼합 단분자층 위 Au(Ag 아님) |
| A-36 | 10.1364/oe.18.005124 | 2010 | Ultra-thin ultra-smooth and low-loss silver films on a germanium wetting layer | Ge seed 계열 — 대표 행 10.1021/nl8027476 note에 통합(초록만) |
| A-37 | 10.1021/nn100466p | 2010 | Enhanced Surface Plasmon Resonance on a Smooth Silver Film with a Seed Growth Layer | seed 성장층 위 평활 Ag의 SPR(ACS Nano 2010) — 플라즈모닉 목적, 두꺼운 Ag |
| A-38 | 10.1016/j.orgel.2010.03.022 | 2010 | E-beam deposited ultra-smooth silver thin film on glass with different nucleation layers: An optimization study for OLED micro-cavity application | e-beam Ag와 여러 핵생성층(OLED 미세공진용, 2010) — 후보로 남김, 본문 미확인 |
| A-39 | 10.1016/j.tsf.2007.05.080 | 2008 | Coalescence of magnetron-sputtered silver islands affected by transition metal seeding (Ni, Cr, Nb, Zr, Mo, W, Ta) and other parameters | 같은 그룹 중복(LBNL Anders) — 대표 행 10.1016/j.ssc.2006.08.027 note에 통합 |
| A-40 | 10.1016/s0039-6028(00)00234-x | 2000 | Determination of condensation coefficients of metals on polymer surfaces | 같은 그룹(Faupel) — 대표 행 10.1103/physrevlett.82.1903 note에 통합 |
| A-41 | 10.1016/s0169-4332(98)00826-5 | 1999 | Condensation coefficients and initial stages of growth for noble metals deposited onto chemically different polymer surfaces | 같은 그룹(Faupel) — 대표 행 10.1103/physrevlett.82.1903 note에 통합(검색 색인 초록 수준) |

### 4.2 B. HATCN 형태·피인용 추적 (66편)
| # | DOI | 연도 | 제목 | 제외 사유 |
|---|---|---|---|---|
| B-1 | 10.1021/acs.chemmater.5c02419 | 2026 | Unveiling the Role of Solvents in DBTTF:HATCN Ternary Cocrystals | 용액에서 키운 벌크 공결정(CT 0.11 e). seed층과 관련된 박막 형태 없음 |
| B-2 | 10.48550/arxiv.2610.34543 | 2026 | Role of site-specific coordination for the structural and electronic properties of a Co-HATCN network on silver | Ag(111) 위 Co 배위 HATCN 2D 네트워크(표면과학). HATCN 막이나 Ag 성장과 무관(arXiv 초록) |
| B-3 | 10.1109/access.2026.3717205 | 2026 | Effects of Charge Generation Layer Parameters on Luminance Balance Between Bottom and Top Emitting Units in Top Emitting Tandem Organic Light Emitting Diodes | HATCN/TAPC 대 HATCN/NPB CGL 최적화(주제 C). 형태 없음 |
| B-4 | 10.1038/s41598-026-35105-4 | 2026 | Micro-cavity-induced optical resonance for performance enhancement in ultra-thin CdTe photovoltaic devices | 로컬 전문 확인. Ag/HATCN/Ag 미세공진을 선행기술로만 언급. 데이터 없음 |
| B-5 | 10.1016/j.nxnano.2025.100132 | 2025 | OLED-on-silicon (OLEDoS) microdisplays: Technology challenges, design considerations, and adaptation in eXtended Reality (XR) ecosystem – Review | [Park&Suh 2018 피인용] OLEDoS 리뷰. 앵커를 일반 사례로 인용. HATCN/Ag 1차 데이터 없음(초록 수준) |
| B-6 | 10.1021/acsami.5c11689 | 2025 | Thermal Cross-linked Electron Transport Polymers for Suppressing Efficiency Roll-off in Green Solution-Processed Inverted OLEDs | [Park&Suh 2018 피인용] 용액공정 역구조 OLED용 ETL 고분자. HATCN 형태, Ag 성장 데이터 없음(초록) |
| B-7 | 10.1109/access.2025.3564910 | 2025 | Top Emission Organic Light Emitting Diodes Fabricated on Stainless Steel Foil | 스틸 기판 위 TEOLED. 거칠기는 평탄화층 값(8.61 → 0.46 nm)이고 HATCN은 HIL로만 씀. Al/Ag 6/8 nm 음극. 주제 F |
| B-8 | 10.1016/j.sse.2025.109066 | 2025 | Towards efficient and stable deep blue-light-emitting devices by employing HAT-CN doped TFB-x as the hole injection material | 용액 HAT-CN:TFB HIL(에너지 준위, 이동도). 확보한 정보에 형태 없음(주제 D) |
| B-9 | 10.1016/j.omx.2025.100419 | 2025 | Comprehensive optical characterization of organic thin films for OLED applications via spectroscopic ellipsometry | OLED 막 19종의 광학 상수(Yamaguchi 2020 피인용). 초록에서 HATCN 데이터가 보이지 않음. 관련성 낮음(단, HATCN nk 확인 후보 — 3절) |
| B-10 | 10.1016/j.mseb.2024.117271 | 2024 | Conjugated polymeric nanocomposite-based light-generating active materials for OLED applications: A review | [Park&Suh 2018 피인용] 고분자 나노복합 발광체 리뷰. HATCN, 박막 금속 성장과 무관 |
| B-11 | 10.1016/j.mseb.2023.117118 | 2024 | Solution processed PANI-rGO-SiO2 nanocomposite with inflated electron injection and transport capabilities as ETL for OLED applications | [Park&Suh 2018 피인용] 용액공정 나노복합 ETL. HATCN, 박막 금속 성장과 무관 |
| B-12 | 10.1021/acsami.3c18697 | 2024 | Robust Dipolar Layers between Organic Semiconductors and Silver for Energy-Level Alignment | Ag 위 BTB 단분자층. HATCN은 위층 시험용(웹 요약상 BTB/Ag 위에서 섬). HATCN 위 Ag 성장이 아님. 관련성 낮음 |
| B-13 | 10.1117/12.3025248 | 2024 | Investigation of transparent inverted OLEDs with bilateral emission using Li-doped BCP as electron injection layer | MoO3/Ag(16 nm)/HATCN 상부 전극(HATCN은 캡핑층), T 79.8%. Ag가 HATCN 위에 있지 않음. 주제 F 후보 |
| B-14 | 10.20944/preprints202410.0483.v1 | 2024 | A Study of the Electro-Optical Properties of Silicon-Based Top-Emission OLED | Si/Al/HATCN 10 nm TEOLED 프리프린트. 거칠기는 양극 이야기이고 HATCN 형태 데이터 없음 |
| B-15 | 10.1021/acsaem.4c02180 | 2024 | HAT-CN-Modified Spiro-OMeTAD Realizes Efficient Perovskite Solar Cells | 페로브스카이트 HTL 개질. 초록상 HAT-CN 버퍼가 spiro-OMeTAD 거칠기를 줄임(HAT-CN이 막을 거칠게 한다는 데 대한 약한 반대 단서). 박막 Ag 데이터 없음 |
| B-16 | 10.1016/j.orgel.2024.107060 | 2024 | Recent progress in semitransparent organic solar cells | 리뷰. Ag/HATCN/Ag 미세공진(Liang 2023)을 언급만 함 |
| B-17 | 10.1021/acs.jpcc.3c08262 | 2024 | Kinetic Trapping of Charge-Transfer Molecules at Metal Interfaces | Cu(111) 위 TCNE(HATCN 아님). 누운 배향에서 선 배향으로 바뀌는 시간이 주로 증착 속도로 정해진다는 일반론(초록). 배향 동역학의 배경 자료일 뿐 |
| B-18 | 10.1002/adfm.202214530 | 2023 | Highly Efficient Top‐Emitting Infrared‐to‐Visible Up‐Conversion Device Enabled by Microcavity Effect | [Park&Suh 2018 피인용] Mg:Ag 상부 음극을 쓴 전면발광 상향변환 소자. HATCN seed, 형태 데이터 없음(초록) |
| B-19 | 10.1002/adma.202300360 | 2023 | Trans‐Reflective Structural Color Filters Assisting Multifunctional‐Integrated Semitransparent Photovoltaic Window | Ag 20/HATCN/Ag 30 MIM 컬러필터(두꺼운 Ag). 초록에 형태 정보 없음. Sung 2021 행과 유사 |
| B-20 | 10.1002/pssr.202300001 | 2023 | Chemical Defects and Energetic Disorder Impact the Energy‐Level Alignment of Functionalized Hexaazatriphenylene Thin Films | HATCN 막의 화학 결함 갭 상태(주제 E). 형태 없음 |
| B-21 | 10.3390/coatings12101583 | 2022 | Transparent Conductive Indium Zinc Oxide Films: Temperature and Oxygen Dependences of the Electrical and Optical Properties | [Park&Suh 2018 피인용] 스퍼터 IZO 투명전도막. HATCN, 박막 Ag와 무관 |
| B-22 | 10.1016/j.vacuum.2022.111481 | 2022 | Nitrogen-modified a-C:H layers applied for light emitting diodes | [Park&Suh 2018 피인용] a-C:H 층. HATCN, 박막 Ag와 무관 |
| B-23 | 10.1039/d2na00862a | 2022 | High-efficiency and high-resolution patterned quantum dot light emitting diodes by electrohydrodynamic printing | 로컬 전문 확인. TE-QLED의 HAT-CN/Ag:Mg(Mg 20%) 10 nm 반투명 상부 음극, 증착 0.5–1 Å/s. 합금이고 형태, Rs 데이터 없음 |
| B-24 | 10.1109/jphot.2022.3159278 | 2022 | Full-Color Quantum Dot Light-Emitting Diodes Based on Microcavities | 로컬 전문 확인. MoO3 8/HATCN 15/Ag 30 nm 상부 전극(두꺼운 Ag). HATCN/Ag 형태 없음 |
| B-25 | 10.1063/5.0087868 | 2022 | Highly temperature-tolerant p-type carbon nanotube transistor doped with 1,4,5,8,9,11-hexaazatriphenylenehexacarbonitrile | CNT의 p-도펀트로 쓴 HATCN(200 °C 안정). 막 형태가 아님 |
| B-26 | 10.1007/s11664-021-09338-0 | 2021 | Organic Devices: Fabrication, Applications, and Challenges | [Park&Suh 2018 피인용] 유기 소자 전반 리뷰. HATCN, 박막 Ag 성장의 1차 데이터 없음 |
| B-27 | 10.1016/j.mtener.2021.100704 | 2021 | Fabrication of an oxide/metal/oxide structured electrode integrated with antireflective film to enhance performance in flexible organic light-emitting diodes | [Park&Suh 2018 피인용] 스퍼터 ZnO/AgOx/ZnO OMO(원자 산소 매개 Ag). HATCN 없음. 앵커는 일반 인용(로컬 전문 확인). 주제 A(산소 매개 Ag 성장)와는 관련 가능 |
| B-28 | 10.1007/s11664-021-08776-0 | 2021 | Realizing 8 cd A−1 Current Efficiency for Solution-Processed Inverted Top-Emitting Polymer Light-Emitting Diodes | [Park&Suh 2018 피인용] 용액공정 역구조 전면발광 PLED. 초록이 없고 HATCN 형태나 Ag 성장 데이터의 단서도 없음. 유료 |
| B-29 | 10.1016/j.synthmet.2021.116764 | 2021 | Ammonia gas sensors using 1,4,5,8,9,11-hexaazatriphenylene hexacarbonitrile semiconductor films | HATCN 가스센서 막. 초록 없음, 접근 불가. 관련성 낮음 |
| B-30 | 10.1016/j.rser.2020.110043 | 2020 | Organic light emitting diode devices: An energy efficient solid state lighting for applications | [Park&Suh 2018 피인용] OLED 일반 리뷰. HATCN 형태나 금속 성장 데이터 없음(S2/OpenAlex/COCI 피인용) |
| B-31 | 10.1002/adom.201901421 | 2020 | Microcavity‐Enhanced Blue Organic Light‐Emitting Diode for High‐Quality Monochromatic Light Source with Nonquarterwave Structural Design | [Park&Suh 2018 피인용] DBR 미세공진 청색 OLED. 앵커는 미세공진 맥락으로 인용. HATCN/Ag 성장 데이터 없음(초록) |
| B-32 | 10.1016/j.orgel.2019.105611 | 2020 | Determination of emitting dipole orientation in organic light emitting diodes | [Park&Suh 2018 피인용] 발광 쌍극자 배향 측정법. 무관 |
| B-33 | 10.1021/acsaelm.0c00823 | 2020 | Substrate-Independent Energy-Level Pinning of an Organic Semiconductor Providing Versatile Hole-Injection Electrodes | 수 nm HAT-CN 막이 갭 상태를 통해 모든 기판에서 WF를 약 5.60 eV로 고정(주제 C/E). 형태 없음 |
| B-34 | 10.1109/eltech.2019.8839348 | 2019 | Light Outcoupling Integrated Photonic Crystals for OLEDs | [Park&Suh 2018 피인용] 광결정 광추출 학회 논문. 무관(OpenAlex/S2 피인용, COCI에는 없음) |
| B-35 | 10.5757/asct.2019.28.6.217 | 2019 | Effect of Deposition Sequence in MoO3 and HAT-CN Dual-anode Buffer Layer on Inverted Organic Photovoltaics | 전문 확인. HAT-CN 5 nm/MoO3/Al 100 nm 양극 버퍼. 전기 특성만 있고 형태 없음. 증착 속도 0.01 nm/s |
| B-36 | 10.1016/j.orgel.2018.08.051 | 2018 | Tailoring Mg:Ag functionalities for organic light-emitting diodes | [Park&Suh 2018 피인용] Mg:Ag 음극 연구. 로컬 전문에 HATCN 언급 0회. 앵커는 Ag Rs·두께 맥락으로 인용. 주제 A/F 후보(Mg:Ag) |
| B-37 | 10.1016/j.orgel.2018.09.009 | 2018 | Thermally stable benzo[f]quinoline based bipolar host materials for green phosphorescent OLEDs | [Park&Suh 2018 피인용] 호스트 재료 합성(Suh 그룹). 앵커는 소자 구조 맥락으로 인용. HATCN/Ag 형태 없음 |
| B-38 | 10.1016/b978-0-12-409547-2.12894-x | 2018 | Kinetics of Ultra-Thin Organic Film Growth | Winkler 백과사전 챕터. Winkler Surf. Sci. 2016 리뷰와 같은 내용(Frank 2010 행에서 다룸) |
| B-39 | 10.1002/sdtp.12831 | 2018 | P‐13.4: Enhancement of top emission white organic light‐emitting diodes with a buffer layer | TCO 음극 아래 HATCN 23 nm 스퍼터 버퍼와 HATCN/Mg:Ag 하이브리드 버퍼. 형태, Rs 데이터 없음(초록) |
| B-40 | 10.1016/j.orgel.2018.02.029 | 2018 | The application of the nanostructure aluminum in the blue organic light-emitting devices | 초록 없음. HAT-CN/Al 키워드만 일치. 평가 불가 |
| B-41 | 10.1088/1361-6528/aa6567 | 2017 | Tuning the electrical property of a single layer graphene nanoribbon by adsorption of planar molecular nanoparticles | 그래핀 나노리본 위에 용액으로 캐스팅한 HAT-CN 나노입자. 증착막과 무관 |
| B-42 | 10.1016/j.cap.2017.04.012 | 2017 | In-situ analysis energy level alignment at solution processed HAT(CN)6/PVK (PVK:TAPC) interface via XPS and UPS | PVK 위 용액공정 HAT(CN)6의 에너지 준위(주제 E). 초록 없음 |
| B-43 | 10.1016/j.cplett.2016.12.007 | 2017 | Energy level alignment at the interface of NPB/HAT-CN/graphene for flexible organic light-emitting diodes | 에너지 준위 정렬(주제 E). 초록 없음, 형태 없음 |
| B-44 | 10.1016/j.vacuum.2017.07.007 | 2017 | Study of injection and transport properties of metal/organic interface using HAT-CN molecules as hole injection layer | 금속/HAT-CN 주입 연구. 초록 없음. PDF needed |
| B-45 | 10.1016/j.orgel.2016.12.021 | 2017 | Highly efficient ITO-free organic light-emitting diodes employing a roughened ultra-thin silver electrode | 로컬 전문 확인. 거칠게 만든 Ag 12 nm 양극 위 HAT-CN 10 nm(Ag 위 HATCN, 하부 전극). HATCN 형태 없음 |
| B-46 | WO2016171406A2 | 2016 | 유기 발광 소자 (Organic light-emitting device), LG Chem | Google Patents 본문 확인. HAT-CN은 실시예 1에서 ITO 위 HIL로만 쓰였고 Ag는 일반 음극 목록에만 있음. HAT-CN/Ag seed 내용 없음 |
| B-47 | 10.1021/am5091066 | 2015 | High-Performance Hybrid Buffer Layer Using 1,4,5,8,9,11-Hexaazatriphenylenehexacarbonitrile/Molybdenum Oxide in Inverted Top-Emitting Organic Light-Emitting Diodes | HATCN/MoO3 스퍼터 버퍼. O2 플라스마 후 버퍼의 AFM(초록에 수치 없음). TCO 상부 전극이라 금속 성장이 아님 |
| B-48 | (없음) | 2014 | ORGANIC/METAL INTERFACES FOR HIGH PERFORMANCE ORGANIC DIODES (SNU dissertation, 강찬모, hdl 10371/118972) | 한국 박사학위논문(초록 확인). 펜타센 다이오드의 전기적 어닐링. HAT-CN은 비교용 HIL로만 등장. HAT-CN 형태 없음 |
| B-49 | 10.1016/j.orgel.2014.09.031 | 2014 | Efficiency enhancement in tandem organic light-emitting devices with a hybrid charge generation layer composed of BEDT-TTF-doped TPBi/mCP/HAT-CN | 하이브리드 CGL(BEDT-TTF:TPBi/mCP/HAT-CN). HATCN:TPBi 혼합이 아님. 전기 특성(주제 C/D) |
| B-50 | 10.1039/c4nr01768d | 2014 | Highly efficient and low voltage silver nanowire-based OLEDs employing a n-type hole injection layer | 두꺼운(>150 nm) HAT-CN HIL을 쓴 AgNW 전극(웹 요약). RSC 403. 증착 박막 Ag가 아님. 관련성 낮음 |
| B-51 | 10.1063/1.4865765 | 2014 | Highly efficient inverted top emitting organic light emitting diodes using a transparent top electrode with color stability on viewing angle | HATCN/IZO 스퍼터 상부 전극. 금속 성장, 형태 없음(초록) |
| B-52 | 10.1063/1.4869183 | 2014 | Stability improvement of organic light emitting diodes by the insertion of hole injection materials on the indium tin oxide substrate | HAT-CN6 HIL이 ITO의 In/O 외부 확산을 막음(XPS). 안정성 연구이고 형태 아님(주제 C) |
| B-53 | 10.1063/1.4902871 | 2014 | A high performance semitransparent organic photodetector with green color selectivity | Kim 2013을 인용한 반투명 OPD. 초록에 HATCN 형태, Ag 성장 데이터 없음 |
| B-54 | 10.1021/jp407173w | 2013 | Linear and Hexagonal Porous Structures of an Organic Charge Acceptor Hexaaza-triphenylene-hexacarbonitrile on Au(111) with CN···CN Dipolar Interactions | 80 K, Au(111) 위 단분자층 이하 STM. 소자 두께에서의 막 결정화와 무관(초록 검토) |
| B-55 | 10.1016/j.orgel.2013.02.024 | 2013 | High-performance inverted top-emitting green electrophosphorescent organic light-emitting diodes with a modified top Ag anode | 앵커 역방향 참고문헌 [27](HAT-CN/Ag 상부 양극, Kippelen 그룹). Crossref/OpenAlex/S2에 초록 없고 Elsevier 본문 접근 불가. PDF needed |
| B-56 | 10.1063/1.4829135 | 2013 | Stacked inverted top-emitting white organic light-emitting diodes composed of orange and blue light-emitting units | HAT-CN/Al/LiF 연결 단위(Kippelen 그룹). 형태 없음(초록) |
| B-57 | 10.1016/j.synthmet.2011.12.028 | 2012 | Low voltage top-emitting organic light emitting devices by using 1,4,5,8,9,11-hexaazatriphenylene-hexacarbonitrile | 앵커 역방향 참고문헌 [26](Suh 그룹 HATCN TEOLED). 초록 없음, Elsevier 접근 불가. PDF needed |
| B-58 | 10.1016/j.orgel.2012.02.008 | 2012 | High current conduction with high mobility by non-radiative charge recombination interfaces in organic semiconductor devices | 앵커 역방향 참고문헌 [25](p-n 접합 계면). 전기 특성이고 형태 아님 |
| B-59 | 10.1002/adfm.201200185 | 2012 | Origin of Enhanced Hole Injection in Inverted Organic Devices with Electron Accepting Interlayer | 상부 양극 아래 NPB 위 HAT-CN/MoO3 중간층의 주입 효율(주제 C/D). 형태 없음 |
| B-60 | 10.1016/j.synthmet.2012.07.021 | 2012 | Electron mobility in hexaazatriphenylene hexacarbonitrile field-effect transistors | HATCN OFET(μ 2.5e-4 cm²/Vs, 검색 스니펫상 'domain-like nanoscale' 형태). 초록 확보 불가. PDF needed |
| B-61 | 10.1063/1.3662878 | 2011 | Interlayer molecular diffusion and thermodynamic equilibrium in organic heterostructures on a metal electrode | Koch 그룹. Bröker 2010 행에 통합(HATCN이 α-NPD 단분자층을 뚫고 Ag로 확산) |
| B-62 | (없음) | 2011 | Dependence of dark spot growth on cathode/organic interfacial adhesion in Organic Light Emitting Devices (U. Waterloo MASc thesis, R. Phatak) | 학위논문 전문 확인. 음극/유기 접착 연구, HATCN 없음 |
| B-63 | 10.18452/16257 | 2011 | Electronic and structural properties of interfaces between electron donor & acceptor molecules and conductive electrodes (Humboldt-Universität zu Berlin dissertation; DataCite DOI, not Crossref) | 박사학위논문(OA, edoc.hu-berlin.de/18452/16909). 여러 전극 위 HATCN, 밀도에 따른 재배향. Bröker PRL 2010 행과 같은 내용이라 통합(초록 확인) |
| B-64 | 10.1016/j.cplett.2009.04.019 | 2009 | Layer growth and desorption kinetics of a discoid molecular acceptor on Au(1 1 1) | Au(111) 위 HATCN(Frank 외, Winkler 그룹). Frank 2010 행에 통합(결과는 Winkler 2016 리뷰로 요약) |
| B-65 | 10.1016/j.tsf.2009.05.011 | 2009 | Efficient hole injection by doping of hexaazatriphenylene hexacarbonitrile in hole transport layer | HAT-CN 도핑 HTL(H4 관련)이나 출판사가 초록을 생략. PDF needed |
| B-66 | 10.1007/s10876-004-5772-6 | 2004 | Experimental and Computational Studies of Charge-Transfer and Reduction Products of 1, 4, 5, 8, 9, 11-Hexaazatriphenylene-Hexacarbonitrile: HAT-(CN)6 | 벌크 결정 구조와 CT 염. 배경 자료일 뿐(Winkler 리뷰가 HATCN T-stacking 패킹 근거로 인용). 접근 불가 |

### 4.3 C·D. 전하생성층·혼합 억셉터 (28편)
| # | DOI | 연도 | 제목 | 제외 사유 |
|---|---|---|---|---|
| CD-1 | 10.1016/j.jiec.2025.06.044 | 2026 | Achieving highly stable tandem organic light-emitting diodes using interfacial property control of charge generation layer | PDF needed: CGL 계면 특성 제어로 탠덤 안정화(2025–26, n-CGL/HAT-CN 계면 가능성) — 초록 미제공 |
| CD-2 | 10.1002/sdtp.70003_519 | 2026 | P‐14.11: Charge Generation Layer Multi‐stacking Strategy for High‐Efficiency Tandem Organic Light Emitting Diodes | 제목만(CGL 다중 적층, SID 2026) — 수치 확인 불가 |
| CD-3 | 10.1038/s41566-025-01679-0 | 2025 | Stable, deep blue tandem phosphorescent organic light-emitting diode enabled by the double-sided polariton-enhanced Purcell effect | 주제 F: 탠덤 CGL NPD/HATCN 12 nm/BPyTP2:Li 5 vol%, 상부 Ag 양극 아래 HATCN 5 nm — CGL 전기 수치 없음 |
| CD-4 | 10.1016/j.mssp.2023.108045 | 2024 | Effects of charge generation layers on multiple guest/host bilayer-based tandem OLEDs | PDF needed: Liq/HAT-CN vs Liq/MoO3 CGL 비교(검색 요약상 Liq/HAT-CN 주입능 우수) — 초록 직접 확인 불가 |
| CD-5 | 10.1016/j.displa.2024.102713 | 2024 | Dual side transparent organic light-emitting diodes with a modified Ag top cathode | HATCN 없음(LiF/Ag 음극 TrOLED) — 주제 F |
| CD-6 | 10.1002/sdtp.17565 | 2024 | 37‐4: Degradation Analysis on Li‐doped Organic Charge Generation Layer | 같은 그룹 중복: Kim et al. J. Inf. Disp. 2024(Li 도핑 n-CGL 열화)와 동일 주제 SID 발표 |
| CD-7 | 10.1016/j.orgel.2021.106353 | 2022 | Lateral current suppression in tandem organic light-emitting diodes by adopting a buffer layer | PDF needed: 탠덤 CGL 측면전류 억제 버퍼(HATCN 측면 전도 관련 가능) — 초록 미제공 |
| CD-8 | 10.1016/j.optmat.2022.112345 | 2022 | Interfacial charge transfer study of hexacarbonitrile-based intermediate connecter in blue tandem organic light-emitting diodes | PDF needed: 헥사카보니트릴 연결층 계면 전하이동(청색 탠덤) — 초록 미제공 |
| CD-9 | 10.1016/j.displa.2022.102306 | 2022 | Highly efficient tandem OLED based on a novel charge generation layer of HAT-CN/CuPc heterojunction | HAT-CN/CuPc 이종접합 CGL(p측) — n-ETL/HATCN 질문과 관련 낮음 |
| CD-10 | 10.1016/j.orgel.2021.106202 | 2021 | Efficient inverted organic light-emitting devices using a charge-generation unit as electron-injection layers | PDF needed: 역구조 OLED에서 CGU를 전자주입층으로 사용(팀 Ag/HATCN/n-ETL과 같은 물리) — 초록 미제공 |
| CD-11 | 10.1016/j.orgel.2020.105979 | 2021 | Investigation on the mechanism of charge generation in organic heterojunctions: Analysis of I–V and C–V characteristics | PDF needed: 유기 이종접합 전하생성 메커니즘 I–V/C–V — 초록 미제공 |
| CD-12 | 10.1002/adpr.202000122 | 2021 | A Deep Blue Strong Microcavity Organic Light‐Emitting Diode Optimized by a Low Absorption Semitransparent Cathode and a Narrow Bandwidth Emitter | Ag 단독 반투명 음극 강공진 TEOLED(Suh 그룹) — HATCN/n-ETL 여부 미확인, 주제 F |
| CD-13 | 10.1088/1361-6641/ab39f7 | 2019 | High efficiency organic light-emitting diodes based on HAT-CN/TAPC heterojunction charge generation layer as charge injectors | 같은 그룹 중복(Ma D.): Yuan 2020 Org. Electron. 행 note에 요약 — HAT-CN/TAPC OHJ 주입체, HAT-CN/ETL 사이 초박막 Al |
| CD-14 | 10.1103/physrevapplied.12.064052 | 2019 | Dipolar Doping of Organic Semiconductors to Enhance Carrier Injection | 개념적 관련만: 극성 분자(ETL)를 비극성 매트릭스에 희석하는 dipolar doping — LUMO 희석(H4)과 다른 메커니즘 |
| CD-15 | 10.1016/j.orgel.2017.10.042 | 2018 | Highly efficient tandem organic light-emitting devices adopting a nondoped charge-generation unit and ultrathin emitting layers | PDF needed: 비도핑 CGU 탠덤(LiF/Al/HAT-CN 추정) — 초록 미제공 |
| CD-16 | 10.1016/j.orgel.2017.11.022 | 2018 | Highly-efficient tandem organic light-emitting device employing bis-4,6-(3,5-di-3-pyridylphenyl)-2-methylpyrimi-dine doped with cesium azide in charge generation unit | PDF needed: B3PYMPM:CsN3 n-CGL — 초록 미제공 |
| CD-17 | 10.1038/srep43130 | 2017 | Stacking multiple connecting functional materials in tandem organic light-emitting diodes | HATCN 미사용(Liq/Al/MoO3 + C60 버퍼) — 연결층 측면전류 사례로만 참고 |
| CD-18 | 10.1016/j.orgel.2015.05.033 | 2015 | The role of gap states on energy level alignment at an α-NPD/HAT(CN) 6 charge generation interface | PDF needed: alpha-NPD/HAT(CN)6 전하생성 계면 갭 상태(p측) — 초록 미제공 |
| CD-19 | 10.1016/j.orgel.2013.11.033 | 2014 | Highly efficient p-dopants in amorphous hosts | PDF needed: 비정질 호스트 내 p-도펀트 효율 — 초록 미제공 |
| CD-20 | 10.1016/j.orgel.2014.10.015 | 2014 | Role of hole injection layer in intermediate connector of tandem organic light-emitting devices | PDF needed: 연결층 내 HIL(HATCN 등) 역할 — 초록 미제공 |
| CD-21 | 10.1039/c4tc02082k | 2014 | A novel intermediate connector with improved charge generation and separation for large-area tandem white organic lighting devices | PDF needed: 전하생성·분리 개선 연결층(Ding 2014) — 초록 한 문장뿐, 본문 RSC 403 |
| CD-22 | 10.1016/j.tsf.2014.09.012 | 2014 | Efficiency enhancement of tandem organic light-emitting devices by a combined charge generation layer and organic n-type bis(ethylenedithio)-tetrathiafulvalene-doped electron transport layer | 같은 그룹(Kim T. W.?) BEDT-TTF 도핑 ETL CGL — 초록 미제공 |
| CD-23 | 10.1016/j.materresbull.2014.04.045 | 2014 | Enhancement of efficiencies for tandem green phosphorescent organic light-emitting devices with a p-type charge generation layer | p형 CGL 탠덤 — 초록 미제공·관련성 낮음 |
| CD-24 | 10.1016/j.jlumin.2014.04.027 | 2014 | Efficient green phosphorescent tandem organic light emitting diodes with solution processable mixed hosts charge generating layer | 용액공정 혼합 호스트 CGL — 초록 미제공·HATCN 관련 불명 |
| CD-25 | 10.1016/j.orgel.2012.06.039 | 2012 | Determination of the interface energy level alignment of a doped organic hetero-junction using capacitance–voltage measurements | PDF needed(최우선): Kim J.-J. 그룹 C–V로 도핑 이종접합(n-ETL/HATCN 추정) 계면 정렬 결정 — 초록·본문 미확보로 행 제외 |
| CD-26 | 10.1117/12.929559 | 2012 | The mechanism of charge generation in charge generation units containing HATCN for high-luminance tandem OLED display | 중복: Lee et al. AFM 2011(10.1002/adfm.201102212)과 같은 내용의 SPIE 발표 |
| CD-27 | 10.1002/adfm.201290023 | 2012 | Organic Light‐Emitting Diodes: The Mechanism of Charge Generation in Charge‐Generation Units Composed of p‐Doped Hole‐Transporting Layer/HATCN/n‐Doped Electron‐Transporting Layers (Adv. Funct. Mater. 4/2012) | 중복: Lee et al. AFM 2011의 표지(frontispiece) 기사 |
| CD-28 | 10.1016/j.orgel.2009.05.007 | 2009 | P-type doping of organic wide band gap materials by transition metal oxides: A case-study on Molybdenum trioxide | PDF needed: MoO3를 넓은 갭 호스트에 도핑(희석)한 고전 연구 — 초록 미제공 |

### 4.4 E. 계면 전하 이동 (23편)
| # | DOI | 연도 | 제목 | 제외 사유 |
|---|---|---|---|---|
| E-1 | 10.1039/d6el00029k | 2026 | Energetics and chemistry at electron selective interfaces for p-i-n perovskite solar cells: an in situ investigation | 초록이 일반적(수치 없음), RSC 전문 자동 접근 불가. 2차 스니펫상 BCP 위 Ag의 착물(N 1s 신규 성분) 보고 → PDF needed |
| E-2 | 10.1021/acs.jpcc.5c04272 | 2025 | Dial It Down: The Effect of Strongly Interacting Adsorbates on the BiAg2 Rashba Surface State | BiAg2 Rashba 표면상태; seed 문제와 무관 |
| E-3 | arXiv:2510.11470 | 2025 | Enhancing the Plasmonic Hotspot Density via Structural Engineering of Multi-layered MoO3-Ag-Au Systems Under Extreme Electronic Excitation Conditions for Ultra-Sensitive SERS Applications | 이온 조사 SERS 다층막; seed 계면 전하와 무관(스니펫만 확인) |
| E-4 | 10.1021/acs.jpcc.4c08104 | 2025 | Integer Charge Transfer Model–PTCDA on MgO(001)/Ag(001) Probing the Transition from Single to Double Integer Charge Transfer | 유기 분자의 산화물/금속 위 전하 이동; Ag 성장과 무관 (로컬 전문 존재) |
| E-5 | 10.1021/acs.jpcc.5c05641 | 2025 | Microscopic Origin of Charge Transfer at the Organic Semiconductor/MoO3 Hybrid Interface | 유기/MoO3 계면; Ag 없음 (주제 C/D 참고 가능) |
| E-6 | 10.1021/acs.jpcc.2c01432 | 2022 | Direct Experimental Evidence for Substrate Adatom Incorporation into a Molecular Overlayer | Au(111)-F4TCNQ(Ag 아님). Ag 대응 결과는 Ryan 2022(Ag(100)), Mousley 2023·Blowey 2018(Ag(111))로 대체 |
| E-7 | 10.3390/nano11051292 | 2021 | Charge-Transfer Induced by the Oxygen Vacancy Defects in the Ag/MoO3 Composite System | SERS용 공스퍼터 Ag/MoO3 복합체; 본문에 Ag 3d/Mo 3d 정량값 없음(그림만), seed/박막 연속성과 무관 (Europe PMC 전문 확인) |
| E-8 | 10.1039/c9nj03769a | 2020 | The adsorption and growth of Ag n (n = 1–4) clusters on cubic, monoclinic, and tetragonal ZrO2 surfaces: a first-principles study | 초록이 한 문장뿐(수치 없음), 원문 미열람; 스니펫상 '큰 전하 이동 없음'만 확인 |
| E-9 | 10.1039/c9ta08438j | 2019 | Substoichiometric ultrathin zirconia films cause strong metal–support interaction | ZrO2/Pt 역촉매 SMSI; Ag 없음 |
| E-10 | 10.1021/acs.jpcc.8b08478 | 2019 | Core–Shell Charge Transfer in Plasmonic Fe@Ag Nanoparticles on MgO Film | Fe@Ag 코어–셸 입자 내부 전하 이동; seed–Ag 계면 아님 |
| E-11 | 10.1016/j.ceramint.2016.09.157 | 2017 | Tailoring transparence in MoOx/Ag/MoOx electrode through Ag by O2/Ar plasma exposure | Ag 플라즈마 처리 MoOx/Ag/MoOx(스퍼터) 전극 Rs/T 연구 → 주제 A(센서스 P052)로 충분, 계면 전하 분석 없음 |
| E-12 | 10.1021/acsnano.7b02449 | 2017 | Charge Transfer and Orbital Level Alignment at Inorganic/Organic Interfaces: The Role of Dielectric Interlayers | Ag(001) vs MgO/Ag(001) 위 펜타센의 분수/정수 전하 이동; 유기 흡착체 쪽 연구로 Ag 성장과 무관 |
| E-13 | 10.1021/acs.jpcc.5b03768 | 2015 | Fingerprint of Fractional Charge Transfer at the Metal/Organic Interface | Au 위 tetrafluoropentacene 물리흡착; 억셉터/Ag 계면 아님 (arXiv 전문 확인) |
| E-14 | 10.1016/j.vacuum.2012.05.031 | 2012 | X-ray photoelectron spectroscopy: Silver salts revisited | Ag 화합물 XPS 기준값 문헌이나 초록 미제공·미열람; Ag 3d 음의 이동 근거는 Chen 2022·Kaspar 2010·Loi 2023으로 대체 |
| E-15 | 10.1038/nchem.591 | 2010 | Charge-transfer-induced structural rearrangements at both sides of organic/metal interfaces | Cu(100) 위 TCNQ/F4TCNQ(금속 쪽 원자 들림)의 고전 논문이나 Ag 아님, 초록 API 미제공 → 노트에서만 언급 |
| E-16 | 10.1063/1.3427430 | 2010 | Charge generation layers comprising transition metal-oxide/organic interfaces: Electronic structure and charge generation mechanism | TMO/유기 CGL(주제 C); Ag 없음 |
| E-17 | 10.1016/j.surfrep.2009.07.004 | 2009 | Properties of oxide thin films and their adsorption behavior studied by scanning tunneling microscopy and conductance spectroscopy | 총설(Nilius). OA 사본(MPG)은 bot 검사로 접근 불가, 미열람 |
| E-18 | 10.1016/j.surfrep.2007.07.001 | 2007 | Interaction of nanostructured metal overlayers with oxide surfaces | 핵심 총설(Fu & Wagner)이나 미열람(초록 elided) → PDF needed |
| E-19 | 10.1063/1.2118587 | 2005 | Pinning mass-selected Agn clusters on the TiO2(110)−1×1 surface via deposition at high kinetic energy | 제목 선별: 고에너지 착지 pinning 연구, 전하 상태·젖음 자료 아님 |
| E-20 | 10.1063/1.1859271 | 2005 | Landing of size-selected Agn+ clusters on single crystal TiO2 (110)-(1×1) surfaces at room temperature | 제목 선별: 크기 선택 Ag 클러스터 착지 형태; 전하 상태 정량 없음(초록 미확인) |
| E-21 | 10.1016/j.ssc.2004.04.007 | 2004 | The effect of heat treatment on bistable Ag-TCNQ thin films | 초록·전문 미확보(Elsevier, 초록 elided). 검색 스니펫에만 'N 1s 새 상태(Ag→TCNQ 전자 이동)' 언급 → 검증 불가. PDF needed |
| E-22 | 10.1016/j.synthmet.2004.04.012 | 2004 | The preparation of nanograin Ag–TCNQ thin films by vacuum evaporation with post-heat treatment | 초록·전문 미확보 → 내용 검증 불가 |
| E-23 | 10.1016/S0167-5729(96)00011-8 | 1997 | Ultrathin metal films and particles on oxide surfaces: structural, electronic and chemisorptive properties | 핵심 총설(Campbell)이나 초록 없음·paywall로 미열람 → 행에 넣지 않음. PDF needed(금속/산화물 성장 모드와 전하 이동·반응성 상관의 표준 인용처) |

### 4.5 F. 박막 Ag 음극 소자 (25편)
| # | DOI | 연도 | 제목 | 제외 사유 |
|---|---|---|---|---|
| F-1 | 10.1016/j.jiec.2026.05.027 | 2026 | The effects of metal interlayer structures on the performance of phosphorescent transparent organic light-emitting diodes | 변수는 하부 ITO/Ag interlayer; 상부 BPhen:Cs 60/Ag 15/NPD 고정, 상부 전극 형태·seed 정보 없음 (DB P092 검토) |
| F-2 | 10.1002/adma.202518212 | 2026 | High‐Performance Flexible Near‐Infrared Organic Light‐Emitting Diodes with a Heterogeneous Alloy Semitransparent Electrode | Mg-Bi 합금 반투명 전극(NIR OLED) — Ag 박막 아님 (abstract 검토) |
| F-3 | 10.1080/15421406.2025.2612133 | 2026 | Plasmon–microcavity coupling and spectral modulation in blue OLEDs employing Ag/ZnS/Ag multilayer electrodes | Ag/ZnS/Ag 하부 양극 플라즈몬-공진 결합 연구 — 상부 박막 Ag 아님 (abstract 검토) |
| F-4 | 10.1093/bulcsj/uoag003 | 2026 | A multifunctional dipyrimidylpyridine electron-transporter combined with Ag for robust and efficient organic light-emitting diodes | Ag 도핑 DPmPy-DBF1 ETL (Ag 배위) — 상부 박막 Ag 전극 자료 없음, 주제 C 참고 |
| F-5 | 10.1016/j.mtchem.2025.102612 | 2025 | High color gamut top-emitting organic light-emitting diode based on metal/transparent conductive oxide composite electrode | 변수는 하부 반사전극(Al/ITO vs Ag/ITO); 상부 Bphen:Ag2O/Ag 20 고정, 형태 정보 없음. 참고: EQEmax 20.74% (Table 2) |
| F-6 | 10.1002/adma.202507556 | 2025 | Increasing the Stability of Deep Blue Phosphor‐Sensitized OLEDs Using the Polariton‐Enhanced Purcell Effect | 스퍼터 Ti/Ag 16 nm 상부 양극의 polariton Purcell 효과로 수명 개선; seed 형태·열화 기구가 아닌 광학적 안정화 (DB P103 검토) |
| F-7 | 10.3390/coatings15121370 | 2025 | Recent Progress in Dielectric/Ag/Dielectric Transparent Electrodes on Flexible Substrates | 리뷰(유연 기판 D/Ag/D 전극) — 참고 리뷰로 notes에 기재 |
| F-8 | 10.1039/d4tc03534h | 2024 | A pyrimidine end-capped electron transport material interacted with silver for improving electron-injection and long-term stability in OLEDs | DPmPy-BP ETL + Ag EIL (Ag–N 상호작용으로 주입·LT50 17,000 h) — 박막 Ag 상부 전극이 아닌 EIL 용도, 주제 C/E 참고 |
| F-9 | 10.1515/ntrev-2023-0578 | 2023 | Recent progress in dielectric/metal/dielectric electrodes for foldable light-emitting devices | 리뷰(DMD 전극, 접이식 LED) — 1차 수치 없음, 참고 리뷰로 notes에 기재 |
| F-10 | 10.3390/nano13071282 | 2023 | Research on the High Light Out-Coupling Efficiency Deep-Blue Top-Emitting Organic Light-Emitting Diode through FDTD Optical Simulation | FDTD 광학 설계 중심 청색 TOLED — 상부 전극 seed·형태 자료 없음 (abstract 검토) |
| F-11 | 10.1063/5.0132938 | 2023 | Highly transparent and stretchable organic light-emitting diodes with ultrathin metal films as double electrodes | 초박막 금속 양 전극 신축 투명 OLED(T 73%) — 초록에 seed·두께·형태 정보 없음 |
| F-12 | 10.35848/1882-0786/abdac5 | 2021 | Optical outcoupling in monochromatic top-emitting organic light-emitting diodes with an Au/Ag semi-transparent electrode | Au 2/Ag 7 nm 상부 전극 TEOLED 수치 시뮬레이션만(실험 수치 없음) |
| F-13 | 10.1063/5.0057184 | 2021 | A highly transparent laminated composite cathode for organic light-emitting diodes | ALD ZnO/Mg:Ag/ZnO 합금 음극 — 순수 Ag 아님 |
| F-14 | 10.1016/j.orgel.2020.105900 | 2020 | Transparent bi-directional organic light-emitting diodes with color-tunable top emission | BPhen:Cs/Ag 15 표준 음극 위 색변환 캡핑층 연구; seed·형태·열화 정보 없음 |
| F-15 | 10.1038/s41528-020-00096-3 | 2020 | Multi-directionally wrinkle-able textile OLEDs for clothing-type displays | 섬유 기판 평탄화·기계적 신뢰성 연구; 상부 MoO3 10/Ag 20 전극 특성 자료 없음 |
| F-16 | 10.1021/acsaelm.0c00156 | 2020 | Ultrathin Ag Transparent Conducting Electrode Structure for Next-Generation Optoelectronic Applications | Kwon 그룹 하부 투명전극: 유기 HIL 사이 Ag 9 nm (Rs 17.6 Ω/sq, T 79.2%, EQE 23.8%) — 상부 전극 아님. HIL 종류가 초록에 없어 HAT-CN 여부 미확인 → 주제 A/B로 이관, PDF needed |
| F-17 | 10.1088/1361-6528/ab51bb | 2019 | Solution-processed, top-emitting, microcavity polymer light-emitting diodes for the pure red, green, blue and near white emission | 용액공정 PLED 위 MoO3 10/Ag 15 상부 양극; 형태·seed 비교 없음, MoO3/Ag 상부 양극 대표 문헌과 중복 |
| F-18 | 10.1039/c8ra00298c | 2018 | Extremely flat metal films implemented by surface roughness transfer for flexible electronics | 평탄 금속막(거칠기 전사) 기판 연구; 상부 LiF/Al 2/Ag 18/WO3 30 음극 고정, 상부 전극 형태 자료 없음 |
| F-19 | 10.1002/adom.201801262 | 2018 | Ultrathin MoO<sub>3</sub> Layers in Composite Metal Electrodes: Improved Optics Allow Highly Efficient Organic Light‐Emitting Diodes | 하부발광 Au/Ag 복합 금속 전극의 초박막 MoO3 seed 광학 효과 — 상부 전극 아님, 주제 A/E(H2) 참고 |
| F-20 | 10.1016/j.orgel.2016.10.007 | 2016 | Effects of transparent bottom electrode thickness on characteristics of transparent organic light-emitting devices | 하부 투명전극 두께 변수(상부 Al/Ag 12 고정); F 주제와 간접 |
| F-21 | 10.1016/j.orgel.2014.08.004 | 2014 | Highly efficient inverted top emitting organic light emitting diodes using a horizontally oriented green phosphorescent emitter | 역구조 TEOLED(수평배향 녹색 인광) — 초록·본문 접근 불가로 상부 전극(HAT-CN/Ag 여부) 미확인, PDF needed |
| F-22 | 10.1016/j.optlastec.2012.07.007 | 2013 | Improving electron injection and microcavity effect for constructing highly efficient inverted top-emitting organic light-emitting diode | CsOx vs LiF 하부 EIL 비교; 상부 MoOx 2/Ag 20 고정 (Zhang 그룹 중복) |
| F-23 | 10.1016/j.orgel.2013.06.004 | 2013 | Stacked inverted top-emitting green electrophosphorescent organic light-emitting diodes on glass and flexible glass substrates | Kippelen 그룹 적층 역구조 TEOLED — 초록·본문 접근 불가, 상부 Ag 양극 구조 미확인, PDF needed |
| F-24 | 10.1016/j.jlumin.2011.07.010 | 2012 | A very simple method of constructing efficient inverted top-emitting organic light-emitting diode based on Ag/Al bilayer reflective cathode | 변수는 하부 반사 음극(Al vs Ag/Al); 상부 MoOx 2/Ag 20 고정, 형태 정보 없음 (같은 그룹 Opt. Laser Technol. 2013 doi:10.1016/j.optlastec.2012.07.007도 동일 사유) |
| F-25 | 10.1063/1.4736573 | 2012 | Highly efficient inverted top-emitting green phosphorescent organic light-emitting diodes on glass and flexible substrates | Kippelen 그룹 역구조 TEOLED이나 상부 양극이 MoO3/Au 20 nm (Ag 아님; local PDF 확인) |

### 4.6 G. 계산 방법 (34편)
| # | DOI | 연도 | 제목 | 제외 사유 |
|---|---|---|---|---|
| G-1 | 10.1021/acs.jpclett.5c00517 | 2025 | Tailoring Charge-Transfer at Metal–Organic Interfaces Using Designer Shockley Surface States | HATCN/Ag 박막 CT 상태 조절 — Topic E 소관(계산 방법론 아님). |
| G-2 | 10.1016/j.apsusc.2024.161774 | 2025 | Comparison of the effect of gaseous and metallic additives on ultrathin Ag layer formation: An experimental and numerical investigation | Yun 그룹 연작 — Yun 2026 행에 통합. |
| G-3 | 10.1016/j.ensm.2023.103035 | 2023 | A systematic study on the metallophilicity of ordered five-atomic-layer MXenes using high-throughput automated workflow and machine learning | 금속 음극용 HT DFT+ML '금속친화도' 스크리닝(가장 가까운 타 분야 선례 중 하나)이나 초록·본문 미확보 → PDF needed. |
| G-4 | 10.3390/ma15020648 | 2022 | Aluminum Doping Effect on Surface Structure of Silver Ultrathin Films | Al 도핑 Ag 표면 확산 DFT/MD — 합금(Ag:Al) 주제로 Topic A/B_alloy 소관; seed 스크리닝 아님. |
| G-5 | 10.1038/s41524-022-00852-1 | 2022 | Data-driven models for ground and excited states for Single Atoms on Ceria | 단원자/세리아 데이터 기반 모델 — O'Connor 2018·Liu 2020 행과 중복(촉매). |
| G-6 | 10.1021/acscatal.0c04525 | 2021 | Open Catalyst 2020 (OC20) Dataset and Community Challenges | Tran & Ulissi 2018 행에 고처리량 흡착 DB 사례로 통합(분자 흡착물/금속 표면, 금속 원자/유기 기판 아님). |
| G-7 | 10.1002/wcms.1493 | 2021 | Extended tight-binding quantum chemistry methods | xTB 리뷰 — GFN2 행 노트로 통합. |
| G-8 | 10.1002/cphc.202001052 | 2021 | Efficient Computation of Geometries for Gold Complexes | Bursch 2019 행 노트로 통합(GFN2 기하, PBEh-3c 대비 ~1 kcal/mol — 초록). |
| G-9 | 10.1063/5.0021955 | 2020 | OrbNet: Deep learning for quantum chemistry using symmetry-adapted atomic-orbital features | OrbNet Denali 2021 행에 전신으로 통합. |
| G-10 | 10.1016/j.apsusc.2020.145515 | 2020 | Strategy for improving Ag wetting on oxides: Coalescence dynamics versus nucleation density | Yun 그룹 연작 — Yun 2026 행에 통합(핵 밀도 vs 합체 동역학 논점은 Topic A에서 상세 검토 권장). |
| G-11 | 10.1016/j.apsusc.2019.05.208 | 2019 | Dynamics of 3D-island growth on weakly-interacting substrates | Gervilla 2020 Sci Rep 행에 통합(R_c 17.3→6.0 Å 수치는 그 행에 기재). |
| G-12 | 10.1103/PhysRevMaterials.3.023804 | 2019 | Active learning of uniformly accurate interatomic potentials for materials simulation | Smith 2018 행에 능동학습 사례(DP-GEN)로 통합. |
| G-13 | 10.1109/IITC.2018.8430425 | 2018 | A First-Principles Density Functional Theory Based Framework for Barrier Material Screening | Cu 배선 배리어 DFT 스크리닝(관련성 높음)이나 OpenAlex 초록이 다른 논문(Ru 배선) 내용으로 불일치 → 내용 미검증. PDF needed. |
| G-14 | 10.1021/acs.jcim.8b00256 | 2018 | High-Throughput Screening Approach for the Optoelectronic Properties of Conjugated Polymers | Wilbraham 2020 행에 원형 방법으로 통합(xTB→(TD)DFT 선형 보정). |
| G-15 | 10.1039/C7CP04913G | 2017 | A look at the density functional theory zoo with the advanced GMTKN55 database for general main group thermochemistry, kinetics and noncovalent interactions | 주족(main-group) 전용 벤치마크로 금속–유기 결합 수치 없음. 방법론 배경 인용 후보(OrbNet Denali·g-xTB가 이 세트로 평가). |
| G-16 | 10.1063/1.4993252 | 2017 | Scaling of elongation transition thickness during thin-film growth on weakly interacting substrates | Lü 2014/Elofsson 2014 행에 동일 그룹 연작으로 통합(θ_Elong 실험 결정). |
| G-17 | 10.1021/acs.jctc.7b00118 | 2017 | A Robust and Accurate Tight-Binding Quantum Chemical Method for Structures, Vibrational Frequencies, and Noncovalent Interactions of Large Molecular Systems Parametrized for All spd-Block Elements (Z = 1–86) | GFN1-xTB — GFN2 행 노트로 통합. |
| G-18 | 10.1016/j.actamat.2017.07.024 | 2017 | Controlling Ag diffusion in ZnO by donor doping: A first principles study | ZnO 벌크 내 Ag 확산 DFT — Wang 2018 JAP(동일 그룹) 행과 중복, 계면 핵생성 아님. |
| G-19 | 10.1038/srep40907 | 2017 | Property database for single-element doping in ZnO obtained by automated first-principles calculations | ZnO 벌크 도핑 고처리량 DB — Ag 계면·핵생성과 무관. |
| G-20 | 10.1063/1.4916983 | 2015 | Coalescence-controlled and coalescence-free growth regimes during deposition of pulsed metal vapor fluxes on insulating surfaces | Lü 2014 APL 행에 동일 그룹 연작으로 통합(펄스 플럭스 kMC). |
| G-21 | 10.1063/1.4866175 | 2014 | Insight into the description of van der Waals forces for benzene adsorption on transition metal (111) surfaces | 벤젠/금속 vdW 벤치마크 — Liu 2015 PRL·Maurer 2016과 중복, 별도 열람 안 함. |
| G-22 | 10.3390/cryst3010079 | 2013 | Metals on Graphene: Interactions, Growth Morphology, and Thermal Stability | Liu 2012 PCCP 행에 동일 그룹 리뷰로 통합(MDPI 자동접근 차단으로 본문 미열람). |
| G-23 | 10.1103/PhysRevB.84.045115 | 2011 | Formation enthalpies by mixing GGA and GGA+U calculations | FERE(Stevanović 2012) 행의 '기준점 보정' 사례로 통합. |
| G-24 | 10.1016/j.apsusc.2009.11.044 | 2010 | Interaction of bathocuproine with metals (Ca, Mg, Al, Ag, and Au) studied by density functional theory | BCP–Ag DFT(팀 ETL E_b와 직접 비교 대상)이나 초록·본문 미확보 → PDF needed. 전신 Bhatt 2009 JJAP 행 참조. |
| G-25 | 10.1088/1367-2630/11/5/053010 | 2009 | Theoretical study of PTCDA adsorbed on the coinage metal surfaces, Ag(111), Au(111) and Cu(111) | PTCDA/동전금속 DFT — vdW-surf·MBD 벤치마크 행과 중복. |
| G-26 | 10.1063/1.3159651 | 2009 | A first principles study of the properties of Al:ZnO and its adhesion to Ag in an optical coating | Wang 2018 JAP 행 노트로 통합(Bristowe 그룹). |
| G-27 | 10.1021/ja802430q | 2008 | Selective Metal Deposition on Photoswitchable Molecular Surfaces | Tsujioka 2016 리뷰 행에 원저로 통합. |
| G-28 | 10.1063/1.2815646 | 2007 | A density functional study of the effect of hydrogen on the strength of an epitaxial Ag/ZnO interface | Wang 2018 JAP 행 노트로 통합(Bristowe 그룹). |
| G-29 | 10.1016/S0379-6779(99)00451-8 | 2000 | The organic–cathode interface in Alq3-based organic light-emitting devices: new insights from ab-initio molecular dynamics | Curioni & Andreoni 1999 행에 통합(초록 미확보). |
| G-30 | 10.1016/S0167-5729(97)00015-0 | 1998 | Microscopic view of epitaxial metal growth: nucleation and aggregation | 핵생성 이론·STM 고전 리뷰(Brune) — 초록 미확보. kMC/핵 밀도 이론 배경 인용 후보. |
| G-31 | 10.1007/BF01112848 | 1991 | Energy-adjusted ab initio pseudopotentials for the second and third row transition elements: Molecular test for M2 (M=Ag, Au) and MH (M=Ru, Os) | Ag2 분자로 ECP를 시험한 논문(팀 Ag2 1.86 eV 검증과 관련) — 초록 미확보, PDF needed. |
| G-32 | 10.1007/BF01114537 | 1990 | Energy-adjusted ab initio pseudopotentials for the second and third row transition elements | def2 계열 Ag ECP 원전 — 방법 인용용, 내용 미열람. |
| G-33 | 10.1088/0034-4885/47/4/002 | 1984 | Nucleation and growth of thin films | Venables 1973 행(10.1080/14786437308219242)에 동일 저자 리뷰로 통합. |
| G-34 | 10.1080/00268977000101561 | 1970 | The calculation of small molecular interactions by the differences of separate total energies. Some procedures with reduced errors | Counterpoise(CP) 보정 원전 — 방법 인용 필수이나 본문/초록 미열람(수치 없음)이라 표 제외. 방법론 절 인용. |

### 4.7 H. HAT 유도체 (44편)
| # | DOI | 연도 | 제목 | 제외 사유 |
|---|---|---|---|---|
| H-1 | 10.1039/d6cc00678g | 2026 | Benzotriquinoxaline-based discotic liquid crystals stabilizing room temperature columnar self-assembly with high charge carrier mobility | 벤조트리퀴녹살린(=HATNA) 원반 액정의 상온 컬럼상·이동도 — 초록에 수치 없음, LC 중심 |
| H-2 | 10.3390/nano13091560 | 2023 | Rational Design of Lewis Base Electron Transport Materials for Improved Interface Property in Inverted Perovskite Solar Cells: A Theoretical Investigation | Lewis 염기 ETM 계산 설계 — 비대칭 PDI로 HAT 아님 |
| H-3 | 10.1038/s41563-022-01402-2 | 2022 | Dynamic molecular switches with hysteretic negative differential conductance emulating synaptic behaviour | 분자 스위치 NDC — HAT 유도체 물성 언급 없음(초록) |
| H-4 | 10.1021/acsaem.2c02377 | 2022 | Effect of a One-Dimensional Columnar Structure on the Cathode Active Material Performance of Single-Component Hexaazatriphenylene Derivatives | HATCNOC2 리튬이온전지 양극 — Hinokimoto 2022 행의 하위 항목(같은 그룹)으로 note에 기재 |
| H-5 | 10.1038/s41586-019-1435-5 | 2019 | High performance from extraordinarily thick organic light-emitting diodes | 초두께 OLED(Matsushima 2019) — HAT-CN 수송층 관련이나 주제 C/D 범위, Yamaguchi 2020 행 note에 기재(미검토) |
| H-6 | 10.1016/j.tetlet.2018.01.004 | 2018 | Recent advances of hexaazatriphenylene (HAT) derivatives: Their applications in self-assembly and porous organic materials | HAT 유도체 리뷰(자기조립·다공성 재료) — 초록 미확보, 다공성/자기조립 중심이라 H5 비핵심 |
| H-7 | 10.1039/c8sc02810a | 2018 | Rational design of red AIEgens with a new core structure from non-emissive heteroaromatics | 적색 AIEgen 설계 — HAT 수치 없음(초록) |
| H-8 | 10.1107/s2056989018000725 | 2018 | Crystal structures of 2,3,8,9,14,15-hexamethyl-5,6,11,12,17,18-hexaazatrinaphthylene and 2,3,8,9,14,15-hexaphenyl-5,6,11,12,17,18-hexazatrinaphthylene dichloromethane disolvate | HATNMe6·HATNPh6 결정구조(π–π 적층 우세) — 광학/EA/증착 데이터 없음 |
| H-9 | 10.1039/c6tc04563d | 2017 | Ground-state charge transfer for NIR absorption with donor/acceptor molecules: interactions mediated via energetics and orbital symmetries | 공여/억셉터 바닥상태 CT 근적외 흡수 — 초록에 HAT 언급·수치 없음(HATCN 혼합 CT 흡수라면 주제 D 관련; PDF 확인 필요) |
| H-10 | 10.1002/anie.201500203 | 2015 | Carbazole Dendrimers as Solution‐Processable Thermally Activated Delayed‐Fluorescence Materials | 카바졸 덴드리머 TADF — HAT 무관(리뷰 인용 목록에만 존재) |
| H-11 | 10.1002/chem.201304071 | 2014 | Hexaazatrinaphthylenes with Different Twists | 꼬인 HATNA(Mateo-Alonso) — 초록에 수치 없음; 필요 시 PDF(uminho 저장소 요청형) |
| H-12 | 10.1016/j.synthmet.2014.06.013 | 2014 | Influence of the peripheral alkyl chains length on the mesomorphic properties of electron acceptor π-extended HAT derivatives | π-확장 HAT 억셉터 액정의 알킬 사슬 효과 — 초록 미확보(LC 중심) |
| H-13 | 10.1002/asia.201301522 | 2014 | Self‐Assembly of Chiral Propeller‐like Supermolecules with Unusual “Sergeants‐and‐Soldiers” and “Majority‐Rules” Effects | 키랄 프로펠러 초분자 — H5 무관 |
| H-14 | 10.1021/jp401919z | 2013 | The Impact of Local Work Function Variations on Fermi Level Pinning of Organic Semiconductors | ITO 위 F4TCNQ/HATCN 일함수·페르미 고정 — 주제 E 범위 |
| H-15 | 10.1002/adma.201302315 | 2013 | 25th Anniversary Article: Recent Advances in n‐Type and Ambipolar Organic Field‐Effect Transistors | n형 OFET 일반 리뷰 — HAT 유도체 수치 없음 |
| H-16 | 10.1039/c2cc36791b | 2013 | Linear and star-shaped naphthalimide-fused pyrazinacenes | 나프탈이미드 축합 피라지나센 — 초록에 수치 없음 |
| H-17 | 10.1039/c2nr30122a | 2012 | Self-assembly of conjugated oligomers and polymers at the interface: structure and properties | 계면 자기조립 일반 리뷰 — HAT 수치 없음 |
| H-18 | 10.1021/ol2030839 | 2012 | Large-Size Linear and Star-Shaped Dihydropyrazine Fused Pyrazinacenes | 별형 디하이드로피라진 피라지나센 — 초록에 수치 없음(응집 경향만) |
| H-19 | 10.1016/j.tetlet.2012.01.137 | 2012 | A hexaazatriphenylene-based organogel that responds to silver(I) with high selectivity under aqueous condition | Ag(I)에 선택 반응하는 HAT 유기겔 — 금속 결합 관련이나 초록 미확보(PDF needed, 우선순위 중) |
| H-20 | 10.1021/cm102772x | 2011 | One-Dimensional Nanostructures of π-Conjugated Molecular Systems: Assembly, Properties, and Applications from Photovoltaics, Sensors, and Nanophotonics to Nanoelectronics | 1D π-나노구조 일반 리뷰 — HAT 유도체 수치 없음 |
| H-21 | 10.1021/cm102117c | 2011 | Discotic Liquid Crystals for Opto-Electronic Applications | 원반형 액정 일반 리뷰(Kaafarani 2011) — HAT 유도체는 Segura 2015 리뷰로 대표 |
| H-22 | 10.1021/ja203776f | 2011 | Triangular Trinuclear Metal-N<sub>4</sub> Complexes with High Electrocatalytic Activity for Oxygen Reduction | 삼핵 금속-N4 ORR 촉매 — 광학/EA/증착 관련 데이터 없음 |
| H-23 | US7914908B2 | 2011 | Organic electroluminescent device having an azatriphenylene derivative | Kodak/Global OLED 특허 — 실시예 화합물이 트리아자트리페닐렌(TAPP, 6-phenyl-1,5,12-triazatriphenylene) ETL로 헥사아자(HAT)가 아님(참고: TAPP ETL 5.0 V, 9.2 cd/A vs Alq 7.3 V, Table 1) |
| H-24 | 10.1039/b925089a | 2010 | Alkyl substituted [6,6]-thienyl-C61-butyric acid methyl esters: easily accessible acceptor materials for bulk-heterojunction polymer solar cells | PCBM 유도체 — HAT 무관 |
| H-25 | 10.1016/j.tet.2010.09.064 | 2010 | Synthesis of polyphilic hexaazatrinaphthylenes and mesomorphic properties | 다친매성 HATNA 액정 — 액정상 중심, H5 관련 수치 없음 |
| H-26 | 10.3390/ijms11030943 | 2010 | Orientational Packing of a Confined Discotic Mesogen in the Columnar Phase | 제한 공간 원반 액정 배향 — H5 무관 |
| H-27 | 10.1002/adma.200801864 | 2009 | Microcontact Printing: Limitations and Achievements | 마이크로컨택 프린팅 리뷰 — HAT 물성 없음 |
| H-28 | 10.1016/j.synthmet.2009.05.004 | 2009 | The spectroelectrochemical, magnetic, and structural characterization of reduced hexaazatriphenylenehexacarbonitrile, HAT(CN)6 | 환원 HAT(CN)6 분광전기화학 — 초록 미확보(Konarev 2020 행 note에 언급) |
| H-29 | 10.1021/nl8021797 | 2008 | “Soft” Metallic Contact to Isolated C<sub>60</sub> Molecules | Segura 리뷰 참고문헌이지만 C60/금속 접촉 연구로 HAT 유도체 아님 |
| H-30 | 10.1063/1.3049616 | 2008 | Gold work function reduction by 2.2eV with an air-stable molecular donor layer | 분자 공여체로 Au 일함수 감소 — HAT 유도체 아님 |
| H-31 | 10.1063/1.2938082 | 2008 | Charge recombination in distributed heterostructures of semiconductor discotic and polymeric materials. | 원반형/고분자 이종구조 재결합 — HAT 유도체 물성 수치 없음 |
| H-32 | 10.1021/ja064934k | 2006 | Patterning of Robust Self-Assembled n-type Hexaazatrinaphthylene-Based Nanorods and Nanowires by Microcontact Printing | HATNA 나노로드 패터닝 — EA/흡수/증착 수치 없음 |
| H-33 | 10.1246/cl.2006.158 | 2006 | An Efficient Energy Transfer Found in Triphenylene/Hexaazatriphenylene System with Electronic Complimentarity and Structural Similarity | 트리페닐렌/HAT 에너지 전달 — 수치 확인 불가, H5 무관 |
| H-34 | 10.1021/ja054740p | 2005 | Dehydroaromatization of Quinoxalines:  One-Step Syntheses of Trinuclear 1,6,7,12,13,18-Hexaazatrinaphthylene Titanium Complexes | HATN 삼핵 Ti 착체 합성(유기금속) — 박막/광학 관련 없음 |
| H-35 | 10.1016/j.tet.2004.01.083 | 2004 | Hexaazatriisothianaphthenes: new electron-transport mesogens? | 헥사아자트리이소티아나프텐(티에노 축합 HAT) — 초록 미확보; Barlow 2007 티에노 HAT와 연결되는 후보(PDF needed, 우선순위 중) |
| H-36 | 10.1021/ja0390956 | 2004 | Charge Transport Properties in Discotic Liquid Crystals:  A Quantum-Chemical Insight into Structure−Property Relationships | 원반 액정 전하수송 계산(이동도) — EA/흡수 비교 없음 |
| H-37 | 10.1002/adma.200305137 | 2003 | Tailoring Discotic Mesophases: Columnar Order Enforced with Hydrogen Bonds | 수소결합 원반 액정 — HAT 물성 수치 없음 |
| H-38 | 10.1002/1439-7641(20020617)3:6<532::aid-cphc532>3.0.co;2-m | 2002 | Electron-Deficient Columnar Plastic Crystals | 전자결핍 컬럼형 플라스틱 결정(아자방향족 에스터) — 초록에 수치 없음, HAT 에스터 계열 여부 본문 확인 필요 |
| H-39 | 10.1039/b107135c | 2001 | Discotic mesogens with potential electron carrier properties | HATNA 원반 액정 최초 보고(초록 수치 없음) — Lehmann 2005 행 note에 언급 |
| H-40 | 10.1016/S0379-6779(01)00342-3 | 2001 | HAT(CN)6: a new building block for molecule-based magnetic materials | HAT(CN)6 분자 자성 빌딩블록 — 초록 미확보, 자성 중심 |
| H-41 | 10.1016/S0040-4039(98)02098-X | 1998 | First hyperpolarizabilities of hexaazatriphenylene derivatives: Octupolar nonlinear optical molecules | 팔중극 NLO HAT — 초록 미확보, 수치 확인 불가(Juárez 2010 행 note에 언급) |
| H-42 | 10.1055/s-1988-27470 | 1988 | Hexadecarboxylative Synthesis of Hexaazatriphenylene | HAT 헥사탈카복실화 합성(Sarma 1988) — 본문 미확보 |
| H-43 | 10.1021/jo00370a032 | 1986 | Improved synthesis of 1,4,5,8,9,12-hexaazatriphenylene | 무치환 HAT 개선 합성(Rogers 1986) — 본문/초록 미확보, 물성 수치 없음 |
| H-44 | 10.1016/S0040-4020(01)82926-4 | 1969 | Spectroscopic and structural studies of some oxocarbon condensation products—I | 옥소카본 축합 생성물 분광(1969) — 본문 미확보, HAT 계열 여부 불명 |

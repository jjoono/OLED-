# HANDOVER — near-unity EQE OLED 논문 작업 인수인계 (2026-09-23)

새 세션이 이 파일 하나로 지금까지의 작업·구조·결정을 파악할 수 있도록 압축한 요약. 상세 근거는 `manuscript/data_audit.md`(감사 기록 1–8차)와 `manuscript/en_v1_to_v2_changes.md`에 있다.

## 0. 사람·규칙

- 저자: Junho Kim (KAIST / Univ. of Cologne, Gather group). PI: Seunghyup Yoo (KAIST), Malte C. Gather (Cologne/St Andrews). 대화는 한국어.
- 저장소 `jjoono/OLED-`, 브랜치 `claude/oled-efficiency-paper-ag-955zsk`만 사용. 다른 브랜치에 push 금지, PR 만들지 않음. `git push -u origin <branch>`; 실패 시 2/4/8/16 s 재시도. 커밋 메시지에 모델명 넣지 않음(Co-Authored-By 줄은 시스템 리마인더대로).
- 결과물은 SendUserFile로 올린다(docx/pptx/xlsx). 사용자는 그림을 Origin으로 직접 그리고, 원고 편집은 Word에서 하며 PDF로 다시 올린다. 이 환경: LibreOffice Writer만 있음(Impress 없음 → pptx 렌더링 불가), Octave 있음, MATLAB 없음, pptxgenjs 없음(python-pptx 사용), 웹은 검색만 가능(science.org/PMC 차단).

## 1. 논문 한 줄 요약과 핵심 모델

- 제목(영문): **Towards near-unity OLEDs: Synergetic photon management**. 저자 Junho Kim¹˒², Seongjeong Lee², Seunghyup Yoo², Malte C. Gather¹˒³. 교신: Yoo, Gather.
- 주장: 외부 광추출 구조가 붙은 OLED에서는 빛이 광추출 구조와 소자 사이를 여러 번 왕복하며 추출되므로(**다중 통과 추출, multi-pass extraction**), 기판에 도달한 빛의 추출 효율 η_ext는 소자의 왕복 손실 A′(반사 전극의 Joule dissipation loss + 투명 전극 흡수)로 결정된다. 저손실 후면 반사체(Ag 또는 DBR) + 저흡수 투명 전극 + n_sub ≈ n_EML + SPP 억제가 설계 규칙. 실험: 일반 유리 + MLA 필름 + Ag → η_ext 91.6 % (EQE 20→55.8 %, 모델 55.9 %); 고굴절 MLA 기판 + Al + DBR 보조 반사체 → 61→77 %. 실제 스택 계산: 등방성 발광체 87–89 %, k_ITO 0.002면 89–91 %(표 1).
- 식: (1) η_EQE = η_int·η_out; (2) η_EQE = η_sub·η_ext; (3) η_ext = p/[p+(1−p)A′] (매 왕복 각도 무작위화 가정, T_LED = 0). 행렬 급수(Methods): P_out = P₀·P_sub·(I − B_R R)⁻¹·B_Tᵀ.
- **급수 vs 식 (3)**: 급수가 정확(측정 0.916 재현). 밀집 반구 MLA는 되돌리는 빛의 절반 이상을 60° 이상으로 보내 식 (3)보다 낮다. 급수는 P_sub(θ)를 그대로 쓰므로 n_sub를 연속으로 바꾸면 도파 모드가 기판에 수평으로 편입되는 굴절률(550 nm: n_eff 1.35/1.52/1.66/1.73, n_sub = n_org 1.8) 직후에 최대 10 %p 계단이 생긴다 → **그림 1, 2(c)–(f)는 식 (3)** (경향), **그림 3(a), 4, 5, 표 1은 급수**. 비교는 Supplementary Fig. 1(xlsx `SI_Fig1_nsub` 시트).
- p = 광추출 구조 고유량: BSDF 투과 성분 B_T(θ)를 cos θ sin θ로 가중 평균(n 1.3/1.5/1.8/2.0: 0.58/0.43/0.29/0.23).
- 터널링: n_sub > n_EML(이방성이면 TM의 n_e)이면 1 < u < n_sub/n_EML 근접장이 기판으로 전파 → η_sub에 포함(솔버 `cps2.solve_pol` 수정, 09-22). 표 1(EML n_e 1.671 < 1.8)과 2(d)–(f)의 n_sub > 1.8에만 영향.

## 2. 마스터 파일과 최신 산출물

| 무엇 | 경로 | 비고 |
|---|---|---|
| 한국어 마스터 | `manuscript/unityEQE_final_ko.md` | 모든 한국어 수정은 여기. 최신 docx `manuscript/unityEQE_v26_ko.docx` |
| 영문 마스터 | `manuscript/unityEQE_en.md` | 본문+Methods+캡션+표 1+참고문헌. 최신 docx `manuscript/unityEQE_en_v3.docx`(그림 포함) |
| 저자 Word 편집본 | 저자가 PDF로 올림(en_v1, en_v2, en_v3 PDF). 저자 Word 파일이 실제 작업본이므로 수정은 "붙여 넣을 문장"으로 전달 | `manuscript/upload_*` |
| docx 빌드 | `python3 manuscript/build_docx.py manuscript/upload_v6_user_20260921.docx <md> <out.docx> figures/final` | 4번째 인자(그림 폴더)를 주면 `**그림 N.**`/`**Figure N.**` 캡션 앞에 `figN.png` 삽입. 표는 8열 고정 너비 |
| 최종 그림 PNG | `figures/final/fig1..fig5.png` | 저자 v22 docx에서 추출(흰 배경 합성) |
| 표 1 슬라이드 | `figures/table1/table1.pptx` (`make_table1.py`) | "다파장 / 550 nm" 병기 |
| 그림 원자료 xlsx | `figures/update/figure_update_rawdata.xlsx` (`make_update_xlsx.py`) | 시트: README, Fig2c/2d/2e/2f(식 3), SI_Fig1_nsub, Fig3a_kITO/nAg, Fig5c, Fig5d, Fig5d_reference_lines, SI_Table1 |
| 감사 기록 | `manuscript/data_audit.md` | 1–8차 전수조사와 결정 |
| 변경 기록 | `manuscript/en_v1_to_v2_changes.md` | v1→v2, 저자 PDF→v2 |

## 3. 시뮬레이션 코드와 데이터 (Python)

- `sim/fig3c/cps2.py`: 1축 이방성 CPS 다이폴 솔버. `Stack(lam, eml(no,ne), d_eml, z0, above, below, n_sub, h)`, `solve_pol`(air/sub/wg/spp/abs + 터널링 이동), `stack_reflectance(S, θ)`, `sub_angular(S, θ, pols)`(단위 각도당 파워). `sim/fig3c/materials.py`: AG McPeak 0.04382+3.818978i, AL 0.958337+6.686780i, ITO König 1.86362+0.0032285i, ETL dict(B3PyMPM (1.8205, 1.6086)).
- `sim/mla/`: 자체 광선 추적기 `raytrace.py`(육방 최밀, cap/ellipsoid, PITCH 2.0), `run_aspect.py`(→ `mla_aspect.csv`, 5(d)), `series.py`(`eta_ext(BT, BR, R_LED, Psub, n_term)`; M = BR·R[:,None], 반사 후 각도에 R), `lt_hemisphere_bsdf.mat`(저자 LightTools 반구 BSDF, `BSDF_MLA` 180×90×15, n_MLA 1.30:0.05:2.00; 슬라이스 5 = 1.5, 11 = 1.8; 행 1–90 투과, 91–180 반사(역순으로 읽음), 열 = 입사각), `compare_lt_bsdf.py`(트레이서 ≡ LightTools, 소수 셋째 자리), `mla_aspect_lt_slice11.csv`.
- `sim/audit/recompute_series.py`: 2(c), 2(d)–(f), 3(a), 5(c), 표 1 재계산 → `fig2c_series.csv`, `fig2def_series.csv`, `fig3a_series.csv`, `fig5c_series.csv`, `table1_series.csv`(열: eta_sub, spp, wg, abs, Aprime, p_cos_sin, p_Psub, eta_ext_series, eta_ext_closed, EQE_series, EQE_closed). `recompute_fig2def_fine.py`(Δn 0.01, `graze70` 열), `fig2def_specavg.py`(스펙트럼 평균 참고용), `post_tunnel_fix.py`, `make_audit_xlsx.py`.
- `sim/table1_full/table1_full_spectrum.py`: 저자 라이브러리 분산(`nk_JH_total.mat`)으로 표 1을 400–700 nm(Ir(ppy)₂acac 광자 수 가중)와 550 nm로 계산 → `table1_full_spectrum.csv`. **표 1은 이 파일이 원천.**
- `sim/table1_author/`: 저자 모델 이식 `table1_author_model.m`(Octave 실행 가능), 저자 TMF 함수 3종(evanescent 분기 수정 포함, `TMF_SIGN` 환경변수 기본 +), `nk_JH_total.mat`, BSDF .mat, `test_psub_shape.m`, 결과 `table1_author_model.csv`.
- 실행: `python3 sim/audit/recompute_series.py`(수 분), `cd figures/update && python3 make_update_xlsx.py`(작업 폴더에서 실행해야 함), `cd figures/table1 && python3 make_table1.py`, 표 1 다파장 `python3 sim/table1_full/table1_full_spectrum.py`(4 프로세스, ~3 분).

## 4. 그림·표 상태 (저자 v3 PDF 기준)

| 그림 | 데이터 출처 | 상태 |
|---|---|---|
| 1(a–c) | 저자 모식도 + 식 (3) | 확인 완료. (b) 막대+누적선, 수렴 0.97/0.87/0.73; (c) p 0.25/0.4/0.6 |
| 2(a),(b) | 저자(A′ 분해) | 완료 |
| 2(c)–(f) | **식 (3)**, xlsx Fig2c/2d/2e/2f | 완료. (f) η_sub는 터널링 반영값 |
| 2(g)–(j) | 저자(각도·파장 1−R; 15.6/4.5/5.1 %) | 완료 |
| 3(a) | 급수, xlsx Fig3a_* | 데이터 완료. **모식도 "Glass (n = 1.77)" → 1.8 라벨 수정 남음(선택)** |
| 3(b),(c) | 저자 | 완료 |
| 4(a),(b) | 실측(저자) | 완료. 4(a) 음극 = Ag 100 nm(저손실 반사체), 발광 면적보다 넓은 Ag 막을 봉지층 위에 추가(전기적 분리) — 본문 문장 EN v3/KO v26에 반영 |
| 5(a),(b) | 저자 MATLAB(400–700 nm) | 완료. (b) 최대 87.4 % |
| 5(c) | 급수, xlsx Fig5c | **완료(31→55 %, 83→85 %; v3에서 확인)** |
| 5(d) | 자체 트레이서 급수, xlsx Fig5d | 완료. Ag 최대 86.6 %@AR 0.55, 0.25–1.5에서 95 % 이상; 평면 Al 21.6/Ag 25.5 % |
| 5(e) | 저자 산란층 MC(HG, S 1–15, g 0.5–1.0, 400–800 nm, PLQY 0.98, Θ 0.76, n_sub 1.8로 기재) | 완료 |
| 표 1 | `table1_full_spectrum.csv` | 완료. 각 칸 "다파장 / 550 nm" |

## 5. 본문에 인용된 핵심 수치 (검증 완료)

- 초록·서론·논의: 91.6 %, 55.8/55.9 %, 61→77 %, "80 % 후반", "90 %".
- 2(a–c): Al 반사체 손실 13 %/왕복, ITO 150 nm는 A′의 15 %; Ag 반사체 1.8 %, ITO가 절반 이상; Ag–Al 11–13 %p; 30→200 nm에서 Ag 97→93 %, Al 84→82 %.
- 2(d–f): p 0.58→0.23(2.5배); 격차 7/12/22 %p(1.3/1.5/2.0); Al 1.8에서 61 % 최대 후 59–61 %; Ag 1.8에서 89 %, 2.0에서 88 %.
- 2(g–j): Al 15.6 %, Ag 4.5 %, DBR 5.1 %(흡수 4.1 + 투과 1.0). Supplementary Note 2: λ/4 스택 9.5 %, 균일 10/15/20쌍 5.5/5.7/6.0 %, chirp ZnS 55→77, LiF 91→127 nm.
- 3(a): k_ITO 0→0.01: η_ext 95→81, EQE 91→75; k 0.08: 45 %; n_Ag 0.2: 55 %. 3(b): 10 % 이하 320 nm, 1 % 500 nm. 3(c): 20 % 이하 두께 110/104/90/59 nm(n_e 1.8/1.7/1.6/1.5).
- 4: 20→55.8 %(2.8배), η_sub 60.9 %, 91.6 %; orange 23/61/77 %, η_sub 88 %, η_ext 87.5 %, DBR 4.5쌍 ZnS 70/LiF 115 nm.
- 5: (b) 87.4 %, ETL ≥ 150 nm에서 최댓값의 90 %; (c) 31→55 %, 83→85 %, 등방성 84 %; 실제 스택 87–89 %(+3–5 %p), Θ 0.8–0.9도 87–89 %, k 0.002면 89–91 %; (d) 86.6 %, 0.25–1.5; (e) 80–90 %.
- 표 1(다파장/550): 등방성 0.870/0.878(ETL 200), 0.883/0.890(300); k 0.002: 0.889–0.896/0.897–0.905(200), 0.901–0.905/0.907–0.912(300). 식 (3)은 550 nm에서 η_ext 0.932–0.945, EQE 0.898–0.931(급수보다 1.4–3.3 %p 높음).
- Methods: p 0.43(1.5)/0.29(1.8); 두 구현(Python vs 저자 MATLAB) 일치 η_sub 0.2 %p, η_ext 1 %p; ITO König 550 nm 1.864+0.0032i, 400–700 nm n 1.79–2.08; 5(a),(b) ETL/HTL 10–300 nm, 400–700 nm, 격자 단위당 1,000점(총 3,000).

## 6. 용어 통일 (한국어 원고)

설계 규칙(디자인 룰·design rule 금지), 반사체(반사판·거울 금지; 금속 반사체의 Joule dissipation loss), 다중 통과 추출(광자 재활용 금지), Joule dissipation loss(ohmic 손실 금지), 광추출(띄어쓰기 없음), 투명 전극, 마이크로캐비티/캐비티 조건·차수, 도파 모드(waveguided mode), 기판 전달 효율 η_sub·기판 모드, sub-to-air 추출 효율 η_ext, 파워, 인광 OLED, microlens array, 발광체, 임계각, 왕복 손실, 허용 범위, escape cone, BTDF/BRDF, evanescent 성분(SPP 모드 포함).

## 7. 저자 스크립트에서 발견한 문제 (저자에게 전달됨, 재실행 요청 상태)

1. TMF_birefringence_whole*.m(2016 원본): evanescent 분기 버그(복소 배열의 −0 허수부) → `imag(n·cosθ) ≥ 0` 강제 필요. 미수정 시 η_sub > 1, SPP < 0.
2. `Planar_sweep22_preprint_MLA.m`의 급수 입력 `Psub_norm = I_sub.*sin089`: P_sub_ang이 이미 단위 각도당 파워(균질 매질에서 sinθ)이므로 sin 이중 가중 → η_ext 약 1.6 %p 과소.
3. 같은 스크립트 282–283행 `ROLED_3 = repmat(ROLED,1,1,90)`: R_LED를 반사 전 각도에 곱함 → `permute(…,[1 3 2])`.
→ 그림 4의 모델 EQE(55.9 %), 5(a),(b), 5(e) 재실행 시 1–2 %p 상승 가능. "0.1 %p 이내 일치" 문장은 재실행 후 갱신 필요. 저자는 아직 재실행 결과를 주지 않음.

## 8. 남은 항목 (preprint 기준 "반드시")

1. 참고문헌 [[ ]] 7곳: [12] 마이크로캐비티 고전, [17] TCO/금속 흡수, [18] CPS/다이폴(Neyts 1998/Furno 2012), [19] Poynting vector, [20] 통합 모델(그룹 선행 논문), [22],[23] 측정 장비. [25] Kim, H. S. et al. Sci. Adv. 9, eadf1388 (2023); [26] Kim, H. S., Cheon, H. J. et al. Sci. Adv. 11, eadr1326 (2025)는 채움.
2. 저자 Word 파일에 EN v3의 green 소자 문장(실험 검증 절 1·2문단, 논의 첫 문장, 그림 4 캡션) 반영 여부 확인.
3. 저자 Word 파일: 표 1 Methods 문장 "at 550 nm" → "at 550 nm and over the 400–700 nm emission spectrum"; "ε_e = 2nk" → "ε₂ = 2nk".
4. 선택: 그림 3(a) 모식도 라벨 1.77 → 1.8; Methods "and its angular distribution"의 its; "normalized" → "normalised".
5. 보충자료(Note 1–2, Fig. 1–6)는 본문에서 인용 중이나 미작성(저자가 "일단 패스"). Note 2 초안과 Fig. 1·2 데이터는 저장소에 있음(`manuscript/unityEQE_final_ko.md` 끝부분, xlsx SI_Fig1_nsub, `sim/mla/lt_hemisphere_bsdf.mat`, `figures/fig2_full/fig2j_rawdata.xlsx`).
6. 저자 재실행(7절) 결과가 오면 그림 4·5 관련 수치와 "0.1 %p" 문장 갱신.

## 9. 대화 흐름 요약 (결정 순서)

1. 저자 v6 docx + 그림 pptx → 최종 한국어 원고 작성, 7항목 리뷰(`v6_final_review.md`).
2. 5(d) AR 곡선: 저자 LightTools 반구 BSDF ≡ 자체 트레이서 검증; 식 (3)은 AR ≥ 1에서 최대, 급수는 AR 0.55 → 급수 채택, 허용 범위 문구 유지.
3. 전수조사(`data_audit.md`): 식 (3)/급수 혼용 정리, p 정의 통일, 3(a) 상수 1.86/3.82, Table 1 신설.
4. 라인 에딧 배치, 용어 변경(다중 통과 추출, Joule dissipation loss), UDC 문장 삭제, McPeak [24].
5. 2(e) 급수 계단 발견 → 원인(모드 컷오프, 수평 편입) 설명 → 저자 결정: 2(c)–(f) 식 (3) 복귀, 비교는 SI Fig. 1.
6. 저자 v20/v22(그림 삽입) 검토, 표 1 pptx, 용어·표현·순서 정리(v21), 터널링 분류 수정(v22), ITO 상수 König 통일.
7. 저자 답변 반영(제목·저자·문헌·Methods 조건·증착), 영문 초안 v1 → 저자 PDF 편집 반영 v2(BTDF/BRDF, 400–700 nm, 캡션 문구) → 표 1 다파장 → v3(green 소자 반사체 설명, ε₂).
8. 저자 라이브러리로 저자 모델을 Octave에서 실행, 문제 3건 발견, 두 구현 비교.

## 10. 새 세션 시작 시 권장 순서

1. `git log --oneline | head`, 이 파일, `manuscript/data_audit.md` 8차 항목 읽기.
2. 저자가 올린 최신 PDF/docx가 있으면 `manuscript/unityEQE_en.md`와 문장 단위 diff(스크립트 예: pdftotext → 문장 분할 → difflib; 수식 개체는 서식 차이로 무시).
3. 수정은 md 마스터(KO/EN 둘 다)에 반영 → build_docx.py로 docx 재생성 → 커밋·푸시 → SendUserFile.

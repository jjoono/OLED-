# 그림·표별 원자료 안내 (Optics Express 투고본)

2026-11-22 기준. 본문과 Supplement 1의 모든 그림·표가 어떤 파일에서 나왔는지 정리한다. 그림에 실제로 찍힌 값은 `Figure_source_data.xlsx`에 패널별 시트로 들어 있고(본문 Table 1–2 포함, 총 55개 시트), 각 시트의 출처는 그 파일의 `README` 시트에 적혀 있다.

- 저장소에 있는 것: 그림 원자료 워크북(`Figure_source_data.xlsx`, `figures_oe/*.xlsx`), 반환 행렬(`multipass_reset/**/*.npz`), 소자·렌즈 모델 결과(`labdev/*.json`, `labdev/kernel*.npz`), LightTools 캠페인 기록(루트의 `*.mat`), 모든 스크립트.
- 저장소에 없는 것: `selectivity/data/`(구조 라이브러리 평가 결과, .gitignore 대상)와 측정 원자료. 이 둘은 별도 묶음 `OE_raw_data.zip`으로 전달한다. `selectivity/data/`는 아래 "다시 만드는 순서"로 재생성할 수 있다.

## 본문 그림

| 그림 | 내용 | 원자료 시트 | 계산 결과 파일 | 스크립트 |
|---|---|---|---|---|
| Fig. 1(b)–(d) | 평면, 반구 MLA, 산란층의 반환 행렬 $M$ | `Fig1b_planar_M`, `Fig1c_hemiMLA_M`, `Fig1d_scatter_M` | `multipass_reset/random_sweep/*.npz`, `selectivity/data/slabs.npz` | `selectivity/oe_figs.py` |
| Fig. 1(e), (f) | 패스별 입사 각도 분포, 패스별 탈출 확률 | `Fig1e_pass_resolved`, `Fig1f_escape_per_pass` | 위와 같음 + `selectivity/data/devices.npz` | `selectivity/oe_figs.py` (`passes()`) |
| Fig. 2(a) | $p_{\mathrm{mix}}$ 대 $p$ | `Fig2a_p_mix` | 구조 라이브러리(위) | `selectivity/oe_figs.py` |
| Fig. 2(b) | 소자 혼합과 탈출면 혼합의 차이 | `Fig2b_placement` | `selectivity/data/placement.json` | `selectivity/placement.py` |
| Fig. 2(c) | 굴절률별 슬랩, 텍스처 출사면 | `Fig2c_slab`, `Fig2c_textured` | `placement.json`(slab), `selectivity/data/highindex.json`, `multipass_reset/highindex/*.npz` | `selectivity/placement.py`, `multipass_reset/run_highindex.py`, `selectivity/highindex_slabs.py` |
| Fig. 2(d) | Henyey–Greenstein 부피 | `Fig2d_HG_volume` | `selectivity/data/scatter_volume.json` | `selectivity/scatter_volume.py` |
| Fig. 2(e) | 완전 혼합의 이득 대 $N_{\mathrm{mix}}$ | `Fig2e_gain` | 구조 라이브러리 | `selectivity/oe_figs.py`, `selectivity/reciprocal_limit.py` |
| Fig. 3(a), (b) | 축약 모델 대 전체 급수(342개) | `Fig3a_reduced_reciprocal`, `Fig3b_error_points`, `Fig3b_error_bars`, `Fig3_all_columns` | 구조 라이브러리 | `selectivity/oe_figs.py`, `selectivity/mixing.py` |
| Fig. 3(c) | 적용 범위, 사전 시험 | `Fig3c_applicability` | `selectivity/data/cutoff_modes.json`, `prospective_predictions.json`, `prospective_results.json`, `multipass_reset/prospective/*.npz`, `multipass_reset/blind/*.npz` | `selectivity/cutoff_modes.py`, `multipass_reset/prospective.py`, `selectivity/prospective_test.py` |
| Fig. 4 | $(p, N_{\mathrm{mix}})$ 지도, 무작위화, $\xi$ | `Fig4a_structures`, `Fig4a_contours`, `Fig4b_randomness`, `Fig4c_xi` | 구조 라이브러리 | `selectivity/oe_figs.py` |
| Fig. 5(b)–(e) | 제작 소자와 유한 기하 모델 | `Fig5bc_profiles`, `Fig5d_angular_factor`, `Fig5e_ratios` | `labdev/gonio_az_orange.json`, `gonio_az_e65h80.json`, `layout_results.json`, `green_patch.json`, `kernelabs_*.npz` | `labdev/fig_experiment.py`, `gonio_az.py`, `trace_abs.py`, `layout.py`, `green_patch.py` |

## 본문 표

| 표 | 내용 | 원자료 | 스크립트 |
|---|---|---|---|
| Table 1 | 혼합 위치별 탈출 확률과 평면 계면의 $\eta_{\mathrm{sta}}$ | 시트 `Table1`, `selectivity/data/placement.json` | `selectivity/placement.py` |
| Table 2 | 제작 소자의 각도 인자(모델 범위, 측정) | 시트 `Table2`, `Fig5d_angular_factor` | `labdev/fig_experiment.py`, `labdev/gonio_az.py` |
| (구 Table 3) | 견고성 점검 요약은 Supplement 1의 Table S13으로 옮김 | 각 행의 출처는 표의 Source 열 | 아래 SI 표 항목 참조 |

## Supplement 1의 그림

| 그림 | 내용 | 원자료 시트 | 계산 결과 파일 | 스크립트 |
|---|---|---|---|---|
| Fig. S1 | 소자 입력(첫 패스 휘도, 반사율) | `FigS1_device_inputs` | `selectivity/data/devices.npz` | `selectivity/fig_devices.py` |
| Fig. S2 | 더 단순한 모델 대 전체 급수 | `FigS2a_…`, `FigS2b_…`, `FigS2c_…` | 구조 라이브러리 | `selectivity/oe_figs.py` |
| Fig. S3 | 유한 기판 | `FigS3a_size`, `FigS3b_pixel`, `FigS3c_mixing_gain` | `selectivity/data/finite_size.json` | `selectivity/finite_size.py`, `selectivity/fig_finite_size.py` |
| Fig. S4 | 휘도 순서 상한, 출력의 각도 구성 | `FigS4a_first_pass`, `FigS4b_composition` | 구조 라이브러리 | `selectivity/oe_figs.py` |
| Fig. S5 | 자유형 대 반구 벤치마크 | `FigS5a_profiles`, `FigS5b_best_EQE`, `FigS5c_searches` | 루트의 LightTools 캠페인 기록(`*.mat`) | `make_figures.py` |

## Supplement 1의 표

| 표 | 내용 | 계산 결과 파일 | 스크립트 |
|---|---|---|---|
| Table S1 | 각도 절단에 따른 축약 모델 오차 | `selectivity/data/cutoff_modes.json` | `selectivity/cutoff_modes.py` |
| Table S2 | 몬테카를로 불확도 | `multipass_reset/random_sweep_unc/*.npz`, `selectivity/data/uncertainty.json` | `selectivity/unc_analysis.py` |
| Table S3 | 소자 적층 | `selectivity/data/devices.npz` | `selectivity/fig_devices.py`, `selectivity/run_devices.py` |
| Table S4, S5 | 38개 구조, 종류별 오차 | `figures_oe/fig4_structure_map.xlsx`, `fig3_reduced_model.xlsx` | `selectivity/oe_figs.py` |
| Table S6 | 직접 광선추적, 표본 외 구조 | `multipass_reset/blind/*.npz`, `results.json` | `multipass_reset/direct_blind.py`, `blind_table.py` |
| Table S7 | 이완 속도와 사전 시험 | `selectivity/data/cutoff_modes.json`, `prospective_predictions.json`, `prospective_results.json` | `selectivity/cutoff_modes.py`, `multipass_reset/prospective.py`, `selectivity/prospective_test.py` |
| Table S8 | HG 부피 | `selectivity/data/scatter_volume.json` | `selectivity/scatter_volume.py` |
| Table S9 | 소자·파장별 완전 혼합 이득 | `selectivity/data/mixing_gain_by_condition.json`, `mixing_gain_spectral.json` | `selectivity/make_note_s12.py` |
| Table S10 | 스펙트럼 가중, 혼합 위치 세 가지 | `selectivity/data/spectral_placement.json` | `selectivity/spectral_placement.py` |
| Table S11 | 텍스처 출사면 슬랩 | `selectivity/data/highindex.json`, `multipass_reset/highindex/*.npz` | `multipass_reset/run_highindex.py`, `selectivity/highindex_slabs.py` |
| Table S12 | 광 포획 점검 | `selectivity/data/solar_trapping.json` | `multipass_reset/solar_trapping.py` |
| Table S13 | 견고성 요약 | 각 행의 출처는 표의 Source 열 | — |
| Table S14 | 자유형 대 반구 | LightTools 캠페인 기록(루트 `*.mat`) | `make_supp_figures_and_data.py` |
| Table S15, S16 | 측정 소자(주황, 녹색) | 측정 원자료(별도 묶음) | 요약값은 `labdev/fig_experiment.py` 상단 |
| Table S17 | 렌즈 소자 EQE 차이의 원인 후보 | `labdev/lens_loss_test.json`, `labdev/layout_variants*.json` | `labdev/lens_loss_test.py`, `labdev/layout_variants.py` |
| Table S18 | 렌즈 패치 폭 | `labdev/layout_sweep.json` | `labdev/layout_sweep.py` |

노트 본문에만 나오는 수치의 출처: 구간 폭·상관계수 `selectivity/data/rev3_checks.json`(`rev3_checks.py`), 편광 `multipass_reset/polcheck.py`, 슬랩 직접 추적 `selectivity/data/slab_direct_check.json`, 소자 혼합과 평면 슬랩 `selectivity/data/placement.json`, 산란층 위치 `selectivity/data/device_side_slabs.json`, 고정 캐비티 렌즈 최적화 `multipass_reset/fixedcav/results.json`, 고니오미터 민감도 `labdev/gonio_az_*.json`.

압축 전 Supplement(표 40개, 그림 10개)에 있던 나머지 표와 그림의 자료는 저장소에 그대로 있으며, 2026-11-26 이전 커밋의 이 문서에 대응표가 있다.

## 다시 만드는 순서

```
cd multipass_reset && python3 run_random_sweep.py && cd ..
cd selectivity && python3 extract_inputs.py && python3 run_devices.py && python3 oe_figs.py && cd ..
python3 make_figures.py && python3 make_supp_figures_and_data.py && python3 make_optimizer_convergence.py
cd labdev && python3 fig_experiment.py && cd ..
python3 make_source_data.py            # Figure_source_data.xlsx
python3 tools/cite.py tools/manuscript_en.src.md manuscript_draft.md
python3 md2optica.py && node md2docx.js en && node md2docx.js ko && node md2docx.js supp
```

## 측정 원자료의 범위

본문과 Supplement에는 측정값의 요약 통계(소자별 EQE, 각도 인자, 정규화 각도 프로파일)만 들어 있다. 계측기에서 나온 원본 파일은 저장소에 두지 않았고, 이 세션에 제공된 자료는 `OE_raw_data.zip`의 `3_측정_원자료/`에 정리했다. 소자 제외 기준(S4D3, G-S5D3, 반구 스캔 7)은 Supplement 1의 Note S13에 적혀 있다.

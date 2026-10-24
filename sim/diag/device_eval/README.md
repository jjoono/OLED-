# 측정 소자 성능 분석: afterAgMLA (2026-10-23, green Ir(ppy)₂acac)

`dev3.py`의 측정값 `meas = [0.16, 0.51, 0.59]` 중 **MLA + Ag 후면 반사체(0.59)** 에 해당하는 소자를
측정 원본에서 다시 평가한 결과다. 평가는 OLED-evaluation(GatherLab)의 계산 함수를 그대로 쓰되,
각도 보정만 −90~90° 전체를 쓰도록 바꿨다.

## 요약

| | 각도 분포 d44_p1 | 각도 분포 d45_p3 |
|---|---|---|
| 최대 EQE (−90~90°) | **57.9 %** (45 cd/m², 3.3 V) | **55.3 %** (46 cd/m², 3.3 V) |
| 최대 EQE (0~+90°만, GUI 방식) | 61.0 % | 57.6 % |
| EQE @ 100 / 1000 cd/m² | 57.5 / 54.7 % | 54.9 / 52.2 % |
| CE 최대 / PE 최대 | 202 cd/A / 254 lm/W | 204 cd/A / 247 lm/W |
| PE @ 100 / 1000 cd/m² | 176 / 107 lm/W | 171 / 104 lm/W |
| 전압 @ 100 / 1000 cd/m² | 3.79 / 5.95 V | 3.78 / 5.93 V |
| 각도 보정 계수 e (+쪽 / −쪽) | 0.535 (0.564 / 0.507) | 0.510 (0.532 / 0.488) |

공통: JVL은 d44_p3 하나, 최대 휘도 약 2,100 cd/m²(7.0 V, 1.13 mA/cm²), turn-on 2.0~2.5 V,
CIE (0.29, 0.65), 피크 523 nm, FWHM 52 nm. roll-off는 100 → 1000 cd/m²에서 약 5 %.

**`dev3.py`의 0.59와의 관계**: 0.59는 GUI 방식(0~+90°만)으로 낸 값과 가깝다(57.6~61.0 %).
−90~90°를 모두 쓰면 이 소자는 **55~58 %** 이고, 어느 각도 측정을 쓰느냐로 약 5 % 달라진다.
모델 비교에 쓸 값으로는 각도 분포 d44_p1 쪽(57.9 %)이 측정 상태가 더 좋다(아래 3).

## 측정 조건과 계산

- 픽셀 4 mm², PD–소자 58 mm, PD gain 50 dB(475 kΩ), PD 면적 75 mm², `pd_cutoff` 1×10⁻⁴ V
- 분광기 보정 `Calibration_HDX_QP600_F950SMA-A_2025-09-26` (스펙트럼 범위 344~931 nm와 일치.
  더 새 파일 `Calibration_NR_QP600-2_F950SMA-C_2026-07-14`는 898~1706 nm 근적외 분광기용)
- JVL: Accurate 모드(chop, 목표 정밀도 3 %). PD σ 열로 EQE 1σ(대부분 ±1~2 %)를 냈다
- 각도 분포: gon-spec-deg(열화 보정본), −90~90° 10° 간격, 0° 스펙트럼을 수직 스펙트럼으로 사용

**GUI와 다른 점**
1. 각도 보정 계수를 0~+90°와 −90~0°에서 같은 함수(`calculate_e/v_correction`)로 각각 계산해 평균했다.
   GUI는 0~+90°만 쓴다. 이 소자는 좌우가 비대칭이라 −90~90° 값이 GUI보다 EQE 4~5 % 낮다.
2. **evaluation 코드 버그**: `calculate_v_correction`의 음의 각도 분기에 `abs()`가 빠져 있어(EQE 쪽
   `calculate_e_correction`에는 있음) 휘도 보정 계수가 음수로 나온다. 여기서는 절댓값으로 바로잡았다.
   PE에만 영향. GUI에서 −90~0°만 있는 파일을 평가하면 PE가 틀린다.

## 데이터 파일에서 확인한 문제

1. **원래 d45_p3로 알려진 JVL은 d44_p3 측정이다.** 측정 PC의 `this!2/..._d45_p3_jvl.csv`와
   `this!/..._d44_p3_jvl.csv`는 바이트까지 같고(20:16 생성, σ 열 삭제), 원본은 `this!/격리/..._d44_p3_jvl.csv`(20:14, σ 포함)다.
   여기에는 원본을 `data/2026-10-23_afterAgMLA_d44_p3_jvl.csv`로 넣었다. d45_p3 각도 측정이 같은 픽셀인지는 저자 확인 필요.
2. **`data/..._d44_p1_jvl.csv`(같은 픽셀 d44_p1의 JVL)는 쓸 수 없다.** 3.7~3.8 V에서 전류가 정확히 10.000 µA로
   고정되고 PD가 거의 0으로 떨어지며, 3.9~4.2 V도 전류가 앞뒤와 맞지 않는다(EQE 32~38 %).
   Accurate 모드가 SMU 레인지를 10 µA로 고정하면서 전류 상한(ILimit)이 낮아지고 복구되지 않은 것으로 보인다.
   측정 코드 수정 확인용 증거로만 넣었다.
3. **각도 측정 상태**
   - d44_p1: 50 µA 설정이지만 모든 각도에서 전압이 7.000 V 한계에 걸렸다(실제 전류는 약 43 µA로 추정).
     측정 중 0° 밝기 변화 1.05×, 0.97×, 좌우 대칭 양호.
   - d45_p3: 30 µA, 6.30~6.32 V. 측정 중 0° 밝기가 1.02×, 1.11×로 올라 덜 안정적이고, 좌우 비대칭이 더 크다.
4. turn-on 전압은 측정 격자(2.0 V 다음이 2.5 V) 때문에 2.0~2.5 V로만 알 수 있다.
5. 절대 EQE는 PD 거리 기준점에 민감하다(±1 mm → 약 ±3.5 %).

## 파일

| 경로 | 내용 |
|---|---|
| `../data/2026-10-23_afterAgMLA_d44_p3_jvl.csv` | 분석에 쓴 JVL 원본 (σ 포함) |
| `../data/2026-10-23_afterAgMLA_d44_p1_gon-*.csv`, `..._d45_p3_gon-*.csv` | 각도 분해 스펙트럼 원본·열화 보정본·각도별 구동 전압 |
| `../data/2026-10-23_afterAgMLA_d44_p1_jvl.csv` | 문제 2의 증거 (분석에 쓰지 않음) |
| `pixel_report.py` | 한 픽셀 평가 → xlsx(README, Summary, JVL, Angular, Spectrum) + png 8패널 |
| `run_afterAgMLA.py` | 위 두 경우를 계산하고 `out/afterAgMLA_summary.csv` 생성 |
| `out/<case>/*_performance.{xlsx,png}` | 결과 (xlsx Summary의 최대값은 JVL 시트 참조 수식) |
| `oled_evaluation/` | GatherLab OLED-evaluation 일부 사본 (MIT, `SOURCE_COMMIT.txt`): `evaluation_functions.py`, 보정 파일 4개, `photodiode_gain.json` |

재실행: `cd sim/diag/device_eval && python3 run_afterAgMLA.py` (numpy, pandas, scipy, openpyxl, matplotlib).
다른 픽셀: `python3 pixel_report.py <jvl.csv> <gon-spec-deg.csv> <out_dir> "<이름>"`.

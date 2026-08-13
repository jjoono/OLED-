# 적용 패치 모음

`../README.md` 진단 결과를 실제 코드에 반영하기 위한 파일들이다.
**아직 아무것도 적용하지 않은 상태**이며, 적용 순서와 전제 조건이 있다.

기준 소스는 GitHub upstream (GatherLab/OLED-jvl-measurement, GatherLab/OLED-evaluation)
최신 커밋이다. 랩 로컬 사본(`<데이터폴더>\Software\`)은 upstream 과 몇 군데 다르므로
아래 "로컬 사본과의 차이" 항목을 먼저 확인할 것.

| 파일 | 대상 저장소 | 관련 단계 |
|---|---|---|
| `hardware.diff` | OLED-jvl-measurement `src/hardware.py` | 2-A(레인지 전환), 4(버스트 판독), 5(SMU 설정) |
| `chopped_measurement.py` | OLED-jvl-measurement `src/` 에 새 파일로 추가 | 4 |
| `evaluation_functions.diff` | OLED-evaluation `src/evaluation_functions.py` | 6 |
| `photodiode_gain.json` | OLED-evaluation `usr/photodiode_gain.json` 교체 | 8 |

---

## 적용 순서

### 지금 바로 적용 가능 (벤치 테스트 결과와 무관)

**1. `photodiode_gain.json`** — 40/60 dB cutoff 가 실측 노이즈의 44~48배로 잘못 들어가 있고
모노토닉하지도 않았다(40 dB 300 uV > 50 dB 100 uV). 다른 항목과 같은 기준(실측 노이즈의 5 sigma)으로 고쳤다.

| Gain | 기존 | 수정 | 근거(실측 노이즈) |
|---|---|---|---|
| 40 dB | 3e-4 | **3e-5** | 6.8 uV x 5 |
| 60 dB | 4e-3 | **4e-4** | 83 uV x 5 |

50 dB(현재 유일한 실사용 gain) 값은 랩에서 2023-03-30 에 손으로 바꾼 `1e-4` 를 그대로 유지했다
(upstream 기본값은 `2e-4`). **파일을 통째로 덮어쓰기 전에 로컬 값이 정말 1e-4 인지 확인할 것.**
80 dB 항목은 저항값이 60 dB 보다 작게 적혀 있어 수상하지만, 근거 데이터가 없어 손대지 않았다.

**2. `hardware.py` 오타** — 로컬 사본에만 있는 문제라 diff 에 넣지 않았다. 직접 고칠 것:
```
self.spectrometerr.close()   ->   self.spectrometer.close()
```
`OceanSpectrometer.close_connection()` 안에 있다(로컬 477행 근처). 지금은 spectrometer close 가
항상 AttributeError 로 실패한다.

### 벤치 테스트(1단계) 결과 — 2026-08, 완료

| 테스트 | 조건 | 평균 | sigma |
|---|---|---|---|
| ① | DMM 입력 단락, 10 V 레인지 | 0.2 uV | **1.02 uV** |
| ② | 앰프 연결 + PD 암전, 10 V 레인지 | 8198 uV | **28.62 uV** |
| ③ | ②와 같은 상태, 100 mV 레인지 | 8218 uV | **28 uV** |

```
앰프 기여분 = sqrt(28.62^2 - 1.02^2) = 28.60 uV
DMM 기여분 =                            1.02 uV   (분산 기준 0.13%)
```

**앰프(PDA200C) 지배로 확정.** 결과:

- **레인지 전환(2-A)은 폐기.** ①②에서 예측한 값이 28.60 uV 였고 ③ 실측이 28 uV 로,
  n=200 에서 sigma 추정 불확실도(약 +-5% = +-1.4 uV) 안에서 구분되지 않는다. 얻는 것이 0 이다.
  `hardware.diff` 의 `set_fixed_range` 는 적용은 하되 쓰지 않는다(앰프 교체 후 재확인용).
  ③에서 8.2 mV 오프셋이 100 mV 레인지에 overload 없이 들어간 것은 확인됐으므로,
  TIA 교체로 노이즈가 내려간 뒤에는 레인지를 내리는 선택지가 살아 있다.
- **저잡음 TIA 교체가 유일한 하드웨어 레버.** 475 kOhm 의 Johnson 노이즈 한계(0.089 uV/sqrt(Hz))보다
  약 42배 위에 있고, PD 산탄노이즈는 완전히 무시 가능하다(1 cd/m2 신호 기준 0.04 pA).
  10~30배 개선 시 측정 시간은 sigma^2 에 비례하므로 100~900배 단축된다.
- **측정 시간 추정 상향**: 원 추정이 sigma = 20 uV 가정이었으므로 (28.62/20)^2 = 2.05배.
  1 cd/m2 +-2% 가 50초 -> 약 100초 (PD 거리 60 mm 적용 후 기준).
- **후속 테스트 추가**: 앰프 노이즈가 백색인지 1/f 드리프트가 섞였는지에 따라
  `baseline_period` 와 `samples_per_burst` 가 달라진다.
  `../bench_test/noise_test2_averaging.py` (테스트 ②와 같은 셋업, 물리 작업 없음, 약 3분).

### 적용

**3. `hardware.diff`** — `git apply hardware.diff` (OLED-jvl-measurement 루트에서). 내용:

- `KeithleyMultimeter.set_fixed_range` / `set_auto_range` **주석 해제 + 문법 수정**.
  원래 주석 코드는 `CONF:VOLTage:DC:RANGe <v>` 를 썼는데 2100 에서 유효하지 않다.
  `SENSe:VOLTage:DC:RANGe <v>` 로 바꿨다(이걸 쓰면 auto range 는 자동으로 꺼진다).
  → 벤치 테스트가 **DMM 지배**로 나올 때 쓰는 무비용 개선(2-A).
- `KeithleyMultimeter.configure_burst` / `read_burst` **신규**.
  `CONFigure` 1회 + `SAMPle:COUNt N` + `READ?` 로 한 번에 N 샘플을 받는다.
  기존 `measure_voltage()` 는 호출마다 `MEASure?` 를 보내 기기를 재설정하므로 chopped 스캔에 쓸 수 없다.
- `KeithleySource.set_low_light_mode(on)` **신규** — 5단계. 저휘도 모드에서만
  `Current:AZero ON` + 전류 NPLC 10, 종료 시 원복. 이걸 빼면 PD 를 고쳐놔도 분모(전류)가 EQE 오차를 지배한다.

**4. `chopped_measurement.py`** — `src/` 에 복사. `hardware.diff` 가 먼저 적용돼 있어야 한다.
기존 `AutotubeMeasurement`(일반 JVL 스캔)은 건드리지 않는다.

- 포인트마다 `V_target(settle -> 버스트) / 0 V(settle -> 버스트)` 를 반복하고 `dPD = mean(on) - mean(off)`
  → background 를 매 사이클 다시 재므로 **원인 (2) 계통 오프셋이 구조적으로 소멸**
- 반복 평균으로 **원인 (1) 랜덤 노이즈가 sqrt(N) 으로 감소**
- **N 고정이 아니라 수렴 조건**: `sigma/dPD <= target_relative_sigma` (기본 0.02 = +-2%).
  밝은 포인트는 자동으로 몇 사이클 만에 끝난다.
- **baseline 공유**: `baseline_period` 초마다 한 번만 off 를 재고 인접 사이클이 공유 → 약 2배 단축.
  드리프트가 의심되면 `baseline_period=0` 으로 두면 매 사이클 측정한다(가장 안전, 2배 느림).
- 포인트별 **sigma 를 파일에 같이 저장**(`Photodiode Voltage Sigma` 컬럼).

불확실도 계산은 baseline 공유 여부에 따라 자동으로 맞게 처리한다(`ChoppedSweep._combine`):
같은 baseline 을 쓴 사이클들을 한 그룹으로 묶고 **그룹 평균들의 산포**로 sigma 를 낸다.
전체 delta 의 std/sqrt(N) 을 그대로 쓰면 공유 baseline 오차를 과소평가하고,
baseline 오차를 따로 더하면 매 사이클 재측정하는 경우 같은 오차를 두 번 세게 되기 때문이다.

하드웨어 없이 로직 확인:
```
python chopped_measurement.py --simulate
```
가짜 DMM(신호 + 20 uV 노이즈 + 앰프 DC 오프셋)으로 수렴 동작과 통계를 검증한다.
1 cd/m2 (PD 40 uV) 포인트가 +-2% 에 231 사이클 = 실기기 기준 약 55초로 수렴하며,
이는 `../README.md` §8 의 시간 추정과 일치한다.

GUI 에 붙일 때는 `LowLuminanceMeasurement` (QThread 래퍼) 를 쓰면 되고,
`measurement_parameters` 에 다음 키가 필요하다:
`voltages`, `scan_compliance`, `samples_per_burst`, `settle_time`, `multimeter_range`,
`target_relative_sigma`, `baseline_period`.

`multimeter_range` 는 **10 (V) 로 둔다**. 2026-08 벤치 테스트에서 노이즈가 앰프 지배로 확정됐다
(아래 참조). 저잡음 TIA 로 교체한 뒤에 다시 판단할 것.

**5. `evaluation_functions.diff`** — `git apply evaluation_functions.diff` (OLED-evaluation 루트에서).

`JVLData.__init__` 의 cutoff 처리를 바꾼다. **기존 일반 스캔 동작은 그대로**이고,
데이터에 `pd_voltage_std` 컬럼이 있을 때(= chopped 스캔)만 다르게 동작한다:

```
pd_cutoff 적용 안 함  ->  dPD < 3 sigma 인 포인트만 NaN
```

핵심은 **0 과 NaN 의 구분**이다.
- `0` = "빛이 없다" (물리적 결론)
- `NaN` = "이 포인트는 측정 불가" (정보 없음)

기존 cutoff 는 후자를 전자로 기록해서, 최근 소자 기준 약 30 cd/m2 이하가 통째로
`L=0.00 / EQE=0.00` 으로 남는 원인이었다.

`sigma_threshold` 는 생성자 인자(기본 3.0)라 호출부 수정 없이 그대로 쓸 수 있다.

---

## 로컬 사본과의 차이 (적용 전 확인)

랩 로컬 사본은 upstream 과 최소한 아래가 다르다. diff 가 깨지면 이것부터 의심할 것.

- `hardware.py` 의 `spectrometerr` 오타 (upstream 에는 없음)
- `photodiode_gain.json` 의 50 dB cutoff `1e-4` (upstream `2e-4`)
- 진단 문서의 행 번호는 로컬 사본 기준이라 upstream 과 2~3행 어긋난다

diff 가 적용되지 않으면 `git apply --3way` 또는 해당 함수를 손으로 옮겨 넣으면 된다.
바꾸는 부분이 함수 단위로 독립적이라 수동 이식도 어렵지 않다.

## 이 패치들이 하지 않는 것

- **PD 거리 167 -> 60 mm (3단계)**: 하드웨어 작업 + 거리 재보정이라 코드로 대체 불가.
  chopped 스캔만으로도 오차는 줄지만, 거리를 줄여 신호를 7.7배 키우기 전에는
  1 cd/m2 에서 필요한 사이클 수가 수십 배로 늘어난다.
- **PD gain 전환**: 50 dB 고정 유지. 70 dB 는 신호 10배 대비 노이즈 12배(실측 245 uV)로 손해다.
- **일반 JVL 스캔 변경**: 그대로 둔다. 위 변경은 전부 저휘도 전용 경로에만 적용된다.

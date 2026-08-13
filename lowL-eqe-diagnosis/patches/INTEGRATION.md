# chopped 스캔 적용 절차

GUI 는 나중에 건드린다. **먼저 단독 스크립트로 돌려서 측정이 실제로 되는지 확인**하고,
그 다음에 GUI 에 붙인다. 순서를 지키면 중간에 문제가 생겨도 원인이 한 군데다.

측정 소프트웨어는 **완전히 종료**한 상태로 진행한다 (장비를 스크립트가 직접 잡아야 함).

---

## 0단계 — 되돌릴 수 있게 만들기 (5분)

```powershell
cd C:\Users\GatherLab\Documents\Software\OLED-jvl-measurement
git status                     # git 저장소면 현재 상태 확인
git checkout -b lowlight-scan  # 작업 브랜치
```

git 저장소가 아니면 폴더를 통째로 복사해 둔다:

```powershell
Copy-Item -Recurse ..\OLED-jvl-measurement ..\OLED-jvl-measurement_backup
```

---

## 1단계 — `hardware.py` 에 메서드 3종 추가 (10분)

`src/hardware.py` 를 연다. 두 군데에 코드를 붙여 넣는다.

### [A] `KeithleyMultimeter` 클래스

`measure_voltage()` 정의 **바로 위**에, 주석 처리돼 있던
`# def set_fixed_range` / `# def set_auto_range` 블록을 **지우고** 그 자리에
`hardware_snippets.txt` 의 [A] 덩어리를 붙여 넣는다.

추가되는 것:
- `set_fixed_range` / `set_auto_range` — 주석 해제 + 문법 수정
  (원래 코드의 `CONF:VOLTage:DC:RANGe` 는 2100 에서 유효하지 않다)
- `configure_burst` / `read_burst` — `CONFigure` 1회 + `SAMPle:COUNt` + `READ?` 버스트 판독.
  **chopped 스캔에 필수다.** 기존 `measure_voltage()` 는 호출마다 `MEASure?` 로
  기기를 재설정해서 느리다.

### [B] `KeithleySource` 클래스

`as_current_source()` 정의 **바로 위**에 [B] 덩어리를 붙여 넣는다.
`set_low_light_mode(on)` 이 추가된다 — 저휘도 모드에서만 전류 autozero ON + NPLC 10,
끄면 원복. 이걸 빼면 PD 를 고쳐놔도 분모(전류)가 EQE 오차를 지배한다.

### 확인

```powershell
python -c "import sys; sys.path.insert(0,'src'); import hardware; print([m for m in dir(hardware.KeithleyMultimeter) if 'burst' in m or 'range' in m]); print(hasattr(hardware.KeithleySource,'set_low_light_mode'))"
```

`['configure_burst', 'read_burst', 'set_auto_range', 'set_fixed_range']` 와 `True` 가 나오면 된다.

> `git apply ..\hardware.diff` 로 한 번에 적용해도 된다. 로컬 사본이 upstream 과
> 달라서 실패하면(오타 `spectrometerr` 등) 위 수동 절차를 쓴다.

---

## 2단계 — 파일 2개 복사 (1분)

`chopped_measurement.py` 와 `run_lowlight_scan.py` 를 `src/` 에 넣는다.

```powershell
Copy-Item chopped_measurement.py, run_lowlight_scan.py src\
```

하드웨어 없이 먼저 로직만 확인해 볼 수 있다:

```powershell
cd src
python chopped_measurement.py --simulate
```

---

## 3단계 — SMU 정착 시간 실측 (10분)

`settle_time` 은 chopped 사이클 시간의 절반 가까이를 차지한다. 기본값 20 ms 가
과한지 실측해서 정한다. **소자 하나를 실제로 켜는 첫 측정이다.**

```powershell
python run_lowlight_scan.py --measure-settle --pixel 1 --voltage 3.5
```

전압 스텝 직후 PD 값을 연속으로 찍어서, 최종값 대비 몇 %인지 표로 보여준다.
목표 정밀도(±2%) 안에 처음 들어오는 시점이 필요한 `settle_time` 이다.
첫 샘플부터 안정적이면 5 ms 까지 줄여도 된다 → 측정 시간이 1.5배 빨라진다.

---

## 4단계 — chopped 스캔 실행 (측정 시간에 따라 수 분)

밝은 포인트 하나로 먼저 짧게 확인한다:

```powershell
python run_lowlight_scan.py --pixel 1 --voltages 4.0 --settle 0.005
```

수십 사이클 안에 수렴하고 `converged=True` 가 나오면 정상이다. 그 다음 본 스캔:

```powershell
python run_lowlight_scan.py --pixel 1 --voltages 2.6,2.8,3.0,3.3,3.6,4.0 --settle 0.005
```

전압값은 소자의 turn-on 에 맞춰 조정한다. 목표는 1 / 2 / 5 / 10 / 30 / 100 cd/m² 근처다.
결과 파일은 `default_saving_path` 에 `lowlight_pixelN_날짜.txt` 로 저장된다.

### 보는 것

| 열 | 의미 |
|---|---|
| `pd_voltage` | chopped 로 구한 순수 신호 (배경 제거됨) |
| `pd_voltage_std` | 그 포인트의 실측 불확실도. **이게 핵심이다** |
| `cycles` | 수렴까지 걸린 사이클 수 |
| `converged` | False 면 신호가 노이즈에 묻힌 것 |

### 기대치 (PD 52 mm 기준)

1 cd/m² 에서 PD 신호가 약 52 µV 이므로, ±2% 에 약 440 사이클 / 약 40초
(settle 5 ms 면 약 26초). 10 cd/m² 이상은 1~2초에 끝난다.

### 비교 실험 (선택)

```powershell
python run_lowlight_scan.py --pixel 1 --voltages 3.0 --pattern ABBA --settle 0.005
```

같은 포인트를 `AB` 와 `ABBA` 로 돌려서 같은 정밀도에 걸리는 시간(`elapsed`)을 비교한다.
드리프트에 선형 성분이 있으면 ABBA 가 빠르다.

---

## 5단계 — 기존 스윕과 비교 검증 (중요)

**밝은 구간에서는 두 방법이 같은 답을 줘야 한다.** 이게 chopped 스캔이 맞게 동작한다는 증거다.

100 cd/m² 이상이 나오는 전압 2~3개를 골라 기존 GUI 스윕과 chopped 스캔을 각각 돌리고
`pd_voltage` 를 비교한다. 몇 % 안에서 일치하면 검증 완료다.
저휘도에서 갈라지는 것은 정상이다 — 그게 원래 고치려던 문제다.

---

## 6단계 — evaluation 쪽 (평가 소프트웨어)

`evaluation_functions.diff` 를 OLED-evaluation 에 적용한다.
`pd_voltage_std` 컬럼이 있는 파일(= chopped 스캔 결과)만 다르게 처리된다.
기존 파일은 종래 동작 그대로다.

```
pd_cutoff 적용 안 함  ->  ΔPD < 3σ 인 포인트만 NaN
```

`photodiode_gain.json` 교체와 `spectrometerr` 오타 수정도 이때 같이 한다 (README.md 참고).

---

## 7단계 — GUI 연결 (나중에, 선택)

4~5단계에서 측정이 검증된 뒤에만 한다. `chopped_measurement.py` 의
`LowLuminanceMeasurement` (QThread) 를 쓰면 되고, `AutotubeMeasurement` 과 시그널 구성이 같다.
`measurement_parameters` 에 필요한 키:

```
voltages, scan_compliance, samples_per_burst, settle_time,
multimeter_range, target_relative_sigma, baseline_period
```

**기존 일반 JVL 스캔은 건드리지 않는다.** 별도 탭/버튼으로 붙인다.

---

## 문제가 생기면

| 증상 | 원인 |
|---|---|
| `AttributeError: configure_burst` | 1단계 [A] 가 안 붙었거나 다른 클래스에 붙음 |
| 판독이 전부 overload | `--range` 를 10 으로 (앰프 DC 오프셋이 8.2 mV) |
| 모든 포인트 `converged=False` | 전압이 너무 낮아 발광이 없음. 전압을 올릴 것 |
| 사이클이 계속 돌고 안 끝남 | 정상일 수 있다. `--max-time` 으로 상한을 건다 |
| Arduino/VISA 연결 실패 | 측정 GUI 가 아직 떠 있는지 확인 |

스크립트는 정상 종료든 Ctrl+C 든 예외든 **출력을 끄고 릴레이를 전부 내린다.**

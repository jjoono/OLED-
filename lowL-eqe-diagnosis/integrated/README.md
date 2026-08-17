# Accurate EQE Mode 통합본

기존 GUI 소프트웨어에 **"Accurate EQE Mode" 토글**을 넣은 전체 코드다.
upstream (GatherLab/OLED-jvl-measurement, OLED-evaluation) 최신 커밋 기반.

- **토글 OFF (기본)**: 기존과 완전히 동일하게 측정된다. 파일 포맷도 그대로 3컬럼.
- **토글 ON**: 같은 전압 리스트를 chopped 방식으로 측정한다.
  포인트마다 소자를 켰다/껐다 반복하며 dPD = mean(on) − mean(off) 를 누적하고,
  불확실도가 **Target Precision (%)** 이하로 떨어지면 다음 포인트로 넘어간다.
  파일에 `Photodiode Voltage Sigma` 4번째 컬럼이 추가된다.

autotube 탭에 위젯 2개가 추가된다: `Accurate EQE Mode` 토글, `Target Precision (%)` 스핀박스(기본 2%).

## 파일

| 경로 | 내용 |
|---|---|
| `jvl-measurement/` | 측정 소프트웨어 수정본 5개 (전체 파일) |
| `evaluation/` | 평가 소프트웨어 수정본 2개 (전체 파일) |
| `diffs/` | upstream 대비 diff (로컬 사본에 git apply 용) |
| `tests/` | Qt/장비 없이 도는 스모크 테스트 |

### jvl-measurement (5개)

- `UI_main_window.py` — 토글 + 정밀도 스핀박스 (autotube 탭 14/15행)
- `main.py` — 기본값 설정 + `read_autotube_parameters` 에 파라미터 2개 추가
- `autotube_measurement.py` — 측정 루프 분기. OFF 경로는 기존 코드 그대로 유지
- `chopped_measurement.py` — **새 파일.** chopped 엔진 (patches/ 의 것과 동일)
- `hardware.py` — 버스트 판독(`configure_burst`/`read_burst`), 레인지 전환 복원, `set_low_light_mode`

### evaluation (2개)

- `core_functions.py` — `read_multiple_files` 가 파일 헤더에서 Sigma 컬럼을 감지.
  **이게 없으면 4컬럼 파일을 3컬럼으로 읽으면서 모든 컬럼이 한 칸씩 밀린다(조용히!).**
  일반 파일과 accurate 파일이 섞여 있어도 된다.
- `evaluation_functions.py` — `pd_voltage_std` 가 있는 데이터는 pd_cutoff 대신
  **ΔPD < 3σ → NaN** 판정. 일반 데이터는 기존 cutoff 동작 그대로.

## ON 상태에서 정확히 무엇이 달라지나

| 항목 | OFF (기존) | ON (accurate) |
|---|---|---|
| background | sweep 시작 시 1샘플, 전 포인트에서 뺌 | 매 사이클 재측정 (계통 오프셋 소멸) |
| PD 판독 | 포인트당 1샘플 | 수렴까지 반복 평균 (기본 ±2%) |
| SMU | NPLC 1, autozero OFF | **NPLC 10, autozero ON** (종료 시 원복) |
| DMM | `MEASure?` 매번 재설정 | `CONFigure` 1회 + `READ?` 버스트 |
| 파일 | 3컬럼 | 4컬럼 (+σ), 헤더에 모드 표기 |
| 부가 진단 | — | 포인트별 trend (측정 중 소자 변화 감지), 수렴 실패 로그 |

중지 버튼은 ON 상태에서도 즉시 듣는다 (진행 중인 포인트가 다음 사이클에서 멈춤).
측정이 어떤 경로로 끝나도(정상/중지/예외) SMU·DMM 설정은 원복된다.

## 겸사겸사 고친 upstream 버그 4개

1. `run()` 첫머리의 `import pydevd` — 개발용 잔재. pydevd 미설치 환경에서 측정 시작 즉시 죽는다
2. 컴플라이언스 도달 시 `self.log_message` 호출 — 존재하지 않는 메서드라 AttributeError
3. **픽셀 간 df_data 미초기화** — 앞 픽셀이 더 길면 남은 꼬리 행들이 다음 픽셀 파일에 섞여 저장됨
4. **저장 시 df 를 제자리에서 문자열로 변환** — 이후 플롯/다음 픽셀이 문자열 위에서 돌게 됨 (신형 pandas 에서는 크래시)

## 설치

측정 GUI 완전 종료 후:

```powershell
cd C:\Users\GatherLab\Documents\Software\OLED-jvl-measurement
git checkout -b accurate-eqe          # 또는 폴더 백업
# 방법 A: 전체 파일 교체 (로컬 수정이 없거나 확인한 경우)
Copy-Item <이폴더>\jvl-measurement\* src\ -Force
# 방법 B: 로컬 수정을 보존해야 하면
git apply <이폴더>\diffs\jvl-measurement.diff   # 실패 시 --3way
```

evaluation 쪽도 동일하게 (`evaluation/` 2개 파일 또는 `diffs/evaluation.diff`).
`patches/photodiode_gain.json` 교체와 로컬 `spectrometerr` 오타 수정은 별도 (patches/README.md).

> **주의**: 로컬 사본은 upstream 과 다르다 (`spectrometerr` 오타 등). 전체 파일 교체(방법 A)는
> 로컬 수정을 덮어쓴다. 로컬에서 뭘 바꿨는지 모르면 먼저 diff 로 확인할 것:
> `git diff --no-index 로컬src upstream클론src`

이 통합본은 `patches/evaluation_functions.diff` 를 **대체한다** (같은 내용 + 리더 수정 포함).
`patches/run_lowlight_scan.py` (단독 실행 스크립트)는 통합본과 무관하게 계속 쓸 수 있다.

## 검증한 것

- 5개 파일 문법 검사 통과
- **스모크 테스트** (`tests/`, Qt/장비 스텁): OFF/ON 모두 2픽셀 전체 플로우 통과 —
  ON 에서 prepare→포인트별 수렴 로그→4컬럼 저장→finish 원복, OFF 에서 기존 3컬럼 저장
- **evaluation 파이프라인**: 3/4컬럼 혼합 읽기에서 컬럼 밀림 없음, accurate 행은 3σ 판정
  (약한 신호 NaN 확인), 일반 행은 기존 cutoff 동작
- 실기기 검증은 아직이다. **INTEGRATION.md 의 5단계(밝은 구간에서 기존 스윕과 일치 확인)를 반드시 할 것**

재실행: `cd tests && python smoke_test.py` (numpy, pandas 만 필요)

## GUI 사용법

1. 저휘도 전압 구간을 촘촘히: Min 2.4 / Changeover 3.2 / low step 0.05~0.1 권장
   (기존 0.5 V 스텝이면 1–10 cd/m² 구간에 포인트가 거의 없다)
2. `Accurate EQE Mode` ON, Target Precision 2%
3. 픽셀 선택 후 Start. 로그에 포인트별 `dPD ± σ | N cycles` 가 찍힌다
4. "not converged" 가 뜨는 포인트는 신호가 노이즈에 묻힌 것 — 정상이며, evaluation 에서 NaN 처리된다
5. 수명 짧은 소자는 max voltage 를 낮게 (측정 시간의 대부분은 저휘도 포인트가 아니라
   고전압 구간이 소자를 상하게 한다)

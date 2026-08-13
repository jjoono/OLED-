# -*- coding: utf-8 -*-
"""
저휘도(1~100 cd/m2) 전용 chopped JVL 스캔.

OLED-jvl-measurement/src/ 에 그대로 넣고 쓰는 파일이다. 기존 AutotubeMeasurement
(일반 JVL 스캔) 는 건드리지 않는다. 이 모드는 저휘도 정밀 측정에만 쓴다.

왜 필요한가 (진단 결과 요약, 자세한 건 ../README.md):
  원인 (1) 포인트당 PD 1샘플 -> sigma = 20 uV 의 랜덤 산포
  원인 (2) sweep 시작 시 background 를 1회만 재고 전 포인트에서 빼기 때문에,
           그 1샘플에 실린 +-20 uV 가 곡선 전체의 계통 오프셋이 된다.
           EQE_측정 = EQE_참값 * (1 + delta/L) 이라 저휘도에서 크게 부푼다.

  이 파일이 하는 일:
    - 전압 포인트마다 소자를 켰다/껐다 반복하며 on-off 차이를 취한다
      -> background 를 매 사이클 다시 재는 셈이므로 원인 (2) 가 구조적으로 사라진다
    - 사이클을 반복 평균한다 -> 원인 (1) 이 sqrt(N) 으로 줄어든다
    - 포인트별 측정 불확실도 sigma 를 파일에 함께 저장한다
      -> evaluation 쪽에서 pd_cutoff 로 0 을 만드는 대신 "유의하지 않음 = NaN" 판정이 가능해진다
       (evaluation_functions.diff 참고)

전제 조건:
  - hardware.diff 가 먼저 적용돼 있어야 한다
    (KeithleyMultimeter.configure_burst / read_burst / set_fixed_range,
     KeithleySource.set_low_light_mode)
  - 벤치 테스트(bench_test/noise_test.py) 결과에 따라 multimeter_range 를 정한다.
    DMM 지배로 나오면 0.1 (100 mV), 앰프 지배로 나오면 10 을 그대로 두고
    앰프 교체 후 다시 내린다.

하드웨어 없이 로직만 확인:
    python chopped_measurement.py --simulate
"""

import time
import datetime as dt

import numpy as np
import pandas as pd

try:
    from PySide6 import QtCore
except ImportError:  # 시뮬레이션/단위테스트용
    QtCore = None

try:
    import core_functions as cf

    log_message = cf.log_message
except ImportError:

    def log_message(message):
        print(message)


# ---------------------------------------------------------------------------
# 측정 엔진 (Qt 비의존, 하드웨어 객체만 있으면 단독 실행 가능)
# ---------------------------------------------------------------------------


class ChoppedSweep:
    """
    한 픽셀에 대한 chopped 스캔.

    keithley_source      : hardware.KeithleySource (voltage source 모드)
    keithley_multimeter  : hardware.KeithleyMultimeter (configure_burst 적용본)
    """

    def __init__(
        self,
        keithley_source,
        keithley_multimeter,
        samples_per_burst=5,
        settle_time=0.02,
        multimeter_range=0.1,
        nplc=1,
        target_relative_sigma=0.02,
        min_cycles=16,
        max_cycles=2000,
        max_time_per_point=240.0,
        baseline_period=5.0,
        scan_compliance=None,
    ):
        """
        target_relative_sigma:
            포인트별 수렴 조건. sigma(dPD) / dPD 가 이 값 이하가 되면 멈춘다.
            0.02 (= +-2%) 가 일상 검증용 권장값. 논문 대표소자만 0.01.
            +-1% -> +-2% 로 완화하면 측정 시간이 4배 줄고 (1 cd/m2 기준 3.5분 -> 50초)
            "1/10/100 cd/m2 효율이 일치하는가" 판정 결론은 달라지지 않는다.

        baseline_period:
            off(baseline) 를 매 사이클 재지 않고 이 시간(초)마다 한 번만 다시 잰다.
            인접 포인트가 baseline 을 공유하므로 약 2배 빨라진다. 드리프트가
            의심되면 0 으로 두면 매 사이클 측정(가장 안전, 대신 2배 느림).

        min_cycles:
            수렴 판정을 시작하기 전 최소 사이클 수. 표본이 적으면 sigma 추정 자체가
            널뛰어서 "운 좋게 조용한 몇 사이클" 로 조기 종료해 버린다. 최소 16 사이클,
            그리고 baseline 그룹 4개 이상을 모은 뒤에만 판정한다.

        max_cycles / max_time_per_point:
            수렴하지 않는 포인트(신호가 노이즈에 묻힌 경우)에서 무한정 도는 것을 막는다.
            이 한계에 걸려 끝난 포인트는 converged=False 로 기록된다.
        """
        self.keithley_source = keithley_source
        self.keithley_multimeter = keithley_multimeter

        self.samples_per_burst = int(samples_per_burst)
        self.settle_time = settle_time
        self.multimeter_range = multimeter_range
        self.nplc = nplc

        self.target_relative_sigma = target_relative_sigma
        self.min_cycles = int(min_cycles)
        self.max_cycles = int(max_cycles)
        self.max_time_per_point = max_time_per_point
        self.baseline_period = baseline_period
        self.scan_compliance = scan_compliance

        # 공유 baseline 상태
        self._baseline_mean = None
        self._baseline_sem = None
        self._baseline_time = None
        self._baseline_id = -1

        self.stop = False

    # -- 저수준 -------------------------------------------------------------

    def prepare(self):
        """저휘도 모드 진입. 반드시 finish() 와 짝으로 쓸 것."""
        self.keithley_multimeter.configure_burst(
            self.samples_per_burst,
            multimeter_range=self.multimeter_range,
            nplc=self.nplc,
        )
        self.keithley_source.set_low_light_mode(True)
        log_message(
            "저휘도 모드 진입: DMM range "
            + str(self.multimeter_range)
            + " V, burst "
            + str(self.samples_per_burst)
            + " samples, SMU autozero ON / NPLC 10"
        )

    def finish(self):
        """저휘도 모드 종료. SMU/DMM 설정을 일반 스캔 조건으로 되돌린다."""
        self.keithley_source.set_low_light_mode(False)
        self.keithley_multimeter.reset()
        log_message("저휘도 모드 종료: SMU/DMM 설정 원복")

    def _read_pd_burst(self):
        """PD 전압 버스트 판독. overload 는 NaN 으로 바꿔서 돌려준다."""
        values = np.asarray(self.keithley_multimeter.read_burst(), dtype=float)
        # 2100 은 overload 를 9.9e37 로 반환한다
        values[np.abs(values) > 1e30] = np.nan

        if np.all(np.isnan(values)):
            raise RuntimeError(
                "PD 판독이 전부 overload 다. multimeter_range 가 앰프 DC 오프셋보다 "
                "작다는 뜻이므로 한 단계 올릴 것 (0.1 -> 1)."
            )
        return values

    def _measure_baseline(self, force=False):
        """
        off 상태(0 V) 의 PD 전압을 측정한다.

        baseline_period 안에 이미 잰 값이 있으면 그대로 재사용한다(baseline 공유).
        """
        now = time.time()
        if (
            not force
            and self._baseline_mean is not None
            and self.baseline_period > 0
            and (now - self._baseline_time) < self.baseline_period
        ):
            return

        self.keithley_source.set_voltage(str(0))
        time.sleep(self.settle_time)
        values = self._read_pd_burst()

        self._baseline_mean = float(np.nanmean(values))
        n = int(np.sum(~np.isnan(values)))
        self._baseline_sem = (
            float(np.nanstd(values, ddof=1)) / np.sqrt(n) if n > 1 else 0.0
        )
        self._baseline_time = now
        self._baseline_id += 1

    @staticmethod
    def _combine(deltas, groups):
        """
        누적한 사이클들에서 (평균, 평균의 불확실도) 를 계산한다.

        baseline 을 공유하면 같은 baseline 을 쓴 사이클들의 오차가 공통 성분을
        갖는다. 그래서 전체 delta 의 표준오차 std/sqrt(N) 을 그대로 쓰면 공유
        baseline 의 오차를 과소평가한다. 반대로 baseline 오차를 따로 제곱합으로
        더하면, baseline 을 매 사이클 다시 재는 경우에는 같은 오차를 두 번 세게 된다.

        그래서 baseline 그룹별 평균을 하나의 독립 측정으로 보고 그룹 평균들의 산포로
        불확실도를 낸다. 이러면 두 경우가 자동으로 맞는다.
          - baseline_period = 0 (매 사이클 재측정) -> 그룹 = 사이클, 표준 결과와 동일
          - baseline 공유                          -> 공통 오차가 그룹 간 산포에 그대로 반영
        그룹이 4개 미만이면 아직 그룹 산포를 못 믿으므로 보수적인 쪽(제곱합)으로 낸다.
        """
        group_means = [float(np.mean(values)) for values in groups.values()]
        k = len(group_means)

        if k >= 4:
            mean = float(np.mean(group_means))
            sigma = float(np.std(group_means, ddof=1)) / np.sqrt(k)
            return mean, sigma

        mean = float(np.mean(deltas))
        n = len(deltas)
        sem = float(np.std(deltas, ddof=1)) / np.sqrt(n) if n > 1 else np.nan
        return mean, float(np.sqrt(np.nan_to_num(sem, nan=np.inf) ** 2))

    # -- 포인트 단위 --------------------------------------------------------

    def measure_point(self, voltage):
        """
        전압 포인트 하나를 chopped 방식으로 측정한다.

        on / off 를 반복하며 dPD = mean(on) - mean(off) 를 누적하고,
        sigma(dPD) 가 목표 상대오차 이하가 되면 멈춘다.

        반환: dict
          voltage, pd_voltage, pd_voltage_std, current, current_std,
          cycles, elapsed, converged
        """
        deltas = []
        groups = {}  # baseline id -> 그 baseline 을 쓴 사이클들의 delta
        currents = []
        start = time.time()
        compliance_hit = False

        while True:
            # --- off (baseline) ---
            self._measure_baseline()

            # --- on ---
            self.keithley_source.set_voltage(str(voltage))
            time.sleep(self.settle_time)
            on_values = self._read_pd_burst()
            oled_current = float(self.keithley_source.read_current())

            delta = float(np.nanmean(on_values)) - self._baseline_mean
            deltas.append(delta)
            groups.setdefault(self._baseline_id, []).append(delta)
            currents.append(oled_current)

            if (
                self.scan_compliance is not None
                and abs(oled_current) >= self.scan_compliance
            ):
                compliance_hit = True
                break

            n = len(deltas)
            elapsed = time.time() - start

            if n >= self.min_cycles:
                delta_mean, sigma = self._combine(deltas, groups)
                if abs(delta_mean) > 0 and (
                    sigma / abs(delta_mean) <= self.target_relative_sigma
                ):
                    break

            if n >= self.max_cycles or elapsed >= self.max_time_per_point:
                break
            if self.stop:
                break

        n = len(deltas)
        delta_mean, sigma = self._combine(deltas, groups)
        converged = (
            not compliance_hit
            and abs(delta_mean) > 0
            and np.isfinite(sigma)
            and sigma / abs(delta_mean) <= self.target_relative_sigma
        )

        return {
            "voltage": voltage,
            "pd_voltage": delta_mean,
            "pd_voltage_std": sigma,
            "current": float(np.mean(currents)) * 1e3,  # mA
            "current_std": (
                float(np.std(currents, ddof=1)) / np.sqrt(n) * 1e3
                if n > 1
                else np.nan
            ),
            "cycles": n,
            "elapsed": time.time() - start,
            "converged": bool(converged),
            "compliance": bool(compliance_hit),
        }

    # -- sweep --------------------------------------------------------------

    def run_sweep(self, voltages, progress_callback=None):
        """전압 리스트를 순서대로 측정하고 DataFrame 으로 반환한다."""
        rows = []
        self.keithley_source.activate_output()
        try:
            for index, voltage in enumerate(voltages):
                if self.stop:
                    break
                row = self.measure_point(voltage)
                rows.append(row)

                log_message(
                    "V = %.3f V | dPD = %.1f uV +- %.1f uV (%.1f%%) | "
                    "%d cycles, %.1f s%s"
                    % (
                        row["voltage"],
                        row["pd_voltage"] * 1e6,
                        row["pd_voltage_std"] * 1e6,
                        (
                            abs(row["pd_voltage_std"] / row["pd_voltage"]) * 100
                            if row["pd_voltage"]
                            else float("nan")
                        ),
                        row["cycles"],
                        row["elapsed"],
                        "" if row["converged"] else "  <- 수렴 실패",
                    )
                )

                if row["compliance"]:
                    log_message("Current compliance reached")
                    break
                if progress_callback is not None:
                    progress_callback(int((index + 1) / len(voltages) * 100))
        finally:
            # 어떤 경로로 빠져나가도 소자에 전압이 남지 않게 한다
            self.keithley_source.set_voltage(str(0))
            self.keithley_source.deactivate_output()

        return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 저장
# ---------------------------------------------------------------------------


def save_data(df_data, file_path, setup_parameters, measurement_parameters):
    """
    기존 autotube 파일과 같은 형식으로 저장하되 sigma 컬럼을 추가한다.

    evaluation 쪽은 'pd_voltage_std' 컬럼이 있으면 chopped 데이터로 인식하고
    pd_cutoff 대신 3 sigma 유의성 판정을 쓴다 (evaluation_functions.diff).
    """
    header_lines = [
        "Chopped low-luminance JVL scan",
        "Date: " + dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "Samples per burst: " + str(measurement_parameters["samples_per_burst"]),
        "Settle time: " + str(measurement_parameters["settle_time"]) + " s",
        "Multimeter range: " + str(measurement_parameters["multimeter_range"]) + " V",
        "Target relative sigma: "
        + str(measurement_parameters["target_relative_sigma"]),
        "Baseline period: " + str(measurement_parameters["baseline_period"]) + " s",
        "SMU low light mode: autozero ON, current NPLC 10",
        "### Measurement data ###",
        "OLEDVoltage\t OLEDCurrent\t Photodiode Voltage\t Photodiode Voltage Sigma",
        "V\t mA\t V\t V",
    ]

    with open(file_path, "w", encoding="utf-8") as handle:
        for key, value in setup_parameters.items():
            handle.write(str(key) + ": " + str(value) + "\n")
        handle.write("\n".join(header_lines) + "\n")
        df_data[["voltage", "current", "pd_voltage", "pd_voltage_std"]].to_csv(
            handle, index=False, header=False, sep="\t", lineterminator="\n"
        )

    log_message("저장 완료: " + str(file_path))


# ---------------------------------------------------------------------------
# GUI 연결용 QThread 래퍼 (AutotubeMeasurement 과 같은 시그널 구성)
# ---------------------------------------------------------------------------

if QtCore is not None:

    class LowLuminanceMeasurement(QtCore.QThread):
        """
        AutotubeMeasurement 과 같은 자리에서 쓰는 저휘도 전용 스캔 스레드.
        선택된 픽셀들을 돌면서 ChoppedSweep 를 실행한다.
        """

        update_plot = QtCore.Signal(list, list, list)
        update_progress_bar = QtCore.Signal(str, float)
        hide_progress_bar = QtCore.Signal()
        reset_start_button = QtCore.Signal(bool)

        def __init__(
            self,
            keithley_source,
            keithley_multimeter,
            arduino,
            measurement_parameters,
            setup_parameters,
            selected_pixels,
            file_path_builder,
            parent=None,
        ):
            super(LowLuminanceMeasurement, self).__init__()

            self.uno = arduino
            self.uno.init_serial_connection()
            self.keithley_source = keithley_source
            self.keithley_source.as_voltage_source(
                measurement_parameters["scan_compliance"]
            )
            self.keithley_multimeter = keithley_multimeter

            self.measurement_parameters = measurement_parameters
            self.setup_parameters = setup_parameters
            self.selected_pixels = selected_pixels
            self.file_path_builder = file_path_builder
            self.parent = parent
            self.stop = False

            self.sweep = ChoppedSweep(
                keithley_source,
                keithley_multimeter,
                samples_per_burst=measurement_parameters["samples_per_burst"],
                settle_time=measurement_parameters["settle_time"],
                multimeter_range=measurement_parameters["multimeter_range"],
                target_relative_sigma=measurement_parameters[
                    "target_relative_sigma"
                ],
                baseline_period=measurement_parameters["baseline_period"],
                scan_compliance=measurement_parameters["scan_compliance"],
            )

        def run(self):
            voltages = self.measurement_parameters["voltages"]

            self.sweep.prepare()
            try:
                for pixel in self.selected_pixels:
                    if self.stop:
                        break
                    log_message("Running on Pixel " + str(pixel))
                    self.uno.trigger_relay(pixel)
                    try:
                        df_data = self.sweep.run_sweep(
                            voltages,
                            progress_callback=lambda value: self.update_progress_bar.emit(
                                "value", value
                            ),
                        )
                    finally:
                        self.uno.trigger_relay(pixel)

                    save_data(
                        df_data,
                        self.file_path_builder(pixel),
                        self.setup_parameters,
                        self.measurement_parameters,
                    )

                    self.update_plot.emit(
                        df_data["voltage"].to_list(),
                        df_data["current"].to_list(),
                        df_data["pd_voltage"].to_list(),
                    )
            finally:
                self.sweep.finish()
                self.hide_progress_bar.emit()
                self.reset_start_button.emit(True)


# ---------------------------------------------------------------------------
# 하드웨어 없이 로직 확인용 시뮬레이션
# ---------------------------------------------------------------------------


def _simulate():
    """
    장비 없이 수렴 로직과 통계 처리를 확인한다.
    가짜 DMM 은 (신호 + 20 uV 노이즈 + 앰프 DC 오프셋) 을 돌려준다.
    """
    rng = np.random.default_rng(0)

    class _FakeDMM:
        def __init__(self):
            self.n = 5
            self.offset = 3.1e-3  # 앰프 DC 오프셋
            self.source = None

        def configure_burst(self, sample_count, multimeter_range=0.1, nplc=1):
            self.n = int(sample_count)

        def read_burst(self):
            signal = self.source.pd_signal()
            return signal + self.offset + rng.normal(0, 20e-6, self.n)

        def reset(self):
            pass

    class _FakeSource:
        """1 cd/m2 당 PD 40 uV (README 의 60 mm 재배치 후 가정)."""

        def __init__(self):
            self.voltage = 0.0

        def set_voltage(self, voltage):
            self.voltage = float(voltage)

        def pd_signal(self):
            luminance = 0.0 if self.voltage <= 0 else 10 ** (self.voltage - 3.0)
            return luminance * 40e-6

        def read_current(self):
            return 1e-8 * 10 ** self.voltage

        def activate_output(self):
            pass

        def deactivate_output(self):
            pass

        def set_low_light_mode(self, on):
            pass

    source = _FakeSource()
    dmm = _FakeDMM()
    dmm.source = source

    sweep = ChoppedSweep(
        source,
        dmm,
        settle_time=0.0,  # 시뮬레이션이므로 대기 없음
        target_relative_sigma=0.02,
        max_cycles=20000,
        max_time_per_point=1e9,
        baseline_period=0.0,
    )
    sweep.prepare()
    # 3.0 / 3.3 / 3.7 / 4.0 V = 대략 1 / 2 / 5 / 10 cd/m2
    df = sweep.run_sweep([3.0, 3.3, 3.7, 4.0])
    sweep.finish()

    print()
    print(df[["voltage", "pd_voltage", "pd_voltage_std", "cycles", "converged"]])
    print()
    print("참값 대비 오차 (목표 +-2%):")
    for _, row in df.iterrows():
        source.voltage = row["voltage"]
        truth = source.pd_signal()
        print(
            "  V = %.2f  참값 %8.1f uV  측정 %8.1f uV  오차 %+6.2f %%"
            % (
                row["voltage"],
                truth * 1e6,
                row["pd_voltage"] * 1e6,
                (row["pd_voltage"] / truth - 1) * 100,
            )
        )


if __name__ == "__main__":
    import sys

    if "--simulate" in sys.argv:
        _simulate()
    else:
        print(__doc__)

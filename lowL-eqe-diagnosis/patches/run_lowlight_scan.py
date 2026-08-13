# -*- coding: utf-8 -*-
"""
저휘도 chopped 스캔 단독 실행 스크립트.

GUI 를 건드리지 않고 chopped 스캔을 검증하기 위한 것이다.
OLED-jvl-measurement/src/ 에 넣고 그 폴더에서 실행한다.
장비 초기화는 랩 소프트웨어의 hardware.py 를 그대로 재사용하므로
주소나 릴레이 동작을 다시 맞출 필요가 없다.

전제:
  - hardware.py 에 configure_burst / read_burst / set_low_light_mode 가 추가돼 있을 것
    (hardware.diff 또는 README 의 수동 적용 절차)
  - chopped_measurement.py 가 같은 src/ 폴더에 있을 것
  - 측정 GUI 는 완전히 종료돼 있을 것 (장비를 이 스크립트가 잡아야 함)

사용 예:

  # 1) SMU 정착 시간 실측 (전압 스텝 후 PD 가 언제 안정되는지)
  python run_lowlight_scan.py --measure-settle --pixel 1 --voltage 3.5

  # 2) 저휘도 chopped 스캔
  python run_lowlight_scan.py --pixel 1 --voltages 2.6,2.8,3.0,3.3,3.6,4.0

  # 3) 드리프트 상쇄 패턴 비교
  python run_lowlight_scan.py --pixel 1 --voltages 3.0 --pattern ABBA

주의: 이 스크립트는 소자에 전압을 인가한다. 어떤 경로로 끝나든
      (정상 종료 / Ctrl+C / 예외) 출력을 끄고 릴레이를 전부 내린다.
"""

import argparse
import datetime as dt
import os
import sys
import time

import numpy as np

import core_functions as cf
import hardware as hw
import chopped_measurement as chop


def parse_arguments():
    parser = argparse.ArgumentParser(description="저휘도 chopped 스캔")
    parser.add_argument("--pixel", type=int, default=1, help="측정할 픽셀 번호 (1-8)")
    parser.add_argument(
        "--voltages",
        type=str,
        default="2.6,2.8,3.0,3.3,3.6,4.0",
        help="쉼표로 구분한 전압 리스트 (V)",
    )
    parser.add_argument(
        "--target", type=float, default=0.02, help="목표 상대 불확실도 (0.02 = +-2%%)"
    )
    parser.add_argument("--settle", type=float, default=0.02, help="전압 스텝 후 대기 (s)")
    parser.add_argument("--samples", type=int, default=1, help="버스트당 샘플 수")
    parser.add_argument("--range", type=float, default=10, help="DMM 레인지 (V)")
    parser.add_argument("--pattern", type=str, default="AB", choices=["AB", "ABBA"])
    parser.add_argument(
        "--compliance", type=float, default=1.0, help="전류 컴플라이언스 (mA)"
    )
    parser.add_argument(
        "--max-time", type=float, default=240.0, help="포인트당 최대 측정 시간 (s)"
    )
    parser.add_argument("--output", type=str, default=None, help="저장 파일 경로")
    parser.add_argument(
        "--measure-settle",
        action="store_true",
        help="스캔 대신 SMU 정착 시간을 실측한다",
    )
    parser.add_argument(
        "--voltage", type=float, default=3.5, help="--measure-settle 에서 쓸 전압 (V)"
    )
    return parser.parse_args()


def connect(settings, compliance_ma):
    """랩 소프트웨어와 같은 방식으로 장비를 연결한다."""
    cf.log_message("Arduino 연결 중...")
    uno = hw.ArduinoUno(settings["arduino_com_address"])

    cf.log_message("Keithley 2450 (SMU) 연결 중...")
    source = hw.KeithleySource(settings["keithley_source_address"], compliance_ma)
    source.as_voltage_source(compliance_ma)

    cf.log_message("Keithley 2100 (DMM) 연결 중...")
    multimeter = hw.KeithleyMultimeter(settings["keithley_multimeter_address"])

    return uno, source, multimeter


def measure_settle_time(source, multimeter, voltage, n_samples=60):
    """
    전압 스텝 직후 PD 전압이 언제 안정되는지 본다.

    settle_time 은 chopped 사이클 시간의 절반 가까이를 차지하므로 실측해서
    줄일 값이다. 출력에서 값이 평평해지는 시점 + 여유를 settle_time 으로 쓴다.
    """
    multimeter.configure_burst(1, multimeter_range=10, nplc=1)

    cf.log_message("0 V 에서 baseline 측정")
    source.set_voltage(str(0))
    time.sleep(0.5)
    baseline = float(np.mean([multimeter.read_burst()[0] for _ in range(10)]))

    cf.log_message("%.2f V 스텝 직후 연속 판독" % voltage)
    source.set_voltage(str(voltage))
    start = time.time()

    rows = []
    for _ in range(n_samples):
        value = multimeter.read_burst()[0]
        rows.append((time.time() - start, value - baseline))

    source.set_voltage(str(0))

    final = float(np.mean([value for elapsed, value in rows[-10:]]))
    print()
    print("  경과(ms)   PD-baseline(uV)   최종값 대비")
    print("  " + "-" * 46)
    for elapsed, value in rows[:25]:
        print(
            "  %8.1f   %14.1f   %9.2f%%"
            % (elapsed * 1e3, value * 1e6, (value / final - 1) * 100 if final else 0)
        )
    print()
    print("최종 안정값: %.1f uV" % (final * 1e6))
    print()
    print("해석: '최종값 대비' 가 목표 정밀도(예: +-2%) 안으로 처음 들어오는 시점이")
    print("      필요한 settle_time 이다. 여유를 조금 주고 --settle 로 쓰면 된다.")
    print("      첫 샘플부터 이미 안정돼 있으면 settle 을 0.005 까지 줄여도 된다.")


def main():
    args = parse_arguments()
    settings = cf.read_global_settings()

    voltages = [float(value) for value in args.voltages.split(",")]
    uno, source, multimeter = connect(settings, args.compliance)

    # 어떤 경로로 끝나든 소자를 끄기 위해 전체를 try/finally 로 감싼다
    try:
        uno.trigger_relay(0)  # 모든 픽셀 OFF 에서 시작
        uno.trigger_relay(args.pixel)
        cf.log_message("픽셀 %d 선택" % args.pixel)

        source.activate_output()

        if args.measure_settle:
            measure_settle_time(source, multimeter, args.voltage)
            return

        sweep = chop.ChoppedSweep(
            source,
            multimeter,
            samples_per_burst=args.samples,
            settle_time=args.settle,
            multimeter_range=args.range,
            target_relative_sigma=args.target,
            baseline_period=0.0,
            chop_pattern=args.pattern,
            max_time_per_point=args.max_time,
            # hardware.py 는 컴플라이언스를 mA 로 받지만 read_current() 는 A 를
            # 돌려주므로, 비교용 값은 A 로 넘긴다
            scan_compliance=args.compliance * 1e-3,
        )

        sweep.prepare()
        try:
            started = time.time()
            df_data = sweep.run_sweep(voltages)
        finally:
            sweep.finish()

        output = args.output
        if output is None:
            output = os.path.join(
                settings["default_saving_path"],
                "lowlight_pixel%d_%s.txt"
                % (args.pixel, dt.datetime.now().strftime("%Y%m%d_%H%M%S")),
            )

        chop.save_data(
            df_data,
            output,
            {"Pixel": args.pixel, "Chop pattern": args.pattern},
            {
                "samples_per_burst": args.samples,
                "settle_time": args.settle,
                "multimeter_range": args.range,
                "target_relative_sigma": args.target,
                "baseline_period": 0.0,
            },
        )

        print()
        print(
            df_data[
                ["voltage", "current", "pd_voltage", "pd_voltage_std", "cycles",
                 "elapsed", "converged"]
            ].to_string(index=False)
        )
        print()
        print("전체 소요: %.1f s" % (time.time() - started))

        failed = df_data[~df_data["converged"]]
        if len(failed):
            print()
            print("수렴하지 못한 포인트가 %d 개 있다 (신호가 노이즈에 묻힘)." % len(failed))
            print("  -> 해당 전압에서는 휘도가 너무 낮다. 더 높은 전압을 쓰거나")
            print("     --target 을 완화하거나 --max-time 을 늘릴 것.")

    finally:
        # 소자 보호: 출력 차단 후 릴레이 전부 내림
        try:
            source.set_voltage(str(0))
            source.deactivate_output()
        finally:
            uno.trigger_relay(0)
            uno.close()
            cf.log_message("출력 차단, 릴레이 전부 OFF")


if __name__ == "__main__":
    sys.exit(main())

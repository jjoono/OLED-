# -*- coding: utf-8 -*-
"""
PD readout 노이즈 원인 규명용 벤치 테스트.

sigma = 20 uV의 범인이 DMM(Keithley 2100)인지 앰프(PDA200C)인지 가린다.
빛이 없는 상태에서 같은 값을 200번 읽고 표준편차를 계산한다. 신호가 0이어야 하는
상황이므로 퍼진 정도가 곧 노이즈다.

실행 전:
  - 측정 소프트웨어를 완전히 종료할 것 (DMM을 이 스크립트가 직접 잡아야 함)
  - 측정 챔버 차광은 평소 측정과 동일하게
  - 측정 중 케이블을 만지거나 움직이지 말 것

테스트 3종:
  (1) DMM 입력에서 앰프 케이블 제거, HI-LO 단락  -> measure(10)  = sigma_1 (DMM 자체)
  (2) 앰프 정상 연결, gain 50 dB, PD 암전         -> measure(10)  = sigma_2 (DMM + 앰프)
  (3) (2)와 같은 상태, 물리 작업 없음              -> measure(0.1) = sigma_3 (100 mV 레인지)

해석 (노이즈는 제곱합: sigma_2^2 = DMM^2 + AMP^2):
  sigma_1 ~= sigma_2  -> DMM 지배. 저휘도 모드에서 100 mV 레인지 전환으로 해결(코드 한 줄).
  sigma_1 <<  sigma_2 -> 앰프 지배. 레인지 전환 무의미. 저잡음 TIA로 교체 검토.
"""

import statistics
import sys

import pyvisa

# 랩 셋업의 Keithley 2100 주소 (OLED-jvl-measurement/usr/global_settings.json 기준)
DMM_ADDRESS = "USB0::0x05E6::0x2100::8011801::INSTR"

rm = pyvisa.ResourceManager()
dmm = rm.open_resource(DMM_ADDRESS)
dmm.timeout = 5000


def measure(range_v, n=200, label=""):
    """range_v 레인지에서 n회 읽고 평균/표준편차를 uV 단위로 출력, sigma를 반환."""
    dmm.write("*RST")
    dmm.write(f"CONF:VOLT:DC {range_v}")
    dmm.write("VOLT:DC:NPLC 1")  # 평소 측정과 동일 조건
    vals = []
    for _ in range(n):
        vals.append(float(dmm.query("READ?")))

    mean = statistics.mean(vals)
    sigma = statistics.stdev(vals)
    ptp = max(vals) - min(vals)
    print(
        f"[{label}] range {range_v} V, n={n}: "
        f"mean {mean * 1e6:9.1f} uV | sigma {sigma * 1e6:7.2f} uV | p-p {ptp * 1e6:8.1f} uV"
    )
    # overload 감지 (2100은 9.9e37을 반환)
    if abs(mean) > 1e30:
        print("  -> OVERLOAD: 앰프 DC 오프셋이 레인지를 넘음. measure(1)로 대체할 것.")
    return sigma


def main():
    print(__doc__)
    print("=" * 78)

    input("[테스트 1] DMM 입력 HI-LO를 단락한 상태인가? 준비되면 Enter...")
    s1 = measure(10, label="1 DMM only, 10V")

    input("\n[테스트 2] 앰프 연결 복구 + PD 암전 상태인가? 준비되면 Enter...")
    s2 = measure(10, label="2 DMM+AMP, 10V")

    print("\n[테스트 3] 같은 상태에서 100 mV 레인지 (물리 작업 없음)")
    s3 = measure(0.1, label="3 DMM+AMP, 100mV")

    print("=" * 78)
    print(f"sigma_1 (DMM only)     = {s1 * 1e6:6.2f} uV")
    print(f"sigma_2 (DMM + AMP)    = {s2 * 1e6:6.2f} uV")
    print(f"sigma_3 (100 mV range) = {s3 * 1e6:6.2f} uV")

    # 앰프 몫을 제곱차로 분리
    amp_sq = s2**2 - s1**2
    if amp_sq > 0:
        amp = amp_sq**0.5
        print(f"\n앰프 기여분 = {amp * 1e6:6.2f} uV,  DMM 기여분 = {s1 * 1e6:6.2f} uV")
        if s1 > amp:
            print("=> DMM 지배: 100 mV 레인지 전환으로 해결 가능 (sigma_3가 개선폭)")
        else:
            print("=> 앰프 지배: 레인지 전환 효과 제한적. 저잡음 TIA 교체 검토")
    else:
        print("\n=> sigma_1 >= sigma_2: DMM이 완전히 지배. 100 mV 레인지 전환 권장")


if __name__ == "__main__":
    try:
        main()
    finally:
        dmm.close()
        sys.exit(0)

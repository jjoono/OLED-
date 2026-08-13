# -*- coding: utf-8 -*-
"""
앰프(PDA200C) 노이즈의 시간 구조 측정 — 벤치 테스트 4.

노이즈 원인은 이미 앰프로 확정됐다 (2026-08 실측: sigma_1 = 1.02 uV DMM 단독,
sigma_2 = 28.62 uV 앰프 연결 -> 앰프 기여분 28.60 uV, DMM 기여분은 분산의 0.13%).
남은 질문은 "그 28.6 uV 가 백색인가, 1/f 드리프트가 섞여 있는가" 다.

이 답이 chopped 스캔의 파라미터를 직접 결정한다:
  - 백색     -> 평균이 1/sqrt(N) 으로 그대로 듣는다. baseline_period 를 길게 잡아
                baseline 공유로 2배 이득을 챙긴다.
  - 1/f 있음 -> 어느 시간 이상 평균해도 sigma 가 안 내려간다(=knee). chopping 주기를
                knee 보다 짧게 가져가야 하고 baseline_period 도 그만큼 짧아야 한다.

셋업: 벤치 테스트 2와 완전히 동일 (앰프 정상 연결, gain 50 dB, PD 암전).
      물리 작업 없음. 데이터만 길게 받는다. 약 2~3분.
      측정 중 케이블을 만지거나 조명을 바꾸지 말 것.

분석 방법:
  받은 시계열을 크기 m 인 블록으로 나눠 블록 평균들의 산포를 본다.
    sigma_avg(m)  : 그냥 m 개를 평균했을 때 남는 노이즈
                    백색이면 sigma_1 / sqrt(m) 를 따라간다. 위로 벗어나면 그게 드리프트다.
    sigma_chop(m) : 인접한 두 블록의 차이를 취했을 때 남는 노이즈 (= chopped 스캔이 하는 일)
                    드리프트는 차분에서 상쇄되므로 sigma_avg 가 꺾여도 이쪽은 계속 내려간다.

  두 곡선이 갈라지기 시작하는 지점이 knee 이고, 그게 chopping 이 실제로 벌어주는 이득이다.

하드웨어 없이 분석 코드만 확인:
    python noise_test2_averaging.py --simulate
"""

import sys
import time

import numpy as np

# 랩 셋업의 Keithley 2100 주소 (noise_test.py 와 동일)
DMM_ADDRESS = "USB0::0x05E6::0x2100::8011801::INSTR"

# 총 샘플 수. NPLC 1 (60 Hz 기준 16.7 ms) + 오버헤드로 샘플당 약 20~30 ms 이므로
# 6000 샘플이면 2~3분. 블록 크기 500 까지 12개 블록이 확보된다.
TOTAL_SAMPLES = 6000
BURST_SIZE = 100  # 한 번의 READ? 로 받는 샘플 수 (시계열 연속성 유지)
BLOCK_SIZES = [1, 2, 5, 10, 20, 50, 100, 200, 500]


def acquire(dmm, total_samples=TOTAL_SAMPLES, burst_size=BURST_SIZE):
    """NPLC 1 로 연속 샘플을 모아 (values, sample_period) 를 반환한다."""
    dmm.write("*RST")
    dmm.write("CONFigure:VOLTage:DC 10")
    dmm.write("SENSe:VOLTage:DC:NPLCycles 1")
    dmm.write("TRIGger:SOURce IMMediate")
    dmm.write("TRIGger:DELay 0")
    dmm.write("SAMPle:COUNt " + str(burst_size))

    values = []
    start = time.time()
    while len(values) < total_samples:
        answer = dmm.query("READ?")
        values.extend(float(v) for v in answer.split(","))
        done = len(values)
        sys.stdout.write("\r  %d / %d 샘플" % (done, total_samples))
        sys.stdout.flush()
    elapsed = time.time() - start
    print()

    values = np.array(values[:total_samples])
    return values, elapsed / len(values)


def analyze(values, sample_period):
    """블록 크기별 sigma_avg / sigma_chop 표를 출력하고 knee 를 찾는다."""
    sigma_1 = float(np.std(values, ddof=1))
    print()
    print("샘플당 시간      : %.1f ms" % (sample_period * 1e3))
    print("전체 시간        : %.1f s" % (sample_period * len(values)))
    print("평균             : %.1f uV" % (float(np.mean(values)) * 1e6))
    print("sigma (1 샘플)   : %.2f uV" % (sigma_1 * 1e6))
    print()
    print("  블록   시간      백색 예측    sigma_avg    sigma_chop    avg/예측")
    print("  " + "-" * 62)

    rows = []
    for m in BLOCK_SIZES:
        n_blocks = len(values) // m
        if n_blocks < 8:
            continue
        blocks = values[: n_blocks * m].reshape(n_blocks, m).mean(axis=1)

        white = sigma_1 / np.sqrt(m)
        sigma_avg = float(np.std(blocks, ddof=1))
        # 인접 블록 차분 -> 드리프트 제거. /sqrt(2) 로 블록 1개 기준 노이즈로 환산
        sigma_chop = float(np.std(np.diff(blocks), ddof=1)) / np.sqrt(2)
        rows.append((m, white, sigma_avg, sigma_chop))

        print(
            "  %5d  %6.2f s   %8.2f uV  %8.2f uV   %8.2f uV     %5.2fx"
            % (m, m * sample_period, white * 1e6, sigma_avg * 1e6,
               sigma_chop * 1e6, sigma_avg / white)
        )

    print()
    print("=" * 70)

    # 1) 입력 개방 의심 -- 이 경우 나머지 해석이 전부 무의미하다
    mean_uv = abs(float(np.mean(values))) * 1e6
    if mean_uv < 1000 and sigma_1 * 1e6 > 100:
        print("경고: 평균이 %.0f uV 로 작은데 sigma 가 %.0f uV 로 크다." % (mean_uv, sigma_1 * 1e6))
        print("  TIA 입력이 개방된 상태일 가능성이 높다. 개방된 입력은 고임피던스")
        print("  안테나로 동작해 주변 픽업을 다 받으므로 이 데이터로는 앰프 노이즈를")
        print("  판정할 수 없다. PD 를 연결한 상태에서 다시 측정할 것.")
        print("=" * 70)
        return sigma_1, None

    # 2) 백색 성분과 드리프트 성분 분리
    #    가장 짧은 시간(블록 1)의 차분 노이즈가 사실상 백색 성분이다.
    sigma_white = rows[0][3]
    drift_ratio = sigma_1 / sigma_white

    # 3) sigma_chop 이 1/sqrt(m) 로 떨어지는가, 평평한가
    #    평평하면 1/f 플리커 (Allan deviation 이 평평한 것이 그 시그니처)
    m_last, _, _, chop_last = rows[-1]
    chop_expected = sigma_white / np.sqrt(m_last)
    flatness = chop_last / chop_expected

    print("백색 성분 (25 ms 차분)     : %.2f uV" % (sigma_white * 1e6))
    print("전체 산포 / 백색 성분      : %.2fx" % drift_ratio)
    print("sigma_chop 평탄도(m=%d)   : %.1fx (1.0 = 백색, 클수록 1/f)" % (m_last, flatness))
    print()

    # 드리프트는 두 가지 방식으로 드러난다. 짧은 시간(25 ms)부터 이미 지배적이면
    # drift_ratio 가 크고, 긴 시간에서만 드러나면 sigma_chop 평탄도가 크다. 둘 중
    # 하나만 걸려도 드리프트 지배로 본다.
    if drift_ratio < 1.3 and flatness < 2.0:
        print("결론: 드리프트 없음 (백색). 측정 구간 전체에서 1/sqrt(N) 이 성립.")
        print("  -> baseline_period 를 넉넉히(5 s 이상) 두고 baseline 공유로 시간을 벌 것.")
        print("  -> samples_per_burst 를 늘리면 그만큼 사이클당 오차가 줄어든다.")
    else:
        if drift_ratio >= 1.3:
            print("결론: 전체 산포의 %.0f%% 가 저주파 드리프트다. 그냥 평균은 무력하다."
                  % ((1 - 1 / drift_ratio) * 100))
        else:
            print("결론: 짧은 시간에서는 백색으로 보이지만 긴 시간에서 드리프트가 드러난다.")
        print("  -> baseline_period = 0 (매 사이클 baseline 재측정). 공유 불가.")
        if flatness > 2.0:
            print("  -> sigma_chop 이 블록 크기에 거의 무관하다 = 1/f 플리커.")
            print("     버스트를 늘려도 사이클당 오차가 안 줄어든다. samples_per_burst = 1 로 두고")
            print("     사이클 수를 최대화할 것 (반복은 1/sqrt(N) 로 듣는다).")
            print("  -> chop_pattern = \"ABBA\" 도 시도해 볼 것 (선형 드리프트 상쇄).")
        else:
            print("  -> sigma_chop 은 1/sqrt(m) 를 대체로 따른다. samples_per_burst 를 늘려도 된다.")
        print()
        print("  chopped 스캔의 사이클당 노이즈 예상치 = %.1f uV (= 백색 x sqrt(2))"
              % (sigma_white * np.sqrt(2) * 1e6))

    print("=" * 70)

    return sigma_1, None


def _simulate():
    """백색 + 1/f 합성 데이터로 분석 코드만 검증한다."""
    rng = np.random.default_rng(1)
    n = TOTAL_SAMPLES
    sample_period = 0.025

    for label, drift_amplitude in [("백색만", 0.0), ("백색 + 1/f 드리프트", 25e-6)]:
        white = rng.normal(0, 28.6e-6, n)
        # 랜덤워크로 저주파 드리프트를 만든다
        drift = np.cumsum(rng.normal(0, drift_amplitude / np.sqrt(n), n))
        values = 8198e-6 + white + drift

        print()
        print("#" * 70)
        print("# 시뮬레이션: " + label)
        print("#" * 70)
        analyze(values, sample_period)


def main():
    import pyvisa

    print(__doc__)
    print("=" * 70)
    input("벤치 테스트 2와 같은 상태(앰프 연결, PD 암전)인가? 준비되면 Enter...")

    rm = pyvisa.ResourceManager()
    dmm = rm.open_resource(DMM_ADDRESS)
    dmm.timeout = 30000
    try:
        values, sample_period = acquire(dmm)
    finally:
        dmm.close()

    np.savetxt("noise_stream.csv", values, delimiter=",")
    print("원본 시계열 저장: noise_stream.csv")
    analyze(values, sample_period)


if __name__ == "__main__":
    if "--simulate" in sys.argv:
        _simulate()
    else:
        main()

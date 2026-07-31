"""
angular_recycling_ideal.py

class B(각도 선별형)의 '이상적 상한' 과 '실현 가능한 DBR' 사이의 간극,
그리고 Al 손실 a 의존성을 계산한다.

핵심 질문: 산란층+각도필터 플랫폼으로 특정 방향에 빛을 얼마나 모을 수 있나?
  - ideal filter : 밴드만 100% 투과, 나머지 100% 반사 (물리적 상한)
  - real DBR     : angular_recycling_bound.py 에서 최적화한 실제 다층막
  - flat         : class A 대표 (벽 = 탈출광 중 33.7%)

결론 구조: ideal 은 손실 a 가 작을 때만 벽을 크게 넘는다. a 가 현실적이면
  재활용 횟수 증가분을 Al 이 먹어버려 이득이 사라진다 -> "손실이 방향성의 실질 한계".
"""
import numpy as np

N_SUB = 1.51
TH_LO, TH_HI = 40., 60.
TH_SUB_LO = np.degrees(np.arcsin(np.sin(np.radians(TH_LO))/N_SUB))
TH_SUB_HI = np.degrees(np.arcsin(np.sin(np.radians(TH_HI))/N_SUB))
TH_C = np.degrees(np.arcsin(1.0/N_SUB))

th = np.linspace(0., 89.5, 180)
w = np.cos(np.radians(th))*np.sin(np.radians(th))
w /= w.sum()
in_band = (th >= TH_SUB_LO) & (th <= TH_SUB_HI)
prop = th < TH_C


def recycle(T, g, a):
    beta_band = float((T*w)[in_band].sum())
    beta_out = float((T*w).sum()) - beta_band
    beta_tot = beta_band + beta_out
    s = (1.0+g)/2.0
    denom_d = 1.0 - (1.0-s)*(1.0-beta_tot)
    Ad_band = (1.0-s)*beta_band/denom_d
    Ad_out = (1.0-s)*beta_out/denom_d
    Au_band = s*(beta_band + (1.0-beta_tot)*Ad_band)
    Au_out = s*(beta_out + (1.0-beta_tot)*Ad_out)
    Au_tot = Au_band + Au_out
    denom = 1.0 - (1.0-Au_tot)*(1.0-a)
    return Au_band/denom, Au_tot/denom


def T_fresnel():
    """평면 계면 (class A 대표)."""
    T = np.zeros_like(th)
    t = np.radians(th[prop])
    c1 = np.cos(t)
    s2 = N_SUB*np.sin(t)
    c2 = np.sqrt(np.clip(1-s2**2, 0, 1))
    rs = ((N_SUB*c1 - c2)/(N_SUB*c1 + c2))**2
    rp = ((N_SUB*c2 - c1)/(N_SUB*c2 + c1))**2
    T[prop] = 1 - 0.5*(rs+rp)
    return T


def T_ideal():
    """이상적 각도필터: 밴드만 투과."""
    T = np.zeros_like(th)
    T[in_band] = 1.0
    return T


if __name__ == '__main__':
    Tf, Ti = T_fresnel(), T_ideal()
    wall = np.sin(np.radians(TH_HI))**2 - np.sin(np.radians(TH_LO))**2
    print(f'[class A 벽] 탈출광 중 밴드 비율 = {wall*100:.1f}%\n')
    print(f'{"a":>6} | {"flat band":>10} {"flat sel":>9} | '
          f'{"ideal band":>11} {"ideal sel":>10} | {"gain":>6}')
    print('-'*62)
    for a in (0.0, 0.01, 0.02, 0.05, 0.10, 0.20):
        fb = max(recycle(Tf, g, a)[0] for g in np.linspace(-0.8, 0.8, 17))
        ft = max(recycle(Tf, g, a)[1] for g in np.linspace(-0.8, 0.8, 17))
        ib = max(recycle(Ti, g, a)[0] for g in np.linspace(-0.8, 0.8, 17))
        it = max(recycle(Ti, g, a)[1] for g in np.linspace(-0.8, 0.8, 17))
        print(f'{a:6.2f} | {fb*100:9.1f}% {fb/ft*100:8.1f}% | '
              f'{ib*100:10.1f}% {ib/it*100:9.1f}% | {ib/fb:5.2f}x')
    print('\n* band = 발생광 대비 목표밴드로 나가는 비율 (phi 창 곱하기 전)')
    print('* sel  = 탈출광 중 밴드 비율 (class A 벽 33.7% 와 비교)')

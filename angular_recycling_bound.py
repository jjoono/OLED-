"""
angular_recycling_bound.py  (vectorized)

외부 광추출층의 "특정 각도 밴드로의 집광" 한계를 재활용까지 포함해 계산.

[이론 구조]
 (1) 재활용 보조정리: 패스 사이에 각도 재랜덤화가 일어나면, 매 패스의 출사 각도
     '모양'이 동일하므로 누적 모양 = 단일 패스 모양. -> 방향성 문제는 "단일 패스에
     어떤 출사 모양을 만들 수 있는가"로 환원된다.
 (2) 클래스 A (에텐듀 보존형: 평면/MLA/freeform/산란층)
     각도를 매핑/랜덤화만 할 뿐 '선별'을 못함 -> 탈출광 중 밴드 비율이
     sin^2(60)-sin^2(40) = 33.7% 로 유계 (평평한 계면이 이미 이를 달성) -> 벽.
 (3) 클래스 B (각도 선별형: DBR/광결정)
     밴드 밖을 '되돌려' 재시도시킴 -> 벽을 넘을 수 있으나 재활용 횟수가 늘어 손실에
     노출. 랜덤화를 '손실 있는 Al'이 아니라 '손실 없는 외부 산란층'에서 시키면
     그 대가를 줄일 수 있다  <-- 본 스크립트가 검증하는 가설.

[모델] 외부영역(산란층) 광자의 마르코프 연쇄. s=(1+g)/2 는 전방산란 확률.
   up  : s -> 상단 구조 도달,        (1-s) -> 기판 복귀
   down: s -> 기판 복귀,             (1-s) -> 상단 구조 도달
   상단 도달 시: 밴드 투과 beta_band / 밴드밖 투과 beta_out / 나머지 반사
   기판 복귀 시: Al 반사에서 (1-a) 로 생존해 up-going 으로 재진입
   * g>0(전방산란)은 초기 추출에 유리하나 재활용에는 불리 -> 상반성 tension.
"""
import numpy as np
from scipy.optimize import differential_evolution

# ---------------- 파라미터 ----------------
N_SUB = 1.51
TH_LO, TH_HI = 40., 60.
LAM = np.linspace(500., 600., 7)
N_HI, N_LO = 2.35, 1.46

TH_SUB_LO = np.degrees(np.arcsin(np.sin(np.radians(TH_LO))/N_SUB))
TH_SUB_HI = np.degrees(np.arcsin(np.sin(np.radians(TH_HI))/N_SUB))
TH_C = np.degrees(np.arcsin(1.0/N_SUB))

th = np.linspace(0., 89.5, 180)
w = np.cos(np.radians(th))*np.sin(np.radians(th))
w /= w.sum()
in_band = (th >= TH_SUB_LO) & (th <= TH_SUB_HI)
prop = th < TH_C                      # 임계각 이내만 공기로 나갈 수 있음


def T_profile(d_list, n_list):
    """각도별 투과율 T(theta) (s,p 평균, LAM 대역 평균). 벡터화 TMM."""
    d_list = np.asarray(d_list, float)
    n_list = np.asarray(n_list, float)
    n_all = np.concatenate(([N_SUB], n_list, [1.0]))
    d_all = np.concatenate(([0.0], d_list, [0.0]))
    L = len(n_all)

    t0 = np.radians(th[prop])                       # (A,)
    kx = N_SUB*np.sin(t0)                           # (A,)
    # cos(theta) in each layer: (A, L)
    cos_t = np.sqrt((n_all[None, :]**2 - kx[:, None]**2).astype(complex))/n_all[None, :]

    Tsum = np.zeros(t0.shape)
    for lam in LAM:
        for pol in ('s', 'p'):
            M11 = np.ones(t0.shape, dtype=complex)
            M12 = np.zeros(t0.shape, dtype=complex)
            M21 = np.zeros(t0.shape, dtype=complex)
            M22 = np.ones(t0.shape, dtype=complex)
            for i in range(L-1):
                n1, n2 = n_all[i], n_all[i+1]
                c1, c2 = cos_t[:, i], cos_t[:, i+1]
                if pol == 's':
                    r = (n1*c1 - n2*c2)/(n1*c1 + n2*c2)
                    tt = 2*n1*c1/(n1*c1 + n2*c2)
                else:
                    r = (n2*c1 - n1*c2)/(n2*c1 + n1*c2)
                    tt = 2*n1*c1/(n2*c1 + n1*c2)
                # M = M @ [[1,r],[r,1]]/t
                a11 = (M11 + M12*r)/tt
                a12 = (M11*r + M12)/tt
                a21 = (M21 + M22*r)/tt
                a22 = (M21*r + M22)/tt
                if i+1 < L-1:
                    dl = 2*np.pi*d_all[i+1]*n_all[i+1]*cos_t[:, i+1]/lam
                    e_m, e_p = np.exp(-1j*dl), np.exp(1j*dl)
                    M11, M12, M21, M22 = a11*e_m, a12*e_p, a21*e_m, a22*e_p
                else:
                    M11, M12, M21, M22 = a11, a12, a21, a22
            t_tot = 1.0/M11
            Tsum += np.real(1.0*cos_t[:, -1]/(N_SUB*cos_t[:, 0]))*np.abs(t_tot)**2
    T = np.zeros_like(th)
    T[prop] = np.clip(Tsum/(2*len(LAM)), 0, 1)
    return T


def recycle(T, g, a):
    """재활용 마르코프 연쇄. 반환 (밴드로 나가는 총비율, 전체 탈출비율)."""
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


def make_stack(x):
    d = np.asarray(x, float)
    n = np.array([N_HI if i % 2 == 0 else N_LO for i in range(len(d))])
    return d, n


def _cost(x, g, a):
    d, n = make_stack(x)
    band, _ = recycle(T_profile(d, n), g, a)
    return -band


def optimize_dbr(npair, g, a, seed=1, maxiter=25):
    bounds = [(30., 200.)]*(2*npair)
    res = differential_evolution(_cost, bounds, args=(g, a), seed=seed,
                                 maxiter=maxiter, popsize=10, tol=1e-4,
                                 polish=True, disp=False)
    d, n = make_stack(res.x)
    return d, T_profile(d, n), -res.fun


if __name__ == '__main__':
    a_loss = 0.10
    print(f'target theta_air [{TH_LO:.0f},{TH_HI:.0f}] -> theta_sub '
          f'[{TH_SUB_LO:.1f},{TH_SUB_HI:.1f}] deg,  critical {TH_C:.1f} deg')
    print(f'lambda {LAM[0]:.0f}-{LAM[-1]:.0f} nm,  Al loss a={a_loss}\n')

    wall = np.sin(np.radians(TH_HI))**2 - np.sin(np.radians(TH_LO))**2
    print(f'[class A 벽] 탈출광 중 밴드 비율 상한 = {wall*100:.1f}%')

    Tf = T_profile([], [])
    print('  --- flat interface (class A 대표) + 재활용 ---')
    for g in (-0.5, 0.0, 0.5):
        band, tot = recycle(Tf, g, a_loss)
        print(f'   g={g:+.1f}: band={band*100:5.1f}%  total_out={tot*100:5.1f}%  '
              f'band/total={band/tot*100:5.1f}%')
    flat_band = max(recycle(Tf, g, a_loss)[0] for g in (-0.5, 0.0, 0.5))

    print('\n[class B] scatterer + DBR:')
    best = None
    for npair in (4, 6):
        for g in (-0.6, -0.3, 0.0, 0.3):
            d, T, band = optimize_dbr(npair, g, a_loss)
            _, tot = recycle(T, g, a_loss)
            print(f'   npair={npair} g={g:+.2f}: band={band*100:5.1f}%  '
                  f'total_out={tot*100:5.1f}%  band/total={band/tot*100:5.1f}%')
            if best is None or band > best[0]:
                best = (band, g, npair, d, T, tot)

    band, g, npair, d, T, tot = best
    print(f'\n>>> best: npair={npair}, g={g:+.2f}, band={band*100:.1f}%'
          f'   (flat 대비 {band/flat_band:.2f}x)')
    print('    d [nm]:', np.round(d, 1))
    np.savez('angular_recycling_result.npz', th=th, T_dbr=T, T_flat=Tf, d=d,
             g=g, a=a_loss, band=band, tot=tot, wall=wall, flat_band=flat_band,
             th_sub_lo=TH_SUB_LO, th_sub_hi=TH_SUB_HI)
    print('\nsaved -> angular_recycling_result.npz')

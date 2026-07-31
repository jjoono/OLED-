"""
angular_recycling_bandwidth.py

핵심 설계 질문: "산란층 + DBR(각도필터)" 플랫폼이 작동하려면 얼마나 단색이어야 하는가?

배경: 이상적 각도필터는 class A 벽(탈출광 중 33.7%)을 크게 넘지만(a=0.1 에서 62%),
  실제 DBR 은 정지대역 가장자리가 파장에 따라 이동하므로 넓은 대역에서는 각도
  선택성이 뭉개진다. 목표 밴드가 기판측 25.2~35.0 deg (폭 ~10 deg) 로 좁아서
  파장폭이 조금만 커져도 필터 기능이 사라진다.

-> 대역폭 dlam 을 스윕하며 최적 DBR 을 설계, 도달 가능한 방향성을 구한다.
   결과는 "이 플랫폼의 실현 조건"을 정량화한다 (Fig 3 후보).
"""
import numpy as np
from scipy.optimize import differential_evolution

N_SUB = 1.51
TH_LO, TH_HI = 40., 60.
N_HI, N_LO = 2.35, 1.46
LAM0 = 550.

TH_SUB_LO = np.degrees(np.arcsin(np.sin(np.radians(TH_LO))/N_SUB))
TH_SUB_HI = np.degrees(np.arcsin(np.sin(np.radians(TH_HI))/N_SUB))
TH_C = np.degrees(np.arcsin(1.0/N_SUB))

th = np.linspace(0., 89.5, 140)
w = np.cos(np.radians(th))*np.sin(np.radians(th))
w /= w.sum()
in_band = (th >= TH_SUB_LO) & (th <= TH_SUB_HI)
prop = th < TH_C


def T_profile(d_list, lams):
    d_list = np.asarray(d_list, float)
    n_list = np.array([N_HI if i % 2 == 0 else N_LO for i in range(len(d_list))])
    n_all = np.concatenate(([N_SUB], n_list, [1.0]))
    d_all = np.concatenate(([0.0], d_list, [0.0]))
    L = len(n_all)
    t0 = np.radians(th[prop])
    kx = N_SUB*np.sin(t0)
    cos_t = np.sqrt((n_all[None, :]**2 - kx[:, None]**2).astype(complex))/n_all[None, :]
    Tsum = np.zeros(t0.shape)
    for lam in lams:
        for pol in ('s', 'p'):
            M11 = np.ones(t0.shape, dtype=complex); M12 = np.zeros(t0.shape, dtype=complex)
            M21 = np.zeros(t0.shape, dtype=complex); M22 = np.ones(t0.shape, dtype=complex)
            for i in range(L-1):
                n1, n2 = n_all[i], n_all[i+1]
                c1, c2 = cos_t[:, i], cos_t[:, i+1]
                if pol == 's':
                    r = (n1*c1 - n2*c2)/(n1*c1 + n2*c2); tt = 2*n1*c1/(n1*c1 + n2*c2)
                else:
                    r = (n2*c1 - n1*c2)/(n2*c1 + n1*c2); tt = 2*n1*c1/(n2*c1 + n1*c2)
                a11 = (M11 + M12*r)/tt; a12 = (M11*r + M12)/tt
                a21 = (M21 + M22*r)/tt; a22 = (M21*r + M22)/tt
                if i+1 < L-1:
                    dl = 2*np.pi*d_all[i+1]*n_all[i+1]*cos_t[:, i+1]/lam
                    e_m, e_p = np.exp(-1j*dl), np.exp(1j*dl)
                    M11, M12, M21, M22 = a11*e_m, a12*e_p, a21*e_m, a22*e_p
                else:
                    M11, M12, M21, M22 = a11, a12, a21, a22
            t_tot = 1.0/M11
            Tsum += np.real(cos_t[:, -1]/(N_SUB*cos_t[:, 0]))*np.abs(t_tot)**2
    T = np.zeros_like(th)
    T[prop] = np.clip(Tsum/(2*len(lams)), 0, 1)
    return T


def recycle(T, g, a):
    bb = float((T*w)[in_band].sum())
    bo = float((T*w).sum()) - bb
    bt = bb + bo
    s = (1.0+g)/2.0
    dd = 1.0 - (1.0-s)*(1.0-bt)
    Adb = (1.0-s)*bb/dd; Ado = (1.0-s)*bo/dd
    Aub = s*(bb + (1.0-bt)*Adb); Auo = s*(bo + (1.0-bt)*Ado)
    Aut = Aub + Auo
    return Aub/(1.0 - (1.0-Aut)*(1.0-a)), Aut/(1.0 - (1.0-Aut)*(1.0-a))


GS = np.linspace(-0.6, 0.6, 7)


def _cost(x, lams, a):
    T = T_profile(x, lams)
    return -max(recycle(T, g, a)[0] for g in GS)


if __name__ == '__main__':
    a = 0.10
    npair = 8
    wall = np.sin(np.radians(TH_HI))**2 - np.sin(np.radians(TH_LO))**2
    print(f'Al loss a={a}, DBR {npair} pairs, lam0={LAM0:.0f} nm')
    print(f'[class A 벽] 탈출광 중 밴드 = {wall*100:.1f}%   '
          f'(flat: band~29%, sel~33.8%)\n')
    print(f'{"dlam(nm)":>9} | {"band":>7} {"sel":>7} | {"vs flat":>8}')
    print('-'*40)
    flat_ref = 0.291     # angular_recycling_ideal.py 결과 (a=0.1)
    rows = []
    for dlam in (0., 10., 30., 60., 100.):
        lams = np.array([LAM0]) if dlam == 0 else np.linspace(LAM0-dlam/2, LAM0+dlam/2, 5)
        res = differential_evolution(_cost, [(30., 200.)]*(2*npair), args=(lams, a),
                                     seed=2, maxiter=30, popsize=10, tol=1e-4,
                                     polish=True, disp=False)
        T = T_profile(res.x, lams)
        band = max(recycle(T, g, a)[0] for g in GS)
        gbest = GS[int(np.argmax([recycle(T, g, a)[0] for g in GS]))]
        tot = recycle(T, gbest, a)[1]
        print(f'{dlam:9.0f} | {band*100:6.1f}% {band/tot*100:6.1f}% | '
              f'{band/flat_ref:7.2f}x')
        rows.append((dlam, band, band/tot, res.x.copy()))
    np.savez('angular_recycling_bandwidth.npz',
             dlam=[r[0] for r in rows], band=[r[1] for r in rows],
             sel=[r[2] for r in rows], wall=wall, a=a, flat_ref=flat_ref)
    print('\nsaved -> angular_recycling_bandwidth.npz')

from omml import *
EQ = []   # (number, title, omml body, Korean note)
def add(n, title, body, note): EQ.append((n, title, body, note))
add('S1', 'BSDF의 방위각 대칭 축약',
    cat(r('BSDF'), par(cat(thi, t(','), sub(r('φ'), r('i')), t(','), sub(r('θ'), r('s')), t(','), sub(r('φ'), r('s')), t(','), r('λ'))), t(' → '), r('B'), par(cat(thi, t(','), sub(r('θ'), r('s'))))),
    '육방 배열 MLA·등방 산란층: 방위각 대칭. θ_s = 0–90°는 투과(BTDF), 90–180°는 반사(BRDF).')
add('S2', '투과 성분과 반사 성분',
    cat(BT(thi), t(' = '), Sum(tht, cat(r('BTDF'), par(cat(thi, t(','), tht)))), t(',      '),
        Sum(thr, BR(thi, thr)), t(' = 1 − '), BT(thi)),
    '두 번째 식은 흡수가 없는 광추출 구조에서 성립(에너지 보존).')
add('S3', '첫 통과의 각도 분포',
    cat(Pk('0'), t(' = '), cat(sub(r('P'), t('sub')), par(th)), t(',      '), Sum(th, cat(sub(r('P'), t('sub')), par(th))), t(' = 1')),
    'P_sub(θ): 쌍극자 모델로 계산한 기판 전달 파워의 각도 분포(η_sub로 정규화). P_sub^(k)는 k번째로 광추출 구조에 도달하는 빛(정규화하지 않음, 합 = 남은 파워).')
add('S4', '통과 사이의 갱신 (광추출 구조에서 반사 → 소자에서 반사)',
    cat(Pk('k+1', thr), t(' = '), RL(thr), Sum(thi, cat(Pk('k', thi), BR(thi, thr)))),
    'R_LED는 되돌아온 각도 θ_r에서 곱한다(편광 평균, 면 방향 균일 가정).')
add('S5', '통과별 추출량과 행렬 급수',
    cat(etak('ext'), t(' = '), Sum(th, cat(Pk(), BT())), t(',      '),
        eta('ext'), t(' = '), Sum(cat(r('k'), t('=0')), etak('ext'), t('∞')), t(' = '),
        sub(r('P'), t('sub')), sup(par(cat(r('I'), t(' − '), sub(r('B'), t('R')), r('R'))), t('−1')), sup(sub(r('B'), t('T')), t('T'))),
    '행벡터 P_sub, B_T; 행렬 B_R(θ_i, θ_r); R = diag[R_LED(θ)]. 소자 흡수 또는 B_T > 0이면 수렴.')
add('S6', 'k번째 통과의 탈출 확률 p_k',
    cat(pk(), t(' = '), frac(Sum(th, cat(Pk(), BT())), Sum(th, Pk()))),
    'k번째로 광추출 구조에 도달한 빛 중 공기로 나가는 분율. p_0은 소자의 P_sub로 정해지고, k ≥ 1은 광추출 구조가 되돌린 분포로 정해진다.')
add('S7', 'k번째 왕복의 손실 A′_k',
    cat(Ak(), t(' = 1 − '), frac(Sum(thr, cat(RL(thr), Sum(thi, cat(Pk('k', thi), BR(thi, thr))))),
                                     Sum(thr, Sum(thi, cat(Pk('k', thi), BR(thi, thr)))))),
    'k번째 통과에서 되돌아온 빛의 각도 분포로 가중한 소자 손실(흡수 + 배면 투과).')
add('S8', '급수의 정확한 통과별 표현',
    cat(eta('ext'), t(' = '), Sum(cat(r('k'), t('=0')), cat(pk(), Prod(cat(r('j'), t('=0')), cat(par(cat(t('1 − '), pk('j'))), par(cat(t('1 − '), Ak('j')))), cat(r('k'), t('−1')))), t('∞'))),
    '(S4)–(S7)에서 바로 나오는 항등식(근사 없음). k = 0 항의 곱은 1.')
add('S9', '각도 무작위화 극한의 p, A′ (본문의 정의)',
    cat(r('p'), t(' = '), frac(Sum(th, cat(BT(), t(' '), cs)), Sum(th, cs)), t(',      '),
        r('A′'), t(' = 1 − '), frac(Sum(th, cat(RL(), t(' '), cs)), Sum(th, cs))),
    'cos θ sin θ: 기판 안 램버시안 분포의 가중치. p는 광추출 구조 고유량(n 1.3/1.5/1.8/2.0 → 0.58/0.43/0.29/0.23).')
add('S10', '식 (3)으로의 축약',
    cat(pk(), t(' = '), r('p'), t(',  '), Ak(), t(' = '), r('A′'), t('   ⇒   '), eta('ext'), t(' = '),
        r('p'), Sum(cat(r('k'), t('=0')), sup(par(cat(par(cat(t('1 − '), r('p'))), par(cat(t('1 − '), r('A′'))))), r('k')), t('∞')), t(' = '),
        frac(r('p'), cat(r('p'), t(' + '), par(cat(t('1 − '), r('p'))), r('A′')))),
    '모든 통과에서 p_k와 A′_k가 같을 때만 정확. 반구 MLA에서는 p_k가 통과마다 줄고 A′_k가 늘어 식 (3)보다 낮다.')
add('S11', '스펙트럼 평균과 식 (2)',
    cat(sub(r('η'), t('EQE')), t(' = '), eta('sub'), eta('ext'), t(',      '), eta('sub'), t(' = '), Sum(r('λ'), cat(r('w'), par(r('λ')), eta('sub'), par(r('λ')))), t(',      '),
        eta('ext'), t(' = '), frac(Sum(r('λ'), cat(r('w'), par(r('λ')), eta('sub'), par(r('λ')), eta('ext'), par(r('λ')))), eta('sub'))),
    'w(λ): 광자 수로 정규화한 발광 스펙트럼. η_int = 1(PLQY = 1, 완전한 전하 균형) 또는 η_int를 η_sub에 포함한 형태.')

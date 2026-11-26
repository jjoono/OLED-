# 반사 혼합과 투과 혼합은 OLED 기판의 광추출에 서로 다르게 작용한다

*(영문 원고 `manuscript_draft.md`의 한국어 대응본. 구성, 수치, 그림·표·식 번호, 참고문헌 번호는 영문본과 같다.)*

**저자** — [미정]

---

## 초록

각도 혼합은 유기발광다이오드(OLED) 기판에 갇힌 빛의 추출을 돕는다고 널리 가정되며, 통계 광선광학의 완전 무작위화 모델이 그 예이다. 이 연구에서는 광학적 상반성을 이용해, 혼합의 이득이 빛을 어디서 섞느냐에 달려 있음을 보인다. 탈출면이 반사로 빛을 섞으면 빛은 주로 탈출할 수 없는 각도로 되돌아가므로 섞인 빛은 덜 탈출하며, 그 감소량은 탈출 확률의 각도 분산이 정한다. 반면 손실 없는 부피를 투과하면서 섞이는 빛에는 이런 감소가 없고, 완전 무작위화 모델이 그대로 성립한다. 세 가지 소자 적층에서 마이크로렌즈, 산란층, 텍스처 등 37개 구조를 평가한 결과, 구조 자체의 완전 혼합은 광추출 효율을 최대 2.2 %p 높이는 데 그친 반면, 전방 산란 부피는 평면 계면 대비 완전 무작위화 이득의 92%에서 99%를 회복한다. 발광체가 다른 두 종류의 마이크로렌즈 소자에서 한 각도 분해 측정은 예측한 각도 재분배를 확인해 주며, 이 결과는 추출 구조에서 각도 혼합을 어디에 둘지에 대한 지침을 준다.

---

## 1. 서론

유기발광다이오드(OLED)는 디스플레이와 조명에 널리 쓰이며, 내부 양자 효율이 1에 가까우므로 이제 효율은 주로 광추출 효율이 제한한다. 평면 OLED는 생성한 빛의 일부만 내보내고, 나머지는 도파 모드와 표면 플라스몬 모드, 전극 흡수, 기판과 공기 사이 계면의 내부 전반사로 잃는다 [1,2]. 내부 전반사로 기판에 갇힌 빛을 꺼내기 위해 기판에 마이크로렌즈 배열(MLA), 산란층, 텍스처 필름을 붙이는 방법이 연구되었다 [3–7]. 그러나 이런 구조로 보고된 외부 양자 효율(EQE)은 대체로 50%와 60% 사이에 머문다 [8,9]. 추출되지 않은 빛이 OLED로 되돌아올 때마다 일부가 흡수되기 때문이다.

이 다중 패스 과정을 설명하는 데에는 두 종류의 모델이 쓰였다. 통계 광선광학은 빛이 매 패스마다 완전히 무작위화된다고 가정한다 [10–12]. 그러면 기판→공기 효율은 $\eta = p/[p+(1-p)A']$이며, 여기서 $p$는 균일 휘도 빛의 탈출 확률, $A'$은 왕복 손실이다 [13]. 각도 분해 모델은 재분배 행렬로 각도 분포를 패스마다 따라간다 [14–17]. 이 모델은 정확하지만 구조를 비교하는 간단한 규칙을 주지 않는다. 두 설명은 실제로 쓰이는 구조에서 가장 크게 어긋난다. 반구 MLA는 추출되지 않은 빛의 대부분을 들어온 각도 근처로 되돌려 보내고, 살아남은 빛은 MLA가 잘 투과시키지 못하는 각도에 모이므로 완전 무작위화 모델은 광추출 효율을 과대평가한다. 이로부터 각도를 더 빨리 섞는 구조를 설계하자는 생각이 나오며, 텍스처 흡수체에서 이런 연구가 있었다 [18]. 그러나 더 빠른 혼합이 실제로 광추출 효율을 높이는지는 밝혀지지 않았다.

이 연구에서는 각도 혼합의 이득이 빛을 섞는 위치에 달려 있고, 광학적 상반성이 이 의존성을 정한다는 것을 보인다. 탈출면이 반사로 빛을 섞으면, 각 각도의 반사율에 비례하여 그 각도로 빛을 되돌려 보내야 한다. 그러면 완전히 섞인 빛은 $p - \sigma_p^2/(1-p)$의 확률로 탈출하며, 여기서 $\sigma_p^2$은 탈출 확률의 각도 분산이다. 반면 손실 없는 부피를 투과하면서 섞이는 빛에는 이런 감소가 없고, 소자에서 섞이는 빛은 소자 반사율이 탈출 확률과 상관되어 있을 때만 감소를 겪는다. 따라서 완전 무작위화 모델은 투과 혼합의 모델이며, 굴절형 추출 구조가 제공하는 반사 혼합을 기술하지 않는다. 또한 혼합 파라미터가 하나인 축약 모델을 유도하고, 세 가지 소자 적층에서 38개 구조를 평가한다. 구조 자체의 완전 혼합은 광추출 효율을 최대 2.2 %p 높이는 데 그치는 반면, 전방 산란 부피는 완전 무작위화 이득의 대부분을 회복한다. 끝으로 모델의 유한 기하 버전을 발광체가 다른 두 종류의 제작 MLA 소자와 비교한다.

---

## 2. 이론과 방법

### 2.1 각도 분해 모델과 상반성

그림 1(a)는 모델을 보여준다. 기판 속 극각을 폭 1°의 구간 $k$로 나누며, 구간의 에텐듀는 $w_k = \sin^2\theta_{k+1}-\sin^2\theta_k$이므로 균일 휘도는 $w_k$에 비례하는 파워에 해당한다. 빛이 구조를 만날 때마다 구간 $k$의 빛 중 $p_k$가 공기로 나가고, 구간 $i$의 빛 중 $B_R(i\to r)$이 구간 $r$로 되돌아온다. 구조는 흡수가 없으므로 $p_i + \sum_r B_R(i\to r) = 1$이다. 반환 행렬을

$$
B_R = M\,\mathrm{diag}(1-\mathbf{p}) \tag{1}
$$

로 쓰면 $M$의 각 열의 합은 1이므로, $M$은 각도 위의 마르코프 연쇄의 전이 행렬이고 $p = \sum_k w_k p_k/\sum_k w_k$는 균일 휘도 빛의 탈출 확률이다.

입사 벡터 $\mathbf{v}^{(0)}$는 발광체가 각 구간에 결합시키는 파워이며, 방위각과 편광에 대해 더하고 합이 1이 되도록 규격화한다. 소자는 되돌아온 빛을 각도 분해 반사율 $R(\theta)$로 반사하고, 따로 언급하지 않으면 각도를 유지한다. $\mathbf{v}^{(j)}$가 구조를 $(j+1)$번째로 만나는 빛이면

$$
\mathbf{v}^{(j+1)} = \mathrm{diag}\,\mathbf{R}\; B_R\,\mathbf{v}^{(j)} , \qquad \eta_{\mathrm{sta}} = \sum_{j\ge 0}\mathbf{p}^{\mathsf T}\mathbf{v}^{(j)} = \mathbf{p}^{\mathsf T}\left(I - \mathrm{diag}\,\mathbf{R}\; B_R\right)^{-1}\mathbf{v}^{(0)} \tag{2}
$$

이고, $\eta_{\mathrm{sta}}$는 기판→공기 광추출 효율, $p_j = \mathbf{p}^{\mathsf T}\mathbf{v}^{(j)}/\sum_k v^{(j)}_k$는 $j$번째 만남의 탈출 확률이다. 스칼라 손실 $A' = 1-\sum_k w_k R_k/\sum_k w_k$는 완전 무작위화 모델과 비교할 때만 쓴다.

수동 상반 구조는 에텐듀에 대한 상세균형

$$
w_i B_R(i\to j) = w_j B_R(j\to i) \tag{3}
$$

을 만족하므로 [19,20], $M$은 대칭 행렬과 닮음이고 고유값은 실수이다. 절댓값이 두 번째로 큰 고유값 $\lambda_2$는 각도 분포가 정상 분포로 이완되는 속도를 정하며, 이를 만남 횟수로 나타낸

$$
N_{\mathrm{mix}} = \frac{1}{1-|\lambda_2|} \tag{4}
$$

를 각도 혼합수라 부른다. 이는 복사 전달의 수송 인자와 가역 마르코프 연쇄의 이완 시간에 대응하는 패스당 양이다 [21,22]. 그림 1(b)–1(d)는 세 구조의 $M$을 보여주며, 평면 계면은 모든 각도를 유지하므로 $N_{\mathrm{mix}} = \infty$이고 램버시안 산란체는 $N_{\mathrm{mix}} = 1$이다. $B_R$의 가장 큰 고유값은 여러 번 순환한 빛을 기술한다. 이 빛은 만남 한 번당 $p_\infty = 1-\rho(B_R)$의 확률로 탈출하며, $\rho$는 스펙트럼 반경이다 [23,24]. 후기 패스 탈출비를 $\xi = p_\infty/\tilde p$로 정의하며, 여기서 $\tilde p$는 균일 휘도 빛에 대한 같은 확률이다. 이 비는 오래 살아남은 빛이 램버시안 빛만큼 쉽게 탈출하면 1이고, 탈출하지도 섞이지도 않는 각도가 있으면 0에 가까워진다. 스펙트럼 양은 Supplement 1의 노트 S1에 설명한 대로 85° 미만의 구간에서 계산한다.

**그림 1.** 각도 위의 마르코프 연쇄로 본 다중 패스 추출. (a) 모델의 개략도. (b)–(d) (b) 평면 계면, (c) 반구 MLA, (d) 산란층의 반환 행렬 $M$. (e) $j$번째 만남에서 반구 MLA에 닿는 빛의 각도 분포. (f) $j$번째 만남의 탈출 확률을 $p$로 나눈 값. 빈 기호는 $p_{\mathrm{mix}}/p$, 점선은 완전 무작위화 모델이다. *(`fig1_mixing_concept.png`)*

### 2.2 축약 모델과 탈출면에서의 완전 혼합

축약 모델은 구조의 $p(\theta)$를 유지하고 $M$을 파라미터가 하나인 족으로 바꾼다. 되돌아온 빛 중 $1-m$은 각도를 유지하고, $m$은 고정된 분포 $q$로 다음과 같이 퍼진다.

$$
B_R(\theta_i\to\theta_r) = \left[1-p(\theta_i)\right]\left[(1-m)\,\delta_{ir} + m\,q(\theta_r)\right]. \tag{5}
$$

식 (3)은 $q \propto w(1-p)$일 때만 성립하며, 이때 어떤 각도로 퍼지는 빛은 그 각도가 되돌려 보내는 빛에 비례한다. 같은 형태가 상반적인 반사율 모델에 나오는 반면 [25], 기존의 다중 패스 추출 모델은 되돌아온 빛을 등방으로 퍼뜨리며 [11,15], 이는 $p$가 각도에 따라 다르면 식 (3)을 어긴다. 식 (5)의 핵은 고유값 1과 $1-m$을 가지므로 $m = 1-\lambda_2$로 둔다.

$m = 1$이면 구조는 되돌아온 빛을 완전히 섞고, 이 빛은 다음 만남에서

$$
p_{\mathrm{mix}} = \sum_r q_r p_r = p - \frac{\sigma_p^2}{1-p}, \qquad \sigma_p^2 = \frac{\sum_k w_k (p_k - p)^2}{\sum_k w_k} \tag{6}
$$

의 확률로 탈출하며, 다중 패스 급수는 닫힌 형태

$$
\eta_{\mathrm{mix}} = \mathbf{p}\cdot\mathbf{v}^{(0)} + \frac{(\mathbf{p}\cdot\mathbf{u})\,[(\mathbf{1}-\mathbf{p})\cdot\mathbf{v}^{(0)}]}{1-(\mathbf{1}-\mathbf{p})\cdot\mathbf{u}}, \qquad u_k = R_k\,q_k \tag{7}
$$

를 가진다. 식 (7)은 $p(\theta)$가 상수일 때만 완전 무작위화 모델과 같다. 분산 항은 굴절면에서 큰데, 굴절면이 임계 원뿔 안에서는 투과시키고 밖에서는 모두 반사하기 때문이다. 예를 들어 평면 계면은 $p$ = 0.40이지만 $p_{\mathrm{mix}}$ = 0.06이다. 완전 무작위화 모델은 섞인 빛에 확률 $p$를 주므로 혼합의 이득을 과대평가한다. 식 (7)은 계수 1의 상반 혼합에 대한 기준값이며 상한은 아니라는 점에 유의해야 한다.

### 2.3 투과에 의한 혼합과 소자에서의 혼합

식 (6)의 감소는 탈출면 자신이 빛을 섞기 때문에 생긴다. 상반성은 반사에 의한 재분배에 반사 요소의 반환 확률을 가중치로 주며, 여기서는 그것이 $1-p$이다. 표 1에 정리한 것처럼 다른 곳에서 섞으면 가중치가 달라진다. 반사가 없는 무손실 층이 투과로 빛을 섞으면, 빛은 패스마다 $w$에 비례하여 퍼진다. 그러면 섞인 빛은 $p$의 확률로 탈출하고, 굴절률이 정합된 기판 안에서 층의 위치와 무관하게 완전 무작위화 모델이 정확히 재현된다.

소자도 반사로 빛을 섞을 수 있으며, 확산 반사 전극이 그 예이다. 이 경우 소자 반사는 열의 합이 $R_i$이고 $w_i D(i\to j) = w_j D(j\to i)$를 만족하는 핵 $D$가 된다. 두 조건을 만족하는 계수 1 핵은

$$
D(i\to j) = \frac{R_i\,R_j\,w_j}{\sum_k w_k R_k} \tag{8}
$$

이다. 이 핵은 빛을 $wR$에 비례하여 퍼뜨리고, 섞인 빛은 다음 만남에서

$$
p_D = \frac{\sum_k w_k R_k\,p_k}{\sum_k w_k R_k} = p + \frac{\mathrm{Cov}_w(R,p)}{\langle R\rangle_w} \tag{9}
$$

의 확률로 탈출한다. 공분산과 평균은 에텐듀 가중이다. 분산 항은 나타나지 않으며, $R$이 각도에 무관하면 $p_D = p$이다. 그러나 강한 마이크로캐비티처럼 소자가 탈출 원뿔 안에서 가장 강하게 흡수하면 공분산이 음이 되어 소자에서의 혼합도 감소를 겪는다(노트 S11).

**표 1.** 혼합 위치 세 가지에서 완전히 섞인 빛의 탈출 확률과, 그에 따른 평면 계면의 $\eta_{\mathrm{sta}}$(Al 소자, 550 nm)

| 빛을 섞는 위치 | 빛이 퍼지는 비례 | 섞인 빛의 탈출 확률 | 평면 계면의 $\eta_{\mathrm{sta}}$ |
|---|---|---|---|
| 탈출면, 반사 | $w(1-p)$ | $p - \sigma_p^2/(1-p)$ | 0.488 |
| 무손실 부피, 투과 | $w$ | $p$ | 0.798 |
| 소자, 반사 | $wR$ | $p + \mathrm{Cov}_w(R,p)/\langle R\rangle_w$ | 0.770 |

### 2.4 광학 시뮬레이션

쌍극자 발광과 적층 반사율은 복굴절 층으로 확장한 Chance–Prock–Silbey 모델과 전달 행렬법으로 계산했고, 측정한 광학 상수를 쓰며 기판은 비간섭으로 다루었다 [26–28]. 기준 적층은 Al (100 nm) / B3PyMPM (61.3 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / TCTA (10 nm) / TAPC (97.5 nm) / ITO (150 nm) / 굴절률 $n$ = 1.51의 유리이다. 두 번째 적층은 Al을 Ag로 바꾼 것이고, 저손실 적층은 참고문헌 [13]을 따랐다. 550 nm에서 세 적층의 $A'$은 각각 0.167, 0.044, 0.026이었다. 결과는 550 nm에서 제시하거나 가시광 영역에 걸쳐 발광 스펙트럼으로 가중한다.

라이브러리는 기판과 굴절률이 정합된 38개 구조로 이루어진다. 평면 계면 하나, 흡수 없는 Henyey–Greenstein(HG) 산란층 6종 [29], 반구 또는 자유형 렌즈의 육각 MLA 9종, 렌즈마다 프로파일, 높이, 또는 위치를 따로 준 무작위 MLA 17종, 가우스 물결면 5종이다(노트 S2). MLA와 텍스처는 높이장 몬테카를로 광선추적기로 계산했고, 반구 MLA의 탈출 확률은 LightTools 결과와 일치했다. 반지름 10 µm의 렌즈에는 광선광학이 적합하다 [30]. HG 산란층은 몬테카를로로 추적하고 평면 출사면과 adding method로 결합했다 [31]. 구간 폭, 광선 수, 추적한 행렬의 잔여 비대칭, 편광은 $\eta_{\mathrm{sta}}$를 최대 0.4 %p 바꾸었다.

### 2.5 소자 제작과 측정

화소가 반지름 1 mm의 원인 주황과 녹색 하부 발광 OLED를 크기 25.4 mm, 두께 0.7 mm의 유리 기판에 제작했다. 주황 적층은 ITO (150 nm) / TAPC (50 nm) / TCTA (10 nm) / TCTA:B3PyMPM:Ir(dmppy-ph)$_2$tmd (25 nm) / B3PyMPM (45 nm) / LiF (1 nm) / Al (100 nm)이고, 녹색 소자는 같은 호스트에 Ir(ppy)$_2$acac를 썼다. MLA 소자에는 같은 두께의 두 번째 유리를 굴절률 정합 오일로 붙였다. 이 유리에는 화소 중심의 3 × 3 mm 패치에 반구 또는 자유형 렌즈의 육각 배열이 있으며, 기판 측면은 검게 칠했다.

각도별 세기 $I(\theta)$는 포토다이오드 고니오미터로 0°에서 80°까지 10° 간격으로 측정했다. 전체 EQE는 램버시안 발광체의 정면 EQE에 각도 인자 $F = 2\int_0^{\pi/2}[I(\theta)/I(0)]\sin\theta\,d\theta$를 곱하여 구했고, 적분은 측정점을 지나는 스플라인으로 했다.

소자는 연쇄를 가로 위치까지 확장하여 모델링했다(노트 S13). 이 모델에서 빛은 두께 $t$인 기판의 앞면을 두 번 만나는 사이에 옆으로 $2t\tan\theta$만큼 이동하고, 후면 반사율은 소자 배치를 따르며, 측면에 닿은 빛은 잃는다. 캐비티와 발광체별 내부 효율 하나는 bare 소자에서 보정했다. 그다음 MLA 소자는 추가 파라미터 없이 예측했고, 각도 인자는 측정값과 같은 방식으로 적분했다. 모델 값은 렌즈 격자 방향과 고니오미터 기하에 대한 범위로 제시한다.

---

## 3. 결과와 논의

### 3.1 반사 혼합과 투과 혼합

그림 1(e)는 Al 소자 550 nm에서 연이은 만남마다 반구 MLA에 닿는 빛의 각도 분포를 보여준다. 첫 만남에서는 빛의 28%가 60° 밖에 있지만 열 번 되돌아온 뒤에는 이 비율이 83%로 늘고, 그림 1(f)에서 보듯이 만남당 탈출 확률은 $p$의 3분의 1로 떨어진다. 즉 패스가 반복되면 빛은 탈출할 수 없는 곳에 모인다. 한두 번의 만남으로 섞는 구조도 $p$로 돌아가지 않으며, HG 산란층과 가파른 물결면은 식 (6)이 예측한 대로 각자의 $p_{\mathrm{mix}}$에 머물렀다.

그림 2(a)는 모든 구조의 $p_{\mathrm{mix}}$를 $p$와 비교한다. 비율 $p_{\mathrm{mix}}/p$는 MLA에서 0.81에서 0.93 사이이고 평면 계면에서는 0.14로 떨어지며, 이는 내부 전반사에 의존하는 면에서 감소가 가장 크다는 것을 뜻한다. 이 연구의 세 적층에서는 소자 반사율이 $p(\theta)$와 거의 상관되어 있지 않으므로 소자에서의 혼합이 이 감소를 피한다. 그래서 그림 2(b)에서 보듯이, 탈출면 혼합에 대한 소자 혼합의 우위는 $p - p_{\mathrm{mix}}$와 함께 커진다. 평면 계면에서 완전 혼합은 $\eta_{\mathrm{sta}}$를 탈출면에서는 0.384에서 0.488로 높이는 데 그쳤고, 소자에서는 0.770, 투과에서는 0.798이었다. 반면 MLA에서 이 우위는 수 %p이다.

그림 2(e)는 각 구조가 자신의 $p(\theta)$를 유지한 채 빛을 완전히 섞을 때의 이득을 보여준다. 평면이 아닌 37개 구조에서 식 (7)은 Al 소자 550 nm에서 $\eta_{\mathrm{sta}}$를 −0.4에서 +1.8 %p 바꾸었고, 스펙트럼 가중을 하면 세 소자 어느 것에서도 최대 2.2 %p 바꾸었다. 반면 완전 무작위화 모델은 같은 구조에서 최대 19 %p의 이득을 예측한다. 반사 혼합은 캐비티가 구조를 잘 투과하지 못하는 각도로 주로 방출할 때 더 도움이 된다. 반구 MLA의 첫 패스 탈출 확률이 $p$보다 훨씬 낮은 500 nm에서는 완전 혼합이 6에서 9 %p를 더하며, 무작위 MLA와 물결면이 이 이득의 대부분을 실현한다(표 S10).

그림 2(d)는 투과 혼합에 이상적인 무작위화기가 필요하지 않음을 보여준다. 평면 출사면 아래에 둔 비대칭 파라미터 $g$ ≥ 0.8, 축소 광학 두께 $\tau^*$ ≈ 0.5의 전방 산란 HG 부피는 Al 소자에서 $\eta_{\mathrm{sta}}$ = 0.764에 이르렀다. 이 값은 평면 계면에서 이상적 투과 혼합 값까지 가는 길의 92%이고, 같은 부피는 Ag와 저손실 소자에서 각각 97%와 99%를 회복했다. 그러나 더 두꺼운 HG 부피는 OLED 쪽으로 더 많은 빛을 반사하므로 이 이점을 잃는다.

이 구분은 OLED에만 해당하지 않는다. 그림 2(c)는 내부 등방 광원과 후면 반사율 0.90을 가진 슬랩의 광추출 효율을 보여준다. $n$ = 1.5에서 평면 출사면의 완전 혼합은 효율을 0.44에서 0.60으로 높인 반면, 이 경우 투과 혼합과 같은 후면 반사체에서의 혼합은 0.87을 주었다. 이 순서는 $n$ = 3.5까지 유지된다. 텍스처 출사면은 $\sigma_p^2$이 작아서 혼합 없이도 완전 무작위화 결과에 가깝다. 같은 순서가 상반 문제인 광 포획에서도 나타난다. 약흡수 슬랩에서 전면이 반사하는 빛만 섞는 경우는 램버시안 후면 반사체보다 경로가 길어진 적이 없었다(표 S28, S39).

**그림 2.** 빛을 섞는 위치에 따른 광추출 효율. (a) 모든 구조의 $p_{\mathrm{mix}}$ 대 $p$. (b) 소자에서의 완전 혼합과 탈출면에서의 완전 혼합 사이의 $\eta_{\mathrm{sta}}$ 차이. (c) 후면 반사체가 있는 슬랩의 광추출 효율 대 굴절률. (d) 평면 출사면 아래 HG 부피의 $\eta_{\mathrm{sta}}$. (e) 구조 자체의 완전 혼합에 의한 $\eta_{\mathrm{sta}}$ 이득(채운 기호)과 완전 무작위화 모델의 이득(빈 기호). (a), (b), (d), (e)는 Al 소자 550 nm에 대한 것이다. *(`fig2_complete_mixing.png`)*

### 3.2 축약 모델과 구조 지도

그림 3(a)는 축약 모델을 식 (2)의 전체 급수와 342개 경우에서 비교하며, 이는 38개 구조, 세 소자, 세 파장으로 이루어진다. $m = 1-\lambda_2$와 상반적 퍼짐을 쓰면 평균 절대 오차는 0.30 %p, 최대 오차는 3.0 %p였고, 그림 3(b)에서 보듯이 특정 소자나 파장에 몰려 있지 않았다. 비교하면 완전 무작위화 모델은 평면 계면에서 최대 60 %p, 등방 퍼짐을 쓴 축약 모델은 얕은 물결면에서 최대 10.5 %p 틀렸다. 또한 극각과 방위각을 모두 따라간 직접 광선추적은 전체 급수와 0.6 %p 이내로 일치했다(노트 S9).

그림 3(c)는 축약 모델의 한계를 보여주며, 모델을 확정한 뒤 고른 21개 추가 구조에서 시험했다. 모델은 연쇄가 가약일 때, 즉 어떤 각도 집합이 반사에 대해 닫혀 있어 $\lambda_2$가 0.999를 넘을 때 실패한다. 옆면 각도가 45°인 원뿔, 피라미드, 프리즘이 이 경우이며, 수직에 가까운 빛이 역반사되고 오차는 약 20 %p였다. 피라미드와 프리즘에서는 급수를 계산하기 전에 $M$에서 이 실패를 예측했다. $\lambda_2$를 처음 되돌아온 빛으로 가중한 이완 속도 $\lambda_{\mathrm{eff}}$와 비교하는 두 번째 진단은 옆면 각도 60°의 원뿔을 찾아냈지만, 사전 시험에서는 신뢰할 만하지 않았다. 나머지 면 구조에서 오차는 최대 7.5 %p였던 반면, 전체 행렬에는 이런 제한이 없다.

그림 4(a)는 구조들을 $(p, N_{\mathrm{mix}})$ 평면에 나타낸 것이며, 구조 종류마다 서로 다른 영역을 차지한다. MLA는 탈출 확률이 0.42에서 0.43으로 가장 높고 $N_{\mathrm{mix}}$는 3에서 8인 반면, HG 산란층은 약 두 번의 만남 안에 이완되지만 탈출 확률이 낮다. MLA의 렌즈 프로파일을 무작위화하면 그림 4(b)에서 보듯이 $p$를 낮추지 않고 $N_{\mathrm{mix}}$가 7.0에서 1.9로 줄었고, 그에 따라 그림 4(c)에서 보듯이 $\xi$는 0.39에서 0.80으로 올랐다. 그러나 더 빠른 혼합의 이득이 함께 생기는 $p(\theta)$의 변화로 상쇄되어 $\eta_{\mathrm{sta}}$는 1 %p 미만으로 변했다. 따라서 효율은 주로 $p$가 정하며, 그림 4(a)의 거의 수직인 등고선이 이를 보여준다.

그림 4(a)는 MLA의 렌즈 형상 최적화가 작은 이득만 주는 이유도 설명한다. 반구 MLA는 이미 균일 휘도의 에텐듀 한계 $1/n^2$ = 0.439에 가까운 탈출 확률을 가진다. 캐비티의 고르지 않은 방출은 첫 패스에서 최대 0.540까지 허용하지만 [32], 굴절 렌즈는 나가는 방향 하나하나를 넓은 범위의 내부 각도로 채우므로 어느 경우에나 약 0.43을 투과시킨다. 이와 일관되게, 참고문헌 [8]의 소자 플랫폼에서 한 LightTools 벤치마크에서 최적화한 자유형 렌즈는 최적화한 반구를 전체 EQE에서 1.3% 넘는 데 그쳤다(노트 S12). 상한에 다가가려면 각도 선택적 투과를 가진 표면이 필요하다.

**그림 3.** 상반적 축약 모델의 정확도. (a) 342개 경우에서 축약 모델의 $\eta_{\mathrm{sta}}$ 대 전체 급수. (b) (a)의 오차를 소자(LL: 저손실)와 파장별로 나타낸 것이며, 막대는 평균 절대 오차이다. (c) 라이브러리 구조와 추가 구조의 최대 오차 대 $\lambda_2 - \lambda_{\mathrm{eff}}$이며, 십자는 가약 연쇄를 표시한다. *(`fig3_reduced_model.png`)*

**그림 4.** 구조의 탈출 확률과 혼합수. (a) 평면이 아닌 모든 구조의 $p$와 $N_{\mathrm{mix}}$, 그리고 반구 MLA의 $p(\theta)$ 형태에 대한 축약 모델의 $\eta_{\mathrm{sta}}$ 등고선(Al 소자, 550 nm). (b) 규격화한 무작위화 파라미터에 대한 $N_{\mathrm{mix}}$. (c) 후기 패스 탈출비 $\xi$ 대 $N_{\mathrm{mix}}$. *(`fig4_structure_map.png`)*

### 3.3 제작 소자

그림 5(a)는 제작 소자의 배치를 보여준다. MLA 패치의 폭은 3 mm에 불과한 반면, 렌즈 유리를 붙인 기판의 두께는 1.4 mm이다. 그 결과 MLA가 되돌려 보낸 빛의 대부분은 패치 밖에 떨어지고, MLA는 주로 첫 패스의 빛에 작용한다. 따라서 이 실험은 모델이 이 빛을 각도에서 어떻게 재분배하는지를 시험한다.

그림 5(b)–5(d)와 표 2는 측정한 각도 거동을 모델과 비교한다. 모델은 두 발광체의 bare, 반구, 자유형 소자를 아무 조정 없이 예측했다. 반구 MLA는 방출을 수직 방향으로 모으고 각도 인자를 약 0.7로 낮춘 반면, 자유형 MLA는 첫 패스의 빛을 퍼뜨려 각도 인자를 1 근처로 유지했다. 측정 평균은 모두 예측 범위 안에 있으며, 주황 자유형 소자의 평균만 범위보다 0.03 높다. 평균 각도 프로파일은 예측 띠의 중심에서 rms로 최대 0.04 벗어난다.

**표 2.** 제작 소자의 각도 인자. 모델 값은 측정 기하에 대한 범위이고, 측정값은 평균 ± 표준편차이다

| 발광체 | 소자 | 모델 | 측정 | 표본 |
|---|---|---|---|---|
| 주황 | Bare | 0.99–1.02 | 1.01 ± 0.05 | 소자 5개 |
| 주황 | 반구 MLA | 0.62–0.68 | 0.66 ± 0.04 | 소자 5개 |
| 주황 | 자유형 MLA | 0.92–1.04 | 1.07 ± 0.07 | 소자 4개 |
| 녹색 | Bare | 1.00–1.02 | 1.02 | 소자 2개 |
| 녹색 | 반구 MLA | 0.65–0.71 | 0.71 ± 0.02 | 소자 1개, 스캔 6회 |
| 녹색 | 자유형 MLA | 0.83–0.91 | 0.87 ± 0.04 | 소자 3개 |

그림 5(e)는 MLA 소자의 정면 EQE와 전체 EQE를 bare 소자와 비교한다. 반구 MLA는 방출을 좁히므로 정면 이득이 광추출 효율을 과대평가하며, 정면 EQE는 bare 소자의 1.2배였지만 전체 EQE는 0.81에서 0.85배에 그쳤다. MLA 소자의 측정 EQE는 무손실 예측의 84%에서 91%였다. 오일로 붙인 렌즈 유리를 지날 때마다 각도에 무관한 투과율 0.925 하나를 두고 주황 반구 소자에 맞추면, 나머지 세 종류 소자에서도 이 차이가 재현되었다.

모델은 MLA 패치에 필요한 폭도 알려 준다. 재순환하는 빛을 추출하려면 패치가 화소 밖으로 기판 두께의 몇 배만큼 뻗어야 하기 때문이다. 10 mm 패치는 20 mm 패치가 주는 이득의 80%를 주는 반면 3 mm 패치는 전체 EQE를 낮추며, 이는 실험실 크기 시료와 화소형 디스플레이에 대한 기존 관찰과 일치한다 [8,17]. 따라서 정면 값과 함께 각도 적분 EQE를 보고해야 한다.

**그림 5.** 제작 소자와 유한 기하 모델. (a) 화소와 MLA 패치가 있는 기판 배치. (b), (c) (b) 주황과 (c) 녹색 소자의 규격화한 각도별 세기이며, 띠는 모델, 기호는 측정이다. (d) 각도 인자. (e) MLA 소자와 bare 소자의 정면 EQE 비와 전체 EQE 비이며, 점선 테두리는 무손실 모델, 채운 상자는 맞춘 렌즈 유리 투과율을 넣은 모델이다. *(`fig5_experiment.png`)*

### 3.4 설계에 대한 시사점

이 결과에서 추출 구조의 설계에 대한 다음 지침이 나온다. 구조를 혼합에 맞추어 최적화하기 전에 식 (6)으로 혼합이 줄 수 있는 이득을 추정할 수 있으며, 이 연구의 MLA에서는 최대 약 2 %p이다. 더 큰 이득을 얻으려면 상반성에 의한 감소가 없는 곳, 즉 전방 산란 부피나 반사율이 $p(\theta)$와 상관되지 않은 소자에서 섞어야 한다. 후자는 Al 소자의 MLA에서 1.5에서 4.1 %p를 더한다. 또한 $N_{\mathrm{mix}}$와 $\xi$는 주로 첫 패스 분포가 $p(\theta)$와 맞지 않을 때 중요하므로, 구조는 먼저 $p$로 순위를 매겨야 한다. 무작위 MLA와 가파른 물결면은 높은 $p$와 빠른 혼합을 함께 가지지만, 주기 MLA에 대한 스펙트럼 가중 우위는 최대 2 %p이다. 이 구조들이 여전히 유용한 것은 주로 각도에 따른 색 변화를 줄이기 때문이다 [33]. 이 결론은 스펙트럼 가중을 해도, 출사 매질의 굴절률이 3.5까지 높아져도, 기판이 유한해도, 전극이 부분 확산이어도 유지된다(표 S40).

---

## 4. 결론

이 연구에서는 광학적 상반성이 각도 혼합의 가치를 혼합의 위치에 묶는다는 것을 보였다. 탈출면이 반사로 섞은 빛은 $p - \sigma_p^2/(1-p)$의 확률로 탈출하는 반면, 손실 없는 부피를 투과하며 섞인 빛은 통계 광선광학의 완전 무작위화 모델이 가정하는 대로 $p$의 확률로 탈출한다. 산란층, 마이크로렌즈, 텍스처 37개 구조에서 구조 자체의 완전 혼합은 스펙트럼 가중 광추출 효율을 최대 2.2 %p 높인 반면, 전방 산란 부피는 평면 계면 대비 완전 무작위화 이득의 92%에서 99%를 회복했다. 파라미터가 하나인 상반적 모델은 곡면 구조와 무작위 구조에서 전체 계산을 평균 절대 오차 0.30 %p로 재현했고, 발광체가 다른 두 종류의 제작 소자는 예측한 첫 패스 빛의 각도 재분배를 확인해 주었다. 이 분석은 기판과 굴절률이 정합되고 광선광학이 성립할 만큼 큰 수동 외부 구조에 적용된다는 점에 유의해야 한다. 이 결과는 추출 구조의 설계에 도움이 될 것으로 기대한다. 광추출 효율을 크게 높이려면 투과 또는 소자에서의 혼합, 기판 두께보다 훨씬 넓은 MLA 패치, 또는 각도 선택 구조가 필요하다.

---

## 연구비

[미정]

## 감사의 글

[미정]

## 이해 상충

저자들은 이해 상충이 없음을 밝힌다.

## 데이터 가용성

이 논문의 결과를 뒷받침하는 데이터는 현재 공개되어 있지 않으나 합당한 요청이 있으면 저자에게서 받을 수 있다.

## 보충 문서

보충 내용은 Supplement 1을 참조한다.

---

## 참고문헌

1. W. Brütting, J. Frischeisen, T. D. Schmidt, et al., "Device efficiency of organic light-emitting diodes: progress by improved light outcoupling," Phys. Status Solidi A **210**, 44–65 (2013).

2. H. Greiner, "Light extraction from organic light emitting diode substrates: simulation and experiment," Jpn. J. Appl. Phys. **46**, 4125–4137 (2007).

3. S. Möller and S. R. Forrest, "Improved light out-coupling in organic light emitting diodes employing ordered microlens arrays," J. Appl. Phys. **91**, 3324–3327 (2002).

4. Y. Qu, J. Kim, C. Coburn, et al., "Efficient, nonintrusive outcoupling in organic light emitting devices using embedded microlens arrays," ACS Photonics **5**, 2453–2458 (2018).

5. J. Song, K.-H. Kim, E. Kim, et al., "Lensfree OLEDs with over 50% external quantum efficiency via external scattering and horizontally oriented emitters," Nat. Commun. **9**, 3207 (2018).

6. W. H. Koo, S. M. Jeong, F. Araoka, et al., "Light extraction from organic light-emitting diodes enhanced by spontaneously formed buckles," Nat. Photonics **4**, 222–226 (2010).

7. M. Kovačič, "Design of external microtextures for efficient light outcoupling in OLEDs with different preferential orientation of emission dipoles," Opto-Electron. Rev. **30**, e141542 (2022).

8. M. Kim, J. Kim, and S. Yoo, "Near-planar light outcoupling structures with finite lateral dimensions for ultra-efficient and optical crosstalk-free OLED displays," Nat. Commun. **16**, 11606 (2025).

9. J. Kim, E. Kim, J. Park, et al., "Toward near-foldable surface light sources with ultimate efficiency: ultrathin substrates embedded with micron-scale inverted lens arrays," ACS Photonics **10**, 1775–1782 (2023).

10. E. Yablonovitch, "Statistical ray optics," J. Opt. Soc. Am. **72**, 899–907 (1982).

11. I. Schnitzer, E. Yablonovitch, C. Caneau, et al., "30% external quantum efficiency from surface textured, thin-film light-emitting diodes," Appl. Phys. Lett. **63**, 2174–2176 (1993).

12. M. A. Green, "Lambertian light trapping in textured solar cells and light-emitting diodes: analytical solutions," Prog. Photovoltaics **10**, 235–241 (2002).

13. J. Kim, S. Lee, S. Yoo, et al., "Minimizing round-trip loss for near-unity quantum efficiency in organic light-emitting diodes," arXiv, arXiv:2610.12065 (2026).

14. J. J. Shiang and A. R. Duggal, "Application of radiative transport theory to light extraction from organic light emitting diodes," J. Appl. Phys. **95**, 2880–2888 (2004).

15. K. Neyts and A. Ullan Nieto, "Importance of scattering and absorption for the outcoupling efficiency in organic light-emitting devices," J. Opt. Soc. Am. A **23**, 1201 (2006).

16. J. Eisenlohr, N. Tucher, O. Höhn, et al., "Matrix formalism for light propagation and absorption in thick textured optical sheets," Opt. Express **23**, A502–A518 (2015).

17. M. Kovačič, D. Samigullina, F. Bouchard, et al., "Analysis and optimization of light outcoupling in OLEDs with external hierarchical textures," Opt. Express **29**, 23701 (2021).

18. H. W. Deckman, C. B. Roxlo, and E. Yablonovitch, "Maximum statistical increase of optical absorption in textured semiconductor films," Opt. Lett. **8**, 491–493 (1983).

19. F. J. J. Clarke and D. J. Parry, "Helmholtz reciprocity: its validity and application to reflectometry," Light. Res. Technol. **17**, 1–11 (1985).

20. E. Veach, "Robust Monte Carlo methods for light transport simulation," Ph.D. thesis, Stanford University (1997).

21. M. C. W. van Rossum and T. M. Nieuwenhuizen, "Multiple scattering of classical waves: microscopy, mesoscopy, and diffusion," Rev. Mod. Phys. **71**, 313–371 (1999).

22. D. A. Levin and Y. Peres, *Markov Chains and Mixing Times*, 2nd ed. (American Mathematical Society, 2017).

23. J. N. Darroch and E. Seneta, "On quasi-stationary distributions in absorbing discrete-time finite Markov chains," J. Appl. Probab. **2**, 88–100 (1965).

24. E. G. Altmann, J. S. E. Portela, and T. Tél, "Leaking chaotic systems," Rev. Mod. Phys. **85**, 869–918 (2013).

25. M. Ashikhmin and P. Shirley, "An anisotropic Phong BRDF model," J. Graph. Tools **5**, 25–32 (2000).

26. R. R. Chance, A. Prock, and R. Silbey, "Molecular fluorescence and energy transfer near interfaces," Adv. Chem. Phys. **37**, 1–65 (1978).

27. C.-K. Moon, S.-Y. Kim, J.-H. Lee, et al., "Luminescence from oriented emitting dipoles in a birefringent medium," Opt. Express **23**, A279 (2015).

28. S. Mladenovski, S. Hofmann, S. Reineke, et al., "Integrated optical model for organic light-emitting devices," J. Appl. Phys. **109**, 083114 (2011).

29. L. G. Henyey and J. L. Greenstein, "Diffuse radiation in the galaxy," Astrophys. J. **93**, 70–83 (1941).

30. K. Tang, R. A. Dimenna, and R. O. Buckius, "Regions of validity of the geometric optics approximation for angular scattering from very rough surfaces," Int. J. Heat Mass Transfer **40**, 49–59 (1997).

31. H. C. van de Hulst, *Multiple Light Scattering: Tables, Formulas, and Applications* (Academic Press, 1980).

32. J. Muschaweck and C. Wiesmann, "Extracting light out of LEDs," Adv. Opt. Technol. **2**, 291–300 (2013).

33. J. Lim, S. S. Oh, D. Y. Kim, et al., "Enhanced out-coupling factor of microcavity organic light-emitting devices with irregular microlens array," Opt. Express **14**, 6564 (2006).

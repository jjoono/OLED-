# 확장 OLED에서 자유형 마이크로렌즈 어레이의 실용적 포화와 후속 광학 설계 경로

**영문 제목 후보**
- *Practical Saturation of Freeform Microlens Arrays for Extended OLED Emitters*
- *When Freeform Microlens Arrays Stop Helping OLEDs: A Design Route Map beyond Hemispheres*
- *Hemispherical Saturation in Freeform OLED Outcoupling*

**저자** — [TBD]

> **본 문서는 영문 원고(`manuscript_draft.md`)의 한국어 대응본입니다.**
> 내용·구조·수식이 동일하며, 국문 발표 및 내부 검토용입니다.

---

## 초록

자유형 마이크로렌즈 어레이(MLA)와 수치 역설계는 OLED의 광추출 효율과 배광을 동시에 제어할 수단으로 널리 기대되어 왔다. 그러나 이 기대는 점광원을 전제한 조명 광학에서 비롯되며, 확장 광원 위에 공면적으로 타일링된 MLA에 그대로 적용되지 않는다. 실험적으로 검증된 dipole-microcavity 광원 모델과 3차원 광선추적을 결합하여, 동일 제약 아래 반구형 렌즈와 축대칭 자유형 렌즈를 비교하고, 총 외부양자효율(EQE)과 네 개의 polar band를 독립 목적함수로 최적화하였다. 최선의 freeform 설계는 동일하게 최적화된 반구를 어떤 목적에서도 7% 넘게 — 총 EQE에서는 1.3% — 앞서지 못하며, 밑돌지도 않는다. 반구 최적해에서 재출발시킨 탐색은 잔여 band 이득이 0.4% 미만임을 확정한다. 표본 전체가 하나의 효율–band power 궤적으로 붕괴하고, 고효율 설계의 배광은 거의 불변이다: 볼록·inverted(오목)·무작위 조립 배열에 걸쳐, 형상 자유도·array 구성·무질서 어느 것으로도 지휘할 수 없는 좁은 창 안에 머문다. 이 실용적 포화를 공면적 타일링에서의 면적 지렛대 부재, 기판 매개 횡방향 혼합, 수동 외부층의 radiance 한계로 설명하고, 그 너머에서 생산적인 설계 경로 — source/cavity engineering, 개구 확장, 각도 선택적 재활용 — 를 제시한다.

---

## 1. 서론

OLED의 외부 광추출은 여전히 소자 효율과 고휘도 구동 수명을 제한하는 핵심 문제다. 유기층 및 기판 도파 모드, 금속 전극과 연관된 손실 모드, 기판/공기 계면의 전반사로 인해 평면 OLED에서 생성된 광자의 상당 부분은 공기 중으로 나오지 못한다 [1,2]. 외부 MLA는 전기적 소자 구조를 건드리지 않으면서 기판 모드를 추출할 수 있어, 대면적 OLED와 유연 OLED에 특히 매력적인 해결책으로 사용되어 왔다 [3–5]. 반구형 또는 준반구형 MLA는 이미 높은 효율을 달성할 수 있으며, embedded hemispherical MLA를 사용한 OLED에서는 매우 높은 EQE도 보고되었다 [4,5]. 나아가 렌즈가 전혀 없는 외부 산란층과 수평 배향 발광체의 조합만으로도 50%를 넘는 EQE가 실증된 바 있다 [13]. 즉 **높은 총 추출 자체는 정교한 렌즈 형상 없이도 도달 가능한 목표**이며, 이 사실은 자유형 형상에 남겨진 고유한 역할이 무엇인지를 되묻게 한다.

그럼에도 더 복잡한 형상에 대한 기대는 계속된다. 자유형 illumination optics와 inverse design은 원하는 배광을 생성하는 강력한 도구이며, 최근에는 micro-LED와 OLED의 광학 패키징에도 적용되고 있다 [6,7]. 직관적으로는 렌즈의 비대칭성, 곡률 분포, 높이, pitch, 또는 microcavity 조건을 동시에 조절하면 총 추출과 특정 방향으로의 방출을 반구형 MLA보다 더 잘 제어할 수 있어 보인다.

이 논문은 이 직관을 정면으로 시험한다. 우리의 질문은 "더 좋은 freeform 렌즈를 찾을 수 있는가?"가 아니라, **확장 OLED 위의 tiled refractive MLA에서 형상 자유도는 실제로 얼마나 가치 있는가?**이다. 렌즈 없는 산란층이 이미 높은 총 EQE에 도달할 수 있다면 [13], 자유형 MLA에 기대할 수 있는 고유한 부가가치는 총 추출 위에 얹히는 **각도 제어**여야 한다. 본 연구는 바로 그 부가가치가 탐색 가능한 설계공간에서 실제로 실현되는지를 정량적으로 조사한다. 이 질문은 실무적으로 중요하다. 자유형 설계와 정밀 성형은 제조 비용과 검증 부담을 높이지만, 얻는 이득이 hemispherical reference를 거의 넘지 못한다면 설계의 다음 단계는 더 복잡한 형상이 아니라 광원, 개구, 또는 재활용 경로를 바꾸는 것이어야 한다.

우리는 제작 소자와의 비교로 검증된 바 있는 동일한 OLED source model [16,17]과 동일한 제조 가능 기하 제약 아래에서, (i) hemispherical MLA, (ii) 축대칭 freeform MLA를 비교하고, 보조 연구로 (iii) 비대칭 3D freeform MLA를 별도 조건에서 시험한다(2.4절). 나아가 inverted, randomly assembled MLA의 두 가지 추가 제조 가능 외부 필름 family에 대해 결과의 일반성을 시험한다(2.5절). 총 EQE뿐 아니라 polar band와 제한된 azimuthal window를 목적함수로 사용하고, 다중 시작점 최적화와 고정밀 재평가로 optimizer dependence를 분리한다. 그 결과를 radiance/étendue 관점 및 기판 내 lateral mixing과 연결한다. 본 연구의 기여는 새로운 MLA 형상을 제안하는 데 있지 않다. 대신 **freeform MLA의 실용적 포화를 진단하고, 그 뒤에 어떤 광학 레버로 넘어가야 하는지를 정량적으로 제시하는 것**이다. 이때의 포화는 보편 정리(universal theorem)가 아니라, 재현 가능한 벤치마크와 불확도를 갖춘 수치적 관측으로 제시된다.

---

## 2. 결과 및 논의

### 2.1 공정하게 정렬된 benchmark: 형상 복잡도 자체의 효과를 분리하다

그림 1은 비교하는 광학계를 보여 준다. OLED microcavity의 dipole emission은 CPS 계산으로 substrate-side angular–spectral distribution $I_{\mathrm{sub}}(\theta,\lambda)$를 얻고, 이를 3차원 광선추적의 source로 사용하였다. 광원은 반경 $r_{\mathrm{OLED}} = 1$ mm의 확장 면광원이고, 기판 두께는 $d_{\mathrm{sub}} = 1.295$ mm이며, 개별 lenslet(반경 ~10 μm)은 2차원 육각배열로 25×25 mm 전면을 타일링한다. 자유형 프로파일은 13개 설계변수(스플라인 제어점; $x_2$–$x_6$ 단조 제약)로 매개변수화하였다. 모든 렌즈는 동일한 기판 굴절률, lens material, pitch, fill factor, 최대 높이 및 최대 draft-angle 제약을 만족하도록 하였다. cavity 두께는 모든 캠페인에서 렌즈와 함께 최적화하였다. 즉 각 렌즈 클래스를 공통 cavity가 아니라 각자의 최적 cavity에서 비교한다.

이러한 정렬은 필수적이다. freeform이 더 넓은 variable range, 더 높은 lens height, 더 큰 aperture, 또는 더 유리한 cavity만 사용하면 형상 효과와 다른 효과가 섞인다. 본 연구에서는 hemispherical MLA도 동일한 outer radius 및 공정 가능한 높이 범위에서 최적화하고, freeform에는 동일하거나 더 불리한 제약을 부과하지 않았다. 따라서 이후의 성능 차이는 가능한 한 형상 자유도 자체의 순수 효과로 해석할 수 있다. 한 가지 스코프를 명확히 해 둔다. 이 통제 비교(controlled comparison)는 위에 명시한 동일 스택·동일 제약의 **축대칭 클래스**(hemisphere vs 축대칭 freeform)에 대해 성립한다. 2.4절의 비대칭·off-axis 탐색은 서로 다른 소자 조건에서 수행된 **별도의 보조 연구(auxiliary study)**이며, 본 절의 통제 비교의 일부가 아니다.

**Fig. 1 | 플랫폼과 비교군.** (a) OLED–substrate–MLA 광학계 개요: 반경 1 mm 확장 광원($r_{\mathrm{OLED}}$), 두께 1.295 mm·$n$=1.51 유리 기판, 반경 ~10 μm lenslet의 육각배열(전체 25×25 mm). (b) 동일 제약 아래 비교하는 두 렌즈 클래스를 렌즈 반경 단위의 같은 축척으로 그린 것: hemispherical reference(자유변수 3개)와 축대칭 freeform(13개). 빈 원은 자유로운 스플라인 제어점 5개이며 양 끝점은 (0,1)과 (1,0)에 고정된다. 그린 freeform은 40–60° 전용 최적해다. (c) 네 polar band와 band 선택성 $S_j$의 정의(전 방위각에 대해). (d) 최적화 흐름: CPS 광원 → MATLAB COM으로 구동하는 LightTools 광선추적 → surrogate 전역탐색 → patternsearch 정련 → 고광선 재평가. *(`fig1_platform.png`)*

### 2.2 목적함수를 바꾸어도 hemispherical reference를 크게 넘지 못한다

그림 2d,e는 총 EQE와 네 개의 polar angular band(0–20°, 20–40°, 40–60°, 60–80°)를 목적으로 했을 때의 최선 성능을 요약한다. 총 EQE와 네 band 각각은 독립적인 전용 단일목적 최적화로, 그 사이 영역은 가중합 스윕으로 다루었다(구간별 상세는 2.3절). 각 목적 $j$에 대해 다음의 relative gain을 정의하였다.

$$
G_j=\frac{\max\left[\mathrm{EQE}_j\mid\mathrm{freeform}\right]}
{\max\left[\mathrm{EQE}_j\mid\mathrm{hemisphere}\right]}.
$$

여기서 최적값은 독립 초기값, 다중 시작점, 그리고 고정밀 ray 재평가를 거친 뒤의 평균으로 정의한다. hemispherical reference는 동일한 제약·스택·patch·fidelity 아래에서, 프로파일 제어점 5개를 사분원 위에 고정하고 cavity 두께와 렌즈 높이를 자유변수로 열어 최적화하였다. 네 polar objective에서 $G_j$ = 1.003, 1.002, 1.070, 1.036 이고 총 EQE에 대해서는 $G$ = 1.013 이다. 총 EQE의 분자는 그 목적 전용 캠페인에서 온다. 고정된 25 × 25 mm 패치에서 독립 시작점 3개(시작점 간 산포 0.8%)로 0.5539에 도달했고, 반구 재출발 대조는 반구 기준값 0.54679와 같은 세션에서 독립적으로 0.5523에 도달했다. 즉 고정 패치에서의 freeform 상한은 서로 다른 두 절차가 0.3% 차이로 두 번 짚어낸 값이다. 가중합 스윕의 더 높은 총합 0.5556은 경쟁하는 값이 아니라 **정규화가 다른 값**이다. 그 캠페인의 최고 설계를 고정 패치에서 재측정하면 0.5428, 100 mm 패치에서 재측정하면 0.5636이 나온다. 아카이브의 0.5556은 실측한 35 mm 값과 100 mm 값 사이에 놓이며, 이는 그 캠페인이 돌던 시점의 유효 패치가 더 컸다는 해석과 부합한다(2.6절; 보충 표 S7). 따라서 그 캠페인의 절대값은 고정 패치의 반구와 직접 비교하지 않는다. freeform 최적해는 동일 조건에서 최적화된 반구를 어떤 목적에서도 7% 넘게 앞서지 못하고, **밑돌지도 않는다.**

이 값에 도달하기까지 하나의 대조 실험이 필요했고, 첫 판본이 그것을 통과하지 못했으므로 전부 보고한다. 무작위 시드로 돌린 원래의 band별 캠페인에서 0–20°와 20–40° gain 은 1보다 **작게**(0.946, 0.966) 나왔다. 이는 광학에 대한 진술로는 성립할 수 없다. 반구 자체가 13변수 feasible set 의 한 점이므로 예산이 무한하면 정의상 $G_j \ge 1$ 이기 때문이다. 따라서 그것은 탐색에 대한 진술이며, 동시에 이 논문 전체에 대한 가장 강한 반론이기도 하다. 자기가 품고 있는 해조차 되찾지 못하는 탐색이 "더 나은 것은 없다" 고 보고할 때 그 말을 믿을 이유가 없다.

그래서 각 arm 의 탐색을 **그 arm 의 반구 최적해에서 재출발**시켜 이 모호함을 제거했다. 반구 점 자체를 시드에 넣고, 그 근방 섭동점들을 함께 넣고, 반구 점에서 국소 pattern search 를 직접 한 번 더 돌린다(4.3절). 반구 기준값은 이전 캠페인에서 읽어오지 않고 같은 세션에서 최종 정밀도로 다시 측정했으며, 두 값은 모든 arm 에서 0.1% 이내로 일치했다. 결과는 두 해석을 깨끗하게 갈라놓는다. 원래 캠페인이 못 미쳤던 두 band 에서는 결손이 회복된 뒤 사실상 거기서 멈춘다. 반구 대비 잔여 이득은 0–20°에서 +0.32%(반복 3회, $t = 1.7$)이고, 20–40°에서는 저장된 두 설계를 탐색 없이 5회씩 새로 재측정한 전용 확정 실험이 잔여 +0.29%가 통계적으로 실재하되($t = 4.1$) 크기가 0.3% 수준임을 보인다. 반면 실제로 여유가 있는 세 목적에서는 동일한 절차가 그것을 압도적 유의도로 찾아낸다. 40–60°에서 +4.70%($t = 115$), 60–80°에서 +3.64%($t = 30$), 총 EQE 에서 +1.01%($t = 69$)이다. +0.29% 잔여까지 분해하는 절차가 그보다 열 배 큰 이득을 놓칠 리는 없다. 앞의 두 band 의 0.4% 미만 잔여는 탐색의 약함이 아니라 그 구간에서 형상 자유도의 값어치가 그만큼 작음을 재는 것이다(그림 2f).

따라서 1 미만의 gain 은 탐색 인공물이었고, 보정된 값은 1.003 과 1.002 다. 이 두 band 에서 형상 자유도는 반구를 회복시킨 뒤 많아야 0.3%를 더 사 올 뿐이다. 한편 60–80° band 에서는 대조 실행이 원래 캠페인의 자체 최적해까지 넘어섰으므로(0.14075 vs 0.13863), 채택된 $G_4$ 는 이전의 1.020이 아니라 1.036이다. 일부 freeform 은 특정 angular bin 에서 국소적인 gain 을 보였지만, total EQE 감소, 다른 bin 으로의 power 이동, Monte-Carlo uncertainty 를 고려하면 hemispherical MLA 를 대체할 정도의 독립 성능 축을 형성하지 못했다.

이 결과는 "freeform이 정확히 반구형이어야 한다"는 뜻이 아니다. 실제 최적 형상들은 기울기와 곡률 분포에서 서로 달랐고, 낮은 효율 영역에는 다양한 형상이 존재했다. 중요한 점은 성능 상단부에서 이 형상 다양성이 유의미한 추가 효율 또는 band power로 연결되지 않았다는 것이다. 따라서 hemisphere는 단순한 convenience baseline이 아니라, 이 클래스의 **practical near-optimum**으로 작동한다.

### 2.3 Pareto 분석: polar shaping은 독립적인 큰 설계 축이 아니다

총 추출과 40–60° band power를 함께 조사하기 위해 가중 목적함수

$$
J_w=w\hat{\eta}_{\mathrm{ext}}+(1-w)\hat{P}_{40\text{–}60}
$$

를 $w=0,0.25,0.5,0.75,1$에서 surrogateopt + patternsearch로 최적화하고, 별도로 유효 무작위 feasible freeform 표본 $N=150$개를 수집하였다. 그 결과가 그림 2c다. 총 EQE는 표본 전반에서 0.12–0.56에 걸쳐 4.5배 변했지만, 모든 (총 EQE, band EQE) 점은 직선이 $R^2 = 0.968$로 설명하는 하나의 **준선형 궤적으로 붕괴**하였다. 즉 가중치를 어느 극단으로 옮겨도 optimizer는 사실상 같은 상단 영역의 설계를 반환하며, 효율과 40–60° band 방출 사이의 트레이드오프는 관측되지 않았다. 이 붕괴는 이 캠페인의 더 큰 유효 패치(2.6절)에 특유한 것도 아니다. 고정 패치 캠페인들이 독립적으로 이를 재현한다 — 볼록 band별 캠페인에서 $R^2 = 0.94$, 오목·무작위 조립 family 에서 0.94 / 0.84 (그림 5 왼쪽 열).

이를 정량화하기 위해 각 band의 선택성(selectivity)을 $S_j = \mathrm{EQE}_{\mathrm{band},j}/\mathrm{EQE}_{\mathrm{total}}$로 정의한다. Lambertian 기준값은 $S(0\text{–}20^\circ)=0.117$, $S(20\text{–}40^\circ)=0.296$, $S(40\text{–}60^\circ)=0.337$, $S(60\text{–}80^\circ)=0.220$이다.

그러나 해석적 Lambertian 값보다는 시뮬레이션 자체에서 얻은 기준이 더 유용하다. 이 플랫폼은 Lambertian 분할을 그대로 내놓지 않기 때문이다 — 우리가 시험한 모든 family가 그 값에서 수 %p 벗어나, 0–20°는 결핍되고 40–60°는 과잉이 된다(2.5절). 따라서 총 EQE 상위 20개 설계의 선택성을 **자연 배광(natural composition)** — 즉 효율만 밀어붙였을 때 설계가 갖게 되는 배광 — 으로 정의한다. 이 분포는 대단히 좁다. 10–90% 폭이 각각 0.092–0.096, 0.277–0.281, 0.358–0.364, 0.233–0.238로 중앙값(0.094 / 0.278 / 0.361 / 0.236)의 수 % 이내다. 이 20개는 서로 다른 네 개의 단일 band 목적함수로 탐색된 설계들이므로, 이 좁음은 목적함수를 공유해서 생긴 인공물이 아니다. 즉 **총 추출이 높아지고 나면 무엇을 최대화하도록 지시했든 배광은 사실상 고정된다.**

이 내부 기준에 대해 band 전용 최적화는 배광을 실제로 이동시키며, 본 연구는 그 이동을 부정하는 대신 정량화한다. 선택성은 0.119(0–20°), 0.300(20–40°), 0.364(40–60°), 0.283(60–80°)으로, 자연 배광 중앙값 대비 각각 +27%, +8%, +1%, +20%이다. 그러나 같은 설계들이 총 추출을 잃는다. 같은 캠페인 자체의 최고 총 EQE 0.548 대비 각각 3%, 2%, 1%, 11% 낮아, 두 효과가 대체로 상쇄된다. 최고 총 EQE에서의 자연 배광을 기준으로 한 각 전용 최적해의 순 band power는 1.22, 1.05, 1.00, 1.07이다(그림 2a,b). 가장 엄격한 비교 — 즉 *다른* band를 최적화하다 우연히 도달한 해당 band 최고값 — 에 대해서는 각각 6.3%, 0.6%, 4.5%, 5.8%만 앞선다.

두 가지가 뒤따른다. 첫째, 지향성 필름의 자연스러운 목표인 40–60° band가 바로 전용 최적화의 이득이 없는 구간이다(선택성 +1%, 순이득 1.00). 총 추출을 최대화하는 배광이 이미 그 band를 최대화하고 있기 때문이다. 둘째, 가장 큰 이동은 0–20°에서 일어나는데, 이 구간은 자연 배광이 Lambertian 대비 오히려 **결핍**되어 있다(0.094 vs 0.117). 형상 자유도는 새로운 방향을 만들어내기보다, 플랫폼이 덜 채워 놓은 방향을 되메우는 데 가장 효과적이다.

다만 선택성이 엄밀히 고정된 것은 아니며, **체계적 편류(systematic drift)**가 존재한다. 전체 표본(기록된 691회 중 유효한 606회, 나머지는 양의 총 EQE를 반환하지 못한 추적)에 대해 총 EQE와 $S_j$의 상관계수를 계산하면 $R(0\text{–}20^\circ)=+0.60$, $R(20\text{–}40^\circ)=+0.55$, $R(40\text{–}60^\circ)=-0.12$, $R(60\text{–}80^\circ)=-0.57$이다. 해석은 명료하다. 효율이 오르면 배광이 저각 쪽으로 약간 기울지만, band 간 순서는 어떤 설계에서도 뒤집히지 않으며 40–60° band가 항상 최대 비중을 유지한다. 이 편류가 Monte-Carlo 노이즈나 협대역 인공물이 아님을 확인하기 위해, 효율 구간에 걸쳐 층화 추출한 설계 20개를 광선 20배(200,000), 3회 독립 반복, 광대역(450–750 nm)으로 재평가하였다. 이 부분표본은 기준 정밀도에서 $+0.59/+0.67/+0.07/-0.70$, 광선 20배에서 $+0.61/+0.68/+0.04/-0.70$ 을 주어 표집에 대해 안정적이다(4.3절 및 보충 표 S4). 부분표본 값이 위의 전체 표본 값과 크기에서 다소 다른 것은 이 표본이 대표 추출이 아니라 층화 추출이기 때문이며, 부호 패턴과 band 순서는 양쪽이 동일하다. 따라서 이 편류는 실재하는 미세 추세이되, 설계자가 활용할 수 있는 독립적인 조향 자유도에는 미치지 못한다. 이 캠페인의 패치 정규화가 만든 인공물도 아니다. 같은 부호 패턴($+,+,\sim 0,-$)이 모든 고정 패치 캠페인에서 재현되고(그림 5 오른쪽 열), 배광 자체는 패치를 4배 키워도 0.24 %p 이상 움직이지 않는다(보충 표 S7). 이 편류와 위에서 제시한 좁은 산포는 서로 모순되지 않는다. 상관계수는 경계에서 멀리 떨어진 설계까지 포함한 전 효율 구간에서 계산한 값이고, 산포는 편류가 이미 진행을 끝내고 배광이 수렴한 최고 효율 설계들 사이에서 계산한 값이기 때문이다.

핵심 결과는 효율과 특정 polar band가 수학적으로 절대 불변이라는 것이 아니다. passive optics는 위치와 각도 사이의 power를 재배분할 수 있으므로, étendue 보존만으로 angular selectivity의 절대적 한계를 주장할 수는 없다. 본 연구가 제시하는 것은 더 제한적이면서도 실용적인 명제다. 즉, **탐색한 제조 가능 freeform class에서는 재배분의 크기가 작아 설계자가 활용할 수 있는 polar-shaping freedom이 사실상 포화한다.** 이 명제는 엄밀한 보편 정리 대신, 재현 가능한 benchmark와 uncertainty를 갖춘 수치적 결과로 제시한다.

**Fig. 2 | Achievable region, Pareto 붕괴, 반구 기준선, 그리고 재출발 대조.** (a) 각 band 전용 최적해의 선택성(빨간 마커)과 총 EQE 상위 20개 설계의 자연 배광(회색 10–90% 띠), Lambertian 분할(점선) 비교. (b) 같은 결과를 선택성 이득 $S_{\mathrm{win}}/S_{\mathrm{nat}}$ × 총 추출 비 $E_{\mathrm{win}}/E_{\max}$ = 순 band 이득(값 표시)으로 분해한 것. 순이득은 최대 1.22를 넘지 않으며 40–60° band 에서는 정확히 1.00이다. (c) 무작위 유효 설계(회색)와 가중합 탐색(초록)의 (총 EQE, 40–60° band EQE) 산점도, 고정밀 가중 최적해 5개(빨간 마름모) 포함. 총 EQE 0.12–0.56 전 구간에서 모든 점이 하나의 준선형 궤적으로 붕괴하며($R^2 = 0.968$), 효율–지향성 트레이드오프가 없다. (d) 네 band 전용 freeform 프로파일과, 자체 최적 높이를 갖는 hemispherical reference의 비교(렌즈 반경 단위). (e) 각 클래스가 해당 목적함수에서 도달한 최고 EQE 와 상대 이득 $G_j$. 모든 값은 고정 패치 캠페인에서 온다(가중합 스윕의 대형 패치 총합은 제외; 보충 표 S7). freeform 값은 정규화가 일치하는 캠페인들 중 더 나은 쪽이며, 이는 탐색 인공물도 정규화 차이도 설계 클래스의 성질로 보고하지 않기 위함이다. freeform 은 어떤 목적에서도 동일 최적화한 반구를 7% 넘게 앞서지 못하고, 밑돌지도 않는다. (f) 반구 재출발 대조: 각 arm 의 탐색을 그 arm 의 반구 최적해에서 재출발시켰을 때의 반구 대비 이득과, 고정밀 3회 반복 기준 유의도 $t$(문턱 2.13). 실제 여유가 있는 세 목적은 $t = 115$, $t = 30$, $t = 69$ 로 회수되고, 원래 캠페인이 못 미쳤던 두 band 는 세션 내 +0.32%, +0.20% 로 문턱 아래다(20–40° 잔여는 전용 5회 재측정이 실재하되 +0.29% 크기임을 확정; 보충 표 S6). *(`fig2_achievable_region.png`)*

**Fig. 3 | 설계공간 전체의 선택성 지도.** 가중합 캠페인의 유효 평가 606회 전체에 대해 band 선택성 $S_j$를 총 EQE에 대해 그린 것으로, 무작위 유효 설계(회색)와 optimizer가 방문한 설계(색상)를 구분하였다. 검은 선은 최소제곱 적합, 빨간 점선은 기판모드 제외 Lambertian 분할(0.117, 0.296, 0.337, 0.220)이다. 표기한 상관계수($R = +0.60$, $+0.55$, $-0.12$, $-0.57$)는 효율이 오를수록 배광이 저각으로 기우는 체계적 편류를 정량화한다. 각 패널의 세로축 폭은 수 %p에 불과하다. band 순서는 어떤 설계에서도 뒤집히지 않고 40–60°가 항상 최대 비중을 유지한다. 광선 20배·3회 반복·450–750 nm 광대역 수렴검사 값은 보충 표 S4에 있다. *(`fig3_selectivity_map.png`)*

### 2.4 비대칭 및 off-axis freeform은 중요한 stress test다

대칭성 자체가 결과를 제한하는지 확인하기 위해, 3D 비대칭 freeform lenslet과 제한된 $(\theta,\phi)$ window를 직접 최적화하였다. 목적함수는

$$
\mathrm{EQE}_{\mathrm{win}}=\int_{\theta_1}^{\theta_2}\int_{\phi_1}^{\phi_2}
I_{\mathrm{air}}(\theta,\phi)\sin\theta\,d\phi\, d\theta
$$

로 정의하였다. 이 시험은 매우 중요하다. 만약 비대칭 freeform이 같은 total EQE에서 window power 또는 contrast를 크게 높인다면, 앞 절의 포화는 단지 축대칭 parameterization의 한계였을 것이기 때문이다.

스코프에 관한 주의가 필요하다. 이 비대칭 탐색은 2.1절의 통제 비교와 **다른 소자 조건**에서 수행된 별도의 보조 연구다. 구체적으로, 반사 전극이 다른 스택(Ag 전극; 통제 비교의 ITO 기반 스택과 상이), 이방성 발광(anisotropic emitter) 셀, 그리고 52개 설계변수(통제 비교의 13개 대비 훨씬 풍부한 매개변수화)를 사용하였다. 따라서 이 결과는 2.1절의 동일 제약 비교에 산입되지 않으며, 다음의 독립적 보조 증거로 읽어야 한다. **훨씬 더 풍부한 매개변수화와 다른 스택 조건에서도, 의미 있는 각도 조향은 관측되지 않았다.**

우리의 계산에서는 비대칭 형상이 far-field의 중심과 세부 분포를 이동시킬 수는 있었지만, hemispherical MLA보다 크게 높은 **절대 window power**를 안정적으로 얻지는 못했다. 이 비교는 의도적으로 비율이 아닌 정성적 진술로 보고한다. 스택·발광 셀·변수 개수가 통제 benchmark와 모두 다르므로, 여기에 수치 이득을 적으면 위에서 분리해 둔 스코프를 무시한 나란한 비교를 유도하게 된다. 이 보조 연구에서 넘어오는 것은 효과의 **부호**이지 그 크기가 아니다. 이 결과는 "주기 array가 어떤 상황에서도 조향할 수 없다"는 뜻이 아니다. 비대칭 prism, diffractive element, metasurface, 또는 충분한 aperture expansion은 방향성 광분배를 만들 수 있다 [8,9]. 다만 균일한 확장 OLED source 위의 coextensive refractive MLA에서는, 탐색한 범위 내에서 그러한 재배분이 hemispherical reference를 크게 넘는 유용한 output channel로 이어지지 않았다는 것이 본 연구의 범위 내 결론이다.

### 2.5 MLA family 전반의 일반성 (Generality across MLA families)

2.1–2.3절의 benchmark는 하나의 array 구성—기판 출광면에 성형된 볼록 lenslet—에 대한 것이다. 관측된 포화가 이 특정 구성의 산물이 아님을 확인하기 위해, 두 가지 추가 제조 가능 MLA family에 동일한 시험을 적용한다. (i) **Inverted(오목형) MLA**: 동일한 freeform profile class의 오목형 대응물로, 초박형 기판 외면에 딤플을 새겨 실현된 바 있다 [16]. (ii) **Randomly assembled MLA**: 비주기 array를 pseudo-random supercell로 표현한 구성으로, lenslet 위치를 육각 격자에서 jitter시키고 각 lenslet에 동일 class에서 독립적으로 추출한 무작위 profile을 부여한 뒤 supercell을 주기적으로 타일링한다(4.5절). 이 supercell 표현은 무질서의 correlation length가 supercell 크기보다 충분히 작을 때 진짜 무작위 array와 통계적으로 동등하다.

family 공간의 범위를 명시해 둔다. 여기서 비교하는 세 구성은 모두 **외부 광추출 필름**이다. 즉 이미 기판에 들어온 빛에 작용하며, 마지막 광학 계면이 필름/공기 경계다. 전극 아래에 놓이는 sub-electrode MLA — 예컨대 기판에 식각한 고굴절 렌즈 array를 고굴절 spacer로 평탄화하고 그 위에 소자를 제작하여 기판 모드가 아닌 도파 모드를 표적하는 구조 [5] — 는 본 연구의 범위 밖이다. 발광 스택 자체를 바꾸므로 본 연구가 일관되게 사용하는 고정 광원 프로토콜로 비교할 수 없기 때문이다. 따라서 본 연구의 결론은 외부 필름에 대한 진술이며, 광추출 일반에 대한 진술이 아니다.

각 family에는 전체 프로토콜의 경량 버전을 적용한다. 즉 achievable region을 채우는 100개의 유효 무작위 설계, 총 EQE에 대한 1회의 전용 최적화, 그리고 최선 후보의 고정밀 재평가다. 각 family에 대해 앞서 확립한 포화의 세 가지 signature를 검사한다. (a) (총 EQE, band EQE) 점들의 준선형 붕괴, (b) 최고 효율 설계들 사이의 좁은 자연 배광, (c) band 전용 최적화의 순 band power 이득이 수십 % 이내에 머무름.

**Inverted(오목형) family**는 시험이 완료되었고, 세 signature가 모두 재현되되 한 가지 시사적인 차이를 보인다. 총 추출은 볼록 family 대비 6% 낮다(최고 총 EQE 0.517 vs 0.548). 돌출이 아니라 파낸 표면이므로 예상된 결과다. 상위 20개 설계의 자연 배광은 0.113/0.303/0.340/0.215이며 10–90% 폭은 각각 0.102–0.114, 0.288–0.305, 0.338–0.345, 0.213–0.225로 역시 좁다. 즉 이 family에서도 효율이 높아지면 배광은 고정된다. band 전용 최적화는 여기서도 총 추출을 대가로 선택성을 산다. 가장 큰 이동인 60–80° band에서 선택성을 35% 끌어올리는 대신 총 EQE를 12% 잃어 순 band power 이득은 1.19이며, 나머지 band의 순이득은 0.98, 0.99, 1.05다. 즉 최대 조향 예산은 볼록 family(1.20)와 사실상 같고, 다만 그것이 나타나는 band가 다르다.

시사적인 차이는 **자연 배광 자체가 구성마다 조금 이동한다**는 점이다. 동일한 patch 크기·광원·스택 아래에서 볼록 family는 0.094/0.278/0.361/0.236으로 고각 쪽으로 기운 반면, 오목 family는 0.113/0.303/0.340/0.215로 Lambertian partition에 더 가깝다. 이 이동은 실재하지만 작다 — 어느 band에서도 0.021을 넘지 않는다 — 그리고 과대해석하지 않도록 주의할 필요가 있다. 같은 볼록 캠페인도 배광을 '상위 20개의 중앙값' 대신 '전체 표본의 평균'으로 잡으면 0.094/0.278/0.361/0.236에서 0.097/0.278/0.350/0.242로 옮겨간다. 즉 family 간 이동폭은 한 family 안에서 평균 방식을 바꿀 때 생기는 이동폭과 같은 크기다. 두 구성이 하지 *못하는* 것은 질적으로 다른 배광에 도달하는 것이다. 이 논점은 2.7절에서 다시 다룬다.

**무작위 조립 family** 역시 세 signature를 재현하며, 그 결과가 위 결론을 더 날카롭게 만든다. 최고 총 EQE는 0.522로 볼록(0.548)과 오목(0.517) 사이에 놓인다. 다른 family와 동일하게 상위 20개 중앙값으로 잡은 자연 배광은 0.111/0.299/0.347/0.213이고, 고정밀 재평가한 승자 실현 3개는 독립적으로 0.110/0.298/0.349/0.213을 준다. 같은 캠페인의 편향 없는 무작위 표집 구간만의 평균은 0.095/0.281/0.360/0.233이다. 두 통계량은 볼록·오목 값 어느 한쪽을 지목하는 대신 그 사이를 감싼다. 이것이 배광에 남은 여유가 얼마나 좁은지를 가장 분명히 보여 준다. lenslet 위치를 격자에서 흩뜨리고 각 lenslet이 독립 추출된 profile을 갖는데도, 이 family는 여전히 주기적인 두 family와 같은 좁은 창 안에 놓인다. 최선 실현을 서로 다른 무질서 seed 3개로 반복하면 0.5216 ± 0.0011(변동계수 0.2%)이므로, 특정 실현은 무의미하고 조립 통계만이 문제가 된다. 선택성–효율 상관 $R = +0.80, +0.85, -0.11, -0.85$ 는 볼록의 전체 표본 패턴($+0.60, +0.55, -0.12, -0.57$)을 바깥 band 에서 더 큰 크기로 재현하며, 0 에 가까운 40–60° 값은 $-0.11$ 대 $-0.12$ 로 근접한다.

오목 결과와 함께 놓으면 배광이 애초에 얼마나 움직일 수 있는지의 상한이 잡힌다. 세 가지 구성, 두 가지 평균 방식, 네 개의 단일 band 목적함수, 총 506회의 유효 평가에 걸쳐 우리가 얻은 모든 고효율 배광은 0.094–0.113(0–20°), 0.278–0.303(20–40°), 0.340–0.361(40–60°), 0.213–0.242(60–80°) 안에 있다. 위치 무질서, lenslet별 형상 무질서, 클래스 내부의 형상 자유도는 배광을 그대로 둔다. 볼록에서 오목으로의 이산적 구성 변경은 배광을 움직이지만, 같은 창의 가장자리까지일 뿐이며 그 크기는 통계량 선택의 효과와 비슷하다. 배광은 플랫폼(확장 광원, 두꺼운 기판, 공면적 외부 굴절 필름)의 성질이지, 그 안에서 쓸 수 있는 어떤 설계 변수의 함수가 아니다. 한 가지 비대칭은 기록해 둘 만하다. 무작위 family에서 60–80° band power는 총 추출과 사실상 무상관($R = +0.04$)인 반면 그 선택성은 가파르게 감소($R = -0.85$)하므로, 이 구간에서 총 추출의 개선분은 전부 저각 band로 귀속된다. 그림 5는 family당 3-panel 일반성 검사로 이 결과를 담는다. 세 signature가 세 family 모두에서 재현되므로, 포화 진술은 단일 array 구성에서 시험한 세 개의 제조 가능 외부 필름 family로 확장된다.

**Fig. 5 | MLA family 일반성 검사.** 세 family — 볼록 freeform(기준), inverted(오목형), randomly assembled(pseudo-random supercell) — 각각에 대해 3개 패널. 모든 통계량은 각 family 자신의 평가 로그에서 **동일한 방식으로** 계산하여, family 간 차이가 통계량 정의를 바꾼 결과로 보이지 않도록 하였다. *(a1–c1)* family의 전체 표본 설계에 대한 총 EQE 대 40–60° band EQE와 최소제곱 직선. 준선형 붕괴가 세 family 모두에서 재현된다($R^2 = 0.94$, 0.94, 0.84). *(a2–c2)* family의 자연 배광(상위 20개 중앙값, 색상 막대)과 볼록 기준(회색), Lambertian 분할(점선 테두리)의 비교. 빈 원은 같은 family의 전체 표본 평균으로, family 간 산포와 통계량 간 산포가 비슷한 크기임을 보여 준다. 세 family 모두 0.09–0.11 / 0.28–0.30 / 0.34–0.36 / 0.21–0.24 창 안에 놓인다. *(a3–c3)* band별 선택성–효율 상관 $R$과 볼록 기준의 비교. 부호 패턴(+, +, ~0, −)이 모든 family에서 재현되며 무작위 조립 family에서 크기가 가장 크다. *(`fig5_families.png`)*

### 2.6 포화의 물리적 원인: 면적, 혼합, 그리고 radiance envelope

관측된 포화는 세 개의 서로 보완적인 관점에서 설명된다.

첫째, 각도 압축에는 일반적으로 output aperture의 확장이 필요하다. 작은 source와 큰 macro-lens의 조합은 source area보다 큰 output area를 사용해 solid angle을 줄일 수 있다. 반면 large-area OLED를 덮는 tiled MLA에서는 각 lenslet의 aperture와 그 lenslet이 담당하는 source patch가 같은 비율로 증가한다. 따라서 lenslet의 절대 크기를 키우는 것만으로는 point-source collimation의 이득을 얻지 못한다.

둘째, substrate thickness와 source size는 lateral mixing을 결정한다. 기판을 따라 전파한 광은 하나의 lenslet 아래에서 발생한 광만 보지 않으며, 인접 cell에서 온 광과 섞인다. 본 연구의 기하(lenslet 반경 ~10 μm, $d_{\mathrm{sub}}=1.295$ mm)에서는 기판 두께가 lenslet pitch보다 두 자릿수 이상 크므로 이 혼합이 특히 강하다. 이 혼합의 크기는 기하만으로 정해진다. 임계각 광선은 1.295 mm 기판을 한 번 가로지르는 동안 $d_{\mathrm{sub}}\tan\theta_c = 1.16$ mm, 재활용 왕복마다 2.32 mm 횡방향으로 이동하며, 이는 렌즈렛 반경 10 μm의 각각 116배와 232배다. 즉 각 lenslet은 자기 자신보다 두 자릿수 넓은 기판 영역에서 출발한 광을 받는다. 이 결과는 freeform surface가 처리할 수 있는 input phase space가 이미 평균화되어 있음을 보여 준다.

전용 패치 크기 연구가 이 혼합 스케일을 직접 잰다(보충 표 S7). 고효율 설계 하나를 고정하고 텍스처 패치만 15 / 25 / 35 / 100 mm 로 바꿔 최종 정밀도로 재평가하면 총 EQE 는 0.5168 / 0.5428 / 0.5513 / 0.5636 이다. 초기 상승은 감쇠길이 약 9 mm — 임계각 왕복 변위 2.32 mm 의 4배 — 로 포화하는 듯 보이지만 먼 꼬리는 그렇지 않다. 35 mm 와 100 mm 사이에서도 총합이 2.2% 더 오르는데, 이는 앞의 세 점에 맞춘 단일 지수 꼬리가 허용하는 양의 약 3배다. 즉 여러 번 재활용된 광은 탈출 전에 수십 회의 왕복에 걸쳐 횡방향으로 이동하며, 이것이 이 절의 혼합 기구가 장거리에서 작동하는 모습이고, 유한 필름의 총 EQE 를 크기에 대해 느리게 증가하는 함수로 만든다. 전 캠페인에 쓴 25 mm 패치는 100 mm 값보다 3.7% 낮으므로 고정 패치의 절대 EQE 는 하한이다. 이 유한 크기 페널티에는 실측 대응물이 있다. Ref. [17]에서 발광 개구 1 mm 위에 반경 2 mm 짜리 통상 MLA 필름을 얹었을 때 EQE 증강이 사실상 없었는데(bare 35.6% 대 35.4%), 광이 추출되기 전에 필름 가장자리 밖으로 새어 나갔기 때문이다 — 본 연구의 패치 시리즈가 정량화한 것과 같은 손실 채널이다. 반면 배광은 패치 수렴 상태다. 25 mm 에서 100 mm 로 총합이 3.8% 오르는 동안, 같은 4배 확대에서 어느 band 선택성도 0.24 %p 이상 변하지 않는다. 이 논문의 모든 비율·상관·클래스 비교는 동일 패치에서 평가되어 영향을 받지 않는다.

셋째, 주어진 source radiance와 출광 면적에서 passive external layer가 특정 solid angle에 공급할 수 있는 power에는 radiance/étendue envelope가 존재한다 [10]. 본 연구에서는 이를 global impossibility theorem이 아니라, numerical frontier를 해석하는 기준으로 사용한다. 평면 계면 또는 hemisphere가 그 envelope에 가까워질수록, 더 복잡한 profile은 새로운 radiance를 만들기보다 기존 power의 작은 재배분만 수행한다. 2.3절에서 관측된 준선형 붕괴와 Lambertian 근방의 선택성은 이 해석과 정합적이다.

마지막으로 이 시뮬레이션 결과의 인식적 무게를 분명히 해 둘 필요가 있다. 본 논문의 핵심 주장은 부정적 주장—특정 profile의 우월성이 아니라 실용적 포화—이므로, 이상화된 시뮬레이션 조건은 형상 자유도가 이득을 드러내기에 가장 유리한 조건이다. 표면은 수학적으로 완벽하고, 정렬은 정확하며, 제작 형상 오차·정렬 오차·lens-to-lens 산포가 전혀 없다. 실제 제작 결함은 freeform profile이 hemispherical reference에 대해 갖는 이득을 줄일 수만 있을 뿐, 새로 만들 수는 없다. 따라서 시뮬레이션에서 얻은 null result는 보수적이다. 제작 공정은 여기서 보고한 작은 이득을 더 줄일 수는 있어도, 결론을 뒤집을 수는 없다.

### 2.7 MLA 포화 이후의 설계 경로

MLA 형상 최적화의 수익이 작다는 결론은 광추출이 더 이상 개선될 수 없다는 뜻이 아니다. 대신 바꾸어야 하는 물리량이 형상이 아니라는 뜻이다. 그림 4는 세 가지 후속 경로를 정리한다.

**Source/cavity engineering.** Microcavity thickness, dipole orientation, 반사 전극 및 resonant structure는 substrate-side source distribution 자체를 바꾼다. 이는 이미 OLED angular emission control에 널리 쓰이는 방법이며 [12,13], 본 연구에서는 MLA가 포화된 뒤 가장 먼저 검토할 source-side lever로 위치시킨다. 특히 Song 등의 결과 [13]는 외부 구조로는 렌즈 없는 산란층만 사용하면서 발광체 배향(source-side lever)으로 50% 이상의 EQE에 도달하였는데, 이는 본 연구의 포화 논지와 일관된다. 즉 높은 총 추출은 외부 형상의 정교화가 아니라 광원 쪽 레버로 확보되었고, 남는 질문—자유형 MLA가 그 위에 각도 제어를 추가로 제공하는가—에 대한 본 연구의 답은 탐색 범위 내에서 부정적이었다.

**Aperture expansion.** Small emitter 또는 pixel-level architecture에서 output/source area ratio를 실제로 크게 만들 수 있다면, macroextractor 또는 비공면 optical element로 angular compression이 가능하다. 이 경로는 tiled MLA와 동일한 평면적 확장성을 갖지 않지만, collimation이 최우선인 응용에는 적합하다.

본 연구의 데이터는 이 경로를 반대편에서 한정한다. 공면적 필름 기하 **안에서** 우리가 바꾼 모든 것 — 형상 자유도, 볼록 대 오목 구성, 위치 및 lenslet별 형상 무질서 — 이 자연 배광을 움직인 폭은 어느 band에서도 0.02 이하였고, 이는 평균 방식을 바꿀 때 생기는 차이와 같은 크기다(2.5절). 즉 공면적 필름 내부에는 의미 있는 지향성 레버가 없다. 각도 압축이 요구된다면 출력/광원 면적비 자체를 바꿔야 하며, 이는 공면적 기하 안에서 재설계하는 것이 아니라 그 밖으로 나가는 것을 뜻한다.

**Angular-selective recycling.** 목표 angular band 밖의 광을 선택적으로 반사하여 재시도시키는 angular filter는 non-selective refractive MLA가 제공하지 못하는 selectivity를 만들 수 있다. 비선별 굴절층의 단일 통과(single-pass) 대역 선택성은 해석 모델에서 약 33.7%, 수치 계산에서 33.8%로, 산란 비대칭 파라미터 $g$에 무관하게 일치한다. 반면 이상적 filter 모델에서 round-trip loss 10%일 때 40–60° band로의 전달은 62.1%에 이를 수 있는 반면, 평면 class-A reference는 29.1%에 머문다. 실제 8-pair DBR의 전달행렬 계산에서는 단색 source에서 48.2%, 100 nm 대역폭에서 32.5%가 얻어져, source bandwidth와 loss가 이 경로의 현실적 제약임을 보인다. 이 기술은 새 구조 제안이 아니라, 이미 알려진 angle-selective OLED 및 photon-recycling 개념 [8,11,14,15]이 MLA 포화 뒤 왜 필요한지를 보여 주는 control이다.

**Fig. 4 | Angular-selective recycling 경로와 design-route map.** (a) Markov 재활용 모델: 왕복 손실 $a$에 대한 40–60° band 전달률(발생광 대비)을 이상적 각도필터와 평면 비선별 reference에 대해 비교. $a=10\%$에서 각각 62.1%와 29.1%지만 $a$가 커지면 두 곡선이 수렴한다. 즉 실질 한계를 정하는 것은 필터 성능이 아니라 손실이다. (b) 실현 가능한 8-pair 유전체 다층막의 전달행렬 계산: 단색 광원 48.2%에서 100 nm 대역폭 32.5%로 떨어져 단일 통과 벽 아래로 다시 내려간다. (c) 평면 계면·최적화 DBR·이상적 필터의 단일 통과 각도 투과율과 40–60°(공기 기준) 수용창(음영). 이로부터 단일 통과 선택성 벽은 33.7%(해석)/33.8%(수치)이며 산란 비대칭 $g$와 무관하다. (d) Design-route map: 반구 기준선으로 포화를 확인한 뒤, 목표 성능지표별로 손을 뻗어야 할 레버. 세 경로 어디에도 렌즈 형상 최적화는 없다. *(`fig4_recycling_routes.png`)*

---

## 3. 결론

우리는 확장 OLED 위의 제조 가능한 tiled refractive MLA에서 freeform 형상 자유도의 실용적 가치를 체계적으로 시험하였다. 총 EQE와 네 polar band 각각을 독립 단일목적으로 최적화하고, 가중합 스윕과 무작위 설계 150개로 achievable region을 지도화했으며, 별도 조건의 비대칭 off-axis 탐색을 stress test로 두었다. freeform MLA는 동일 조건의 hemispherical MLA를 크게 넘는 안정적인 성능 이득을 보이지 않았다. 모든 탐색을 반구 최적해에서 재출발시킨 대조 실험은 이것이 탐색이 아니라 설계공간의 성질임을 확인해 준다. 문제가 된 두 band 에서 많아야 +0.3% 를 주는 바로 그 절차가, 실제로 여유가 있는 모든 목적에서는 +1.0%에서 +4.7% 를 되찾는다. 모든 설계가 준선형 효율–band power 경계로 붕괴했고, 최고 효율 설계들 사이에서는 어떤 목적함수로 찾았든 배광이 수 % 이내로 고정되었다. 특정 band를 전용 목적으로 삼으면 그 배광이 이동하기는 한다 — 선택성 기준 최대 27% — 그러나 총 추출을 최대 11% 잃어, 다른 목적을 추구하다 얻은 설계 대비 순 band power 이득은 6%를 넘지 않는다. 이 결과는 hemisphere가 모든 광학계에서 최적이라는 보편 명제가 아니다. coextensive tiling, extended source, substrate-mediated lateral mixing, 그리고 주어진 radiance/area 조건이 함께 성립할 때 형상 복잡도만으로 확보할 수 있는 추가 자유도가 작다는 정량적 결론이다. 두 개의 추가 제조 가능 외부 필름 family에 대한 일반성 검사는 이를 '배광이 애초에 얼마나 움직일 수 있는가'에 대한 상한 진술로 날카롭게 만든다. 볼록·오목·무작위 조립의 세 구성과 각각의 두 평균 방식 전체에 걸쳐, 모든 고효율 배광은 0.09–0.11(0–20°), 0.28–0.30(20–40°), 0.34–0.36(40–60°), 0.21–0.24(60–80°) 창 안에 있다. lenslet 위치를 무작위화하고 각 lenslet에 독립 profile을 부여해도 배광은 그 창 안에 머무르며(무질서 실현 간 변동계수 0.2%), 오목으로의 이산적 변경은 배광을 창의 가장자리(0.113/0.303/0.340/0.215)로 옮길 뿐 그 너머로 보내지 못한다. 따라서 배광은 플랫폼의 성질이며, 형상 자유도·array 구성·위치 및 형상 무질서 어느 것으로도 지휘할 수 있는 설계 변수가 아니다. 포화 진술은 시험한 세 개의 제조 가능 외부 필름 family에 대한 것이며 그 밖으로 확장되지 않는다. 전극 아래의 내부 MLA는 구성상 본 범위 밖이다.

따라서 본 연구는 freeform MLA를 포기하라는 결론이 아니라, 설계 결정을 앞당기는 기준을 제공한다. hemisphere benchmark에서 포화가 확인되면, 다음 개선은 더 많은 lens-shape degrees of freedom이 아니라 source/cavity engineering, aperture expansion, 또는 angular-selective recycling에서 찾아야 한다. 이 design-route map은 OLED뿐 아니라 확장 박막 광원을 사용하는 PeLED, QLED, micro-LED optical packaging에도 적용 가능한 실용적 출발점이 된다.

---

## 4. 방법

### 4.1 Trans-scale source model

OLED stack의 dipole emission은 CPS formalism으로 계산하였다. 파장 및 substrate-side angle에 따른 $I_{\mathrm{sub}}(\theta,\lambda)$는 macroscopic ray tracing의 spectral/angular source로 사용하였다. 광원은 반경 $r_{\mathrm{OLED}}=1$ mm의 확장 면광원으로, 두께 $d_{\mathrm{sub}}=1.295$ mm의 기판 하면에 배치된다. CPS input, material dispersion, dipole orientation, internal radiative efficiency, emission spectrum, layer thickness 및 wavelength sampling은 표 S1에 정리한다. 동일한 trans-scale 파이프라인—CPS dipole microcavity source와 LightTools 3차원 광선추적의 결합—은 본 연구 그룹의 선행 논문에서, 본 연구의 플랫폼과 매우 가까운 구성으로 제작 소자와 비교하여 실험적으로 검증된 바 있다. Ref. [17]에서는 발광 개구 반경 1 mm, 기판 두께 1.4 mm, $n_{\mathrm{sub}} = 1.51$ — 사실상 본 연구의 기하 — 의 OLED에 대해, bare·MLA 부착·구조화 outcoupler 부착 소자의 실측 각도별 발광 강도 분포와 전류효율이 세 구성 모두에서 trans-scale 시뮬레이션과 잘 일치했다(ref. [17]의 Fig. 3 및 Table 2; 실측 최대 EQE 각각 35.6%, 35.4%, 48.0%). Ref. [16]에서는 같은 파이프라인에 회절 수준의 BSDF 를 결합하여, inverted MLA 소자의 실측 발광 패턴 — 육각 격자임에도 원통 대칭인 출력 — 과 근사 Lambertian 각도 EL 을 재현했으며, 그 소자의 실측 EQE 는 평면 대조군 30.5% 대비 58.0%였다. 즉 검증은 광원 모델, 미세 텍스처 출사면의 광선광학 처리, 그리고 부착형과 inverted 두 렌즈 구성을 모두 포괄하며, 본 연구와 같은 emitter 플랫폼·같은 기판 기하에서 이루어졌다.

### 4.2 Lens classes와 제조 제약

Hemispherical reference와 축대칭 freeform을 동일 제약 아래 비교하고, 3D 비대칭 freeform은 별도 조건의 보조 연구로 다루었다. 통제 비교의 모든 구조는 동일한 lens material, pitch, fill factor, nominal height, boundary continuity 및 최대 draft angle을 사용하였다. 개별 lenslet(반경 ~10 μm)은 2차원 육각배열로 25×25 mm 전면을 덮는다. 축대칭 freeform surface는 13개 설계변수의 스플라인 제어점으로 매개변수화하였고($x_2$–$x_6$ 단조 제약), LightTools native FreeformEntity로 구현하였다. geometry self-intersection, negative thickness, 제조 불가능한 draft angle을 갖는 후보는 최적화 전에 제거하였다. 보조 연구의 비대칭 freeform은 52개 설계변수를 사용하며, Ag 반사 전극 스택과 이방성 발광 셀 등 통제 비교와 다른 소자 조건에서 수행되었다(2.4절).

### 4.3 Optimization과 검증

목적함수는 total EQE, polar-band EQE, 또는 $(\theta,\phi)$ window EQE로 정의하였다. 각 objective에 대해 surrogate 기반 전역 최적화(surrogateopt)와 patternsearch 정련, 다중 시작점을 사용하고, low-ray search 결과는 independent high-ray repeats로 재평가하였다. achievable region의 편향 없는 표본화를 위해 유효 무작위 feasible 설계 $N=150$개를 수집하고, 가중 목적함수 $J_w$를 $w\in\{0,0.25,0.5,0.75,1\}$에서 최적화하였다. 최적화 예산, ray number, independent run 수, feasible sample 수 및 표준편차는 표 S2에 보고한다. freeform superiority는 single best run이 아니라 고정밀 평균과 hemisphere reference의 차이로 판단하였다.

2.2절의 반구 재출발 대조는 목적함수·제약·기하·정밀도를 band별 캠페인과 동일하게 두고 탐색의 출발점만 바꾼다. 각 arm 에 대해 그 arm 의 반구 최적해에 해당하는 13변수 벡터 — 사분원 위에 고정된 제어점과 그 arm 자신의 최적 cavity 두께·렌즈 높이 — 를 시드에 넣고, 각 변수 범위의 8% 이내에서 뽑은 섭동점 8개를 함께 넣는다. 그 다음 세 갈래로 민다. 시드 집합에서 출발하는 `surrogateopt`, surrogate 승자에서 출발하는 `patternsearch`, 그리고 반구 점에서 직접 출발하는 `patternsearch` 다. 마지막 것이 그 점이 13차원에서 국소 최적인지를 가장 직접적으로 시험한다. 반구 자신을 포함한 모든 후보를 탐색 정밀도로 선별한 뒤 승자를 최종 정밀도로 재평가한다. 반구 기준값은 이전 캠페인에서 가져오지 않고 같은 세션에서 다시 측정하여 스크립트 간 설정 차이가 차이값에 섞이지 않게 했으며, 아카이브 값은 교차 검증에만 쓴다. 개선 판정은 차이가 자유도 $2(N_{\mathrm{rep}}-1) = 4$ 의 pooled 단측 95% $t$ 값, 즉 표준오차의 2.13배를 넘을 때로 한다. 유일한 경계 사례였던 20–40°($t = 2.1$)는 저장된 두 설계를 탐색 없이 각 5회씩 새로 재측정하여 확정했다. 5회 기준 문턱 1.86에 대해 $t = 4.1$, 잔여 크기 +0.29%다. 이 대조는 구조상 한쪽으로만 열려 있음에 유의해야 한다. 반구가 후보에 포함되므로 결과가 그보다 의미 있게 낮아질 수 없고, 승자를 노이즈 섞인 다수의 탐색 평가 중 최댓값으로 뽑으므로 절차 자체가 이득을 보고하는 쪽으로 기울어 있다. 따라서 이 대조에서 나온 null 은 보수적인 결과다.

2.3절의 선택성 편류가 수치 인공물이 아님을 확인하기 위해 별도의 수렴검사를 수행하였다. 광선 수를 기본 대비 20배(200,000)로 늘리고, 3회 독립 반복하며, 협대역 대신 광대역(450–750 nm) 스펙트럼으로, 효율 구간에 걸쳐 층화 추출한 설계 20개를 재평가하였다. 총 EQE와 band 선택성 $S_j$의 상관계수는 세 반복에서 각각 $R(0\text{–}20^\circ)=0.59,\,0.61,\,0.64$; $R(20\text{–}40^\circ)=0.67,\,0.67,\,0.72$; $R(40\text{–}60^\circ)=0.07,\,0.04,\,0.05$; $R(60\text{–}80^\circ)=-0.70,\,-0.71,\,-0.76$으로 재현되어, 관측된 편류가 Monte-Carlo 노이즈나 협대역 인공물이 아닌 실재 추세임을 보였다(보충자료).

### 4.4 Angular-selective recycling model

각도별 transmission/reflection $T(\theta,\lambda)$를 가진 external layer와 round-trip loss $a$를 갖는 reflective electrode를 Markov recycling model로 표현하였다. 이상적 angular filter는 target band에서 $T=1$, 그 밖에서 $T=0$으로 정의했다. 실현 가능한 filter는 alternating high/low-index dielectric multilayer(8-pair DBR)의 transfer-matrix calculation(Python)으로 계산했다. 비선별 굴절층의 단일 통과 대역 선택성은 해석 모델에서 33.7%, 수치 계산에서 33.8%로 산란 비대칭 파라미터 $g$에 무관하였다. 이 모델은 MLA와 DBR를 결합한 소자 제안이 아니라, non-selective refractive MLA의 포화와 selective recycling의 차이를 분리하기 위한 reference calculation이다.

### 4.5 일반성 검사를 위한 MLA family 기하

Inverted MLA는 동일 freeform profile class의 오목형 대응물로, 4.2절의 profile 매개변수화·제약·기각 규칙을 오목 표면에 적용한다 [16]. Randomly assembled MLA는 pseudo-random supercell로 표현한다. 즉 여섯 개의 통계 하이퍼파라미터(충전율, 반경 jitter, 위치 jitter, 평균 aspect, aspect jitter, 반구↔무작위 profile 혼합비)로 무질서를 규정하고 — 단일 렌즈 형상 대신 이 통계량을 최적화한다 —, 각 lenslet에 동일 매개변수화 class에서 독립 추출한 무작위 profile을 부여한 뒤, supercell을 array 전면에 주기적으로 타일링한다. 이 구성은 무질서의 correlation length가 supercell 크기보다 충분히 작을 때 무작위 array와 통계적으로 동등하다 (슈퍼셀당 렌즈렛 6 × 6, 높이 격자 141 × 141 이므로 무질서 correlation length는 렌즈렛 피치 1개, 즉 슈퍼셀의 1/6). 두 family는 각자가 비교되는 대상에 맞추어 서로 다른 프로토콜을 따른다. Inverted family는 볼록 benchmark와 동일한 band별 프로토콜(단일 band arm 4개, arm당 surrogate 60회 + polish 15회, 고정밀 재평가)을 사용한다. 그래야 자연 배광을 볼록과 똑같은 방식으로 생성된 설계 모집단에서 뽑게 된다. Randomly assembled family는 설계변수가 렌즈 프로파일이 아니라 무질서 통계량이므로, 무작위 실현 50개를 표집한 뒤 여섯 개의 조립 통계를 최적화한다. 예산과 ray 수는 표 S2에 보고한다. 무작위 family만 표본 수를 줄였다 — 탐색은 5,000 ray·16 파장, 최종 재평가는 10,000 ray·151 파장 — 슈퍼셀 추적 비용이 단일 렌즈 텍스처보다 훨씬 크기 때문이다. 고정 기하 하나에 대한 보정 결과, 이 축소는 4.3배 빠르면서 band 선택성을 최대 0.39 퍼센트포인트만 이동시켰고 이는 Monte-Carlo 산포와 같은 크기이므로 보고하는 비율에는 영향이 없다(보충 표 S3).

---

## Funding

[TBD — 과제/기관 정보 입력 필요]

## Acknowledgments

[TBD]

## Disclosures

저자들은 이해상충이 없음을 밝힌다.

## Data availability

본 논문의 그림과 표의 기반이 되는 모든 결과 아카이브 — 각 최적화 캠페인의 평가
로그와 최적해(`pareto_front_result.mat`, `opt_4band_result_25by25.mat`,
`freeform_EQEtotal_result.mat`, `opt_hemisphere_result.mat`,
`opt_4band_inverted_result.mat`, `stress_random_result.mat`,
`warmstart_hemisphere_result.mat`, `convergence_check_result.mat`,
`calibrate_random_cost.mat`, `patch_convergence_result.mat`, `patch_convergence_100.mat`, `reeval_confirm_2040_result.mat`) 및 재활용 모델 출력(`angular_recycling_result.npz`,
`angular_recycling_bandwidth.npz`) — 는 이를 생성한 스크립트, 그리고 아카이브에서
모든 그림을 직접 재생성하는 `make_figures.py` 와 함께 합리적 요청 시 교신저자로부터
제공받을 수 있다. 2.4절의 별도 조건 비대칭 탐색 데이터도 합리적 요청 시 제공된다.
CPS 광원 모델과 LightTools 프로젝트 파일의 실행에는 LightTools 라이선스가 필요하다.

---

## 참고문헌

1. Brütting, W.; Frischeisen, J.; Schmidt, T. D.; Scholz, B. J.; Mayr, C. **Device Efficiency of Organic Light-Emitting Diodes: Progress by Improved Light Outcoupling.** *Phys. Status Solidi A* **2013**, *210*, 44–65.

2. Yablonovitch, E. **Statistical Ray Optics.** *J. Opt. Soc. Am.* **1982**, *72*, 899–907.

3. Möller, S.; Forrest, S. R. **Improved Light Out-Coupling in Organic Light Emitting Diodes Employing Ordered Microlens Arrays.** *J. Appl. Phys.* **2002**, *91*, 3324–3327.

4. Wrzesniewski, E.; et al. **Enhancing Light Extraction in Top-Emitting Organic Light-Emitting Devices Using Molded Transparent Polymer Microlens Arrays.** *Small* **2012**, *8*, 2647–2651. https://doi.org/10.1002/smll.201102662.

5. Qu, Y.; Kim, J.; Coburn, C.; Forrest, S. R. **Efficient, Nonintrusive Outcoupling in Organic Light Emitting Devices Using Embedded Microlens Arrays.** *ACS Photonics* **2018**, *5*, 2453–2458. https://doi.org/10.1021/acsphotonics.8b00255.

6. Kim, S.; Shin, J. M.; Lee, J.; Park, C.; Lee, S.; Park, J.; Seo, D.; Park, S.; Park, C. Y.; Jang, M. S. **Inverse Design of Organic Light-Emitting Diode Structure Based on Deep Neural Networks.** *Nanophotonics* **2021**, *10*, 4533–4541. https://doi.org/10.1515/nanoph-2021-0434.

7. Ni, Y.; Feng, D.; Ma, D. **Design of Freeform Microlens Arrays with Prescribed Luminance Distributions for MicroLED Optical Packaging.** *Appl. Opt.* **2025**, *64*, 7875–7884. https://opg.optica.org/ao/abstract.cfm?uri=ao-64-27-7875.

8. Buhl, M.; et al. **Resonance-Based Directional Light Emission from Organic Light-Emitting Diodes.** *Adv. Photonics Res.* **2023**, *4*, 2200143. https://doi.org/10.1002/adpr.202200143.

9. Abdelkhalik, M. S.; Garcia-Santiago, X.; van Raaij, T.-J.; López, T.; Berghuis, A. M.; de Jong, L. M. A.; Gómez Rivas, J. **Enhanced and Directional Electroluminescence from MicroLEDs Using Metallic or Dielectric Metasurfaces.** *Commun. Eng.* **2025**, *4*, 63. https://doi.org/10.1038/s44172-025-00401-w.

10. Winston, R.; Jiang, L.; Ricketts, M. **Nonimaging Optics: A Tutorial.** *Adv. Opt. Photon.* **2018**, *10*, 484–511.

11. Rau, U. **Reciprocity Relation between Photovoltaic Quantum Efficiency and Electroluminescent Emission of Solar Cells.** *Phys. Rev. B* **2007**, *76*, 085303.

12. Xiang, C.; Koo, W.; So, F.; Sasabe, H.; Kido, J. **A Systematic Study on Efficiency Enhancements in Phosphorescent Green, Red and Blue Microcavity Organic Light-Emitting Devices.** *Light: Sci. Appl.* **2013**, *2*, e74. https://doi.org/10.1038/lsa.2013.30.

13. Song, J.; et al. **Lensfree OLEDs with over 50% External Quantum Efficiency via External Scattering and Horizontally Oriented Emitters.** *Nat. Commun.* **2018**, *9*, 3207. https://doi.org/10.1038/s41467-018-05671-x.

14. Liao, P.-H.; Lee, W.-K.; Lee, C.-C.; Huang, C.-W.; Wen, S.-W.; Chen, Y.-T.; Chen, C.-C.; Lin, W.-Y.; Kwak, B. L.; Visser, R. J.; Wu, C.-C. **Using Angle-Selective Optical Film to Enhance the Light Extraction of a Thin-Film Encapsulated 3D Reflective Pixel for OLED Displays.** *Opt. Express* **2022**, *30*, 46435–46449. https://doi.org/10.1364/OE.477797.

15. Kim, H.-J.; et al. **High Efficient OLED Displays Prepared with the Air-Gapped Bridges on Quantum Dot Patterns for Optical Recycling.** *Sci. Rep.* **2017**, *7*, 43063. https://doi.org/10.1038/srep43063.

16. Kim, J.; Kim, E.; Park, J.; Song, J.; Kim, S.; Moon, H.; Yoo, S. **Toward Near-Foldable Surface Light Sources with Ultimate Efficiency: Ultrathin Substrates Embedded with Micron-Scale Inverted Lens Arrays.** *ACS Photonics* **2023**, *10*, 1775–1782. https://doi.org/10.1021/acsphotonics.3c00017.

17. Kim, M.; Kim, J.; Yoo, S. **Near-Planar Light Outcoupling Structures with Finite Lateral Dimensions for Ultra-Efficient and Optical Crosstalk-Free OLED Displays.** *Nat. Commun.* **2025**, *16*, 11606. https://doi.org/10.1038/s41467-025-66538-6.

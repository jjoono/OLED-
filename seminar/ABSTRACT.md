# 강연 초록 / Seminar Abstract

---

## 국문

### 제목

**늘어나도 왜곡되지 않는 디스플레이 — 하이브리드 오그제틱 메타물질 플랫폼과 등방 팽창 계수**

부제: *Flat에서 Stretchable까지, 그리고 화면의 크기를 조작하는 시대의 새로운 잣대*

### 초록

평평한 화면은 TV와 스마트폰에는 충분했다. 그러나 인체 부착형·인터랙티브 응용이 요구하는 폼팩터는 커브드, 플렉시블, 폴더블, 롤러블을 거쳐 최근 CES에서 시연된 스트레처블로 이어지고 있다. 여기에는 지금까지 없었던 단절이 하나 있다. 커브드부터 롤러블까지는 화면의 **형태**만 바뀌었을 뿐 물리적 면적과 픽셀 간격은 불변이었던 반면, 스트레처블은 폼팩터 역사상 처음으로 **화면 면적 자체가 커진다**. 즉 스트레처블은 기계 구조물이 아니라 여전히 '디스플레이'이며, 화면이 커지는 순간 이미지 왜곡이 부차적 문제가 아닌 1차 문제가 된다. 대부분의 엘라스토머는 양의 푸아송비를 가지므로 가로로 15% 늘리면 세로는 7.5% 수축해 종횡비가 1.24배 일그러진다. 두 손가락 확대(**u** = *k***r**)와 달리, 그냥 늘리는 것은 확대가 아니다.

음의 푸아송비(*ν* = −1)를 갖는 오그제틱 메타물질이 자연스러운 해답으로 제시되어 왔고, 2022년 한국기계연구원이 발표한 첫 오그제틱 메타 디스플레이가 그 출발점이었다. 그러나 키리가미 기반 구조는 구멍이 뚫린 시트가 곧 소자 기판이 되어 공정 호환성과 해상도에 제약이 있었다. 더 근본적으로는, **전체 푸아송비가 −1이라는 사실이 내부의 모든 점이 등방으로 커진다는 것을 전혀 보장하지 않는다**는 문제가 남아 있었다. 회전 기반 단위셀은 국소적으로 등방 스케일링이 아니며, 엘라스토머에 전면 접합하면 회전 궤적이 그대로 전사되고 반대로 엘라스토머는 회전을 구속한다. 게다가 실제로 잡고 당기는 그립 근처에서는 회전장이 지수적으로 감쇠해 회전축이 중앙만이 아니게 된다.

본 세미나에서는 2026년 *Nature Communications*에 발표된 하이브리드 오그제틱 메타물질 플랫폼(Kim, Kim *et al.*, **17**, 7389)을 중심으로 이 문제가 어떻게 풀리는지를 다룬다. 핵심은 두 단계다. 그립의 존재를 목적함수에 포함시켜 프레임 기하를 먼저 최적화하고(*S*(*W*,*L*,*D*) 최소화, *W*/*D* ≈ 35), 그 위에서 **u** = *k***r**을 만족하는 좌표에만 **선택적으로 접합**한다. 프레임은 소자를 얹는 기판이 아니라 등방 팽창의 경계조건을 제공하는 외골격이 되고, 연속 엘라스토머 기판이 그 사이를 매끄럽게 보간한다. 또한 바깥쪽 푸아송비로는 판별할 수 없는 국부 왜곡을 정량화하기 위해 통계의 상관계수로부터 **등방 팽창 계수 *ρ***를 정의하고 *ρ* ≥ 0.99를 지각 임계로 제시한다. FEA와 DIC, 스프레이 프린팅 패턴, 그리고 20% 인장 하의 패시브 매트릭스 LED 디스플레이로 이를 검증하며, 단순 축소가 아닌 strain shielding을 통해 1000 PPI 이상으로 확장하는 경로를 논의한다.

끝으로 이 접근을 비판적으로 검토한다. *ρ*가 전역 상관계수이기 때문에 고주파 국부 왜곡에 둔감하다는 점, 샘플링 영역의 종횡비에 따라 값이 달라진다는 점, 지지 이론이 관측된 최적점을 예측하지 못한다는 점, 그리고 오그제틱 프레임을 제거한 대조군이 빠져 있다는 점을 함께 논의한다. 스트레처블에서 익스팬더블로 넘어가며 화면의 크기를 실제로 조작하게 된 지금, 우리에게 필요한 왜곡의 잣대는 무엇이어야 하는가를 화두로 남긴다.

**키워드:** 스트레처블 디스플레이, 오그제틱 메타물질, 음의 푸아송비, 등방 팽창, 이미지 왜곡, 선택적 접합, 키리가미, 디스플레이 폼팩터

### 공지용 단축본 (약 200자)

커브드부터 롤러블까지 디스플레이 폼팩터는 계속 바뀌어 왔지만, 화면의 물리적 면적이 실제로 커지는 것은 스트레처블이 처음이다. 그 순간 이미지 왜곡이 새로운 1차 문제가 된다. 본 세미나는 음의 푸아송비 오그제틱 구조가 왜 필요조건일 뿐인지—전체가 등방 팽창해도 내부의 모든 점이 그런 것은 아니다—를 짚고, 2026년 *Nature Communications*에 보고된 선택적 접합 기반 하이브리드 오그제틱 플랫폼과 새 평가 지표인 등방 팽창 계수 *ρ*를 다룬다. 끝으로 이 지표와 설계 방법론의 한계를 비판적으로 검토한다.

---

## English

### Title

**Displays That Grow Without Distorting — Hybrid Auxetic Metamaterial Platforms and the Isotropic Expansion Factor**

Subtitle: *From flat to stretchable, and what it takes to measure distortion once screen area itself becomes a variable*

### Abstract

Flat panels were good enough for televisions and phones. But the form factors demanded by skin-attachable and interactive applications have marched from curved to flexible, foldable, and rollable, arriving at the stretchable prototypes recently demonstrated at CES. There is a discontinuity hidden in that sequence. From curved through rollable, only the **shape** of the screen changed; the physical area and the pixel pitch never did. Stretchable displays are the first form factor in which **the screen area itself grows**. A stretchable display is therefore not a mechanical structure but still a *display*, and the moment the screen grows, image distortion stops being a secondary concern and becomes the primary one. Most elastomers have a positive Poisson's ratio, so a 15% stretch along one axis contracts the other by 7.5% and skews the aspect ratio by a factor of 1.24. Unlike a two-finger pinch-zoom, which enforces **u** = *k***r** at every point, simply pulling on a substrate is not magnification.

Auxetic metamaterials with a negative Poisson's ratio (*ν* = −1) have been the natural answer, beginning with the first auxetic meta-display reported by KIMM in 2022. Yet kirigami-based architectures make the void-rich sheet itself the device substrate, which constrains process compatibility and resolution. More fundamentally, **a global Poisson's ratio of −1 guarantees nothing about whether every interior point expands isotropically**. Rotation-based unit cells are not locally isotropic; bonding an elastomer across the entire frame copies the rotational trajectory onto the device surface, while the elastomer in turn suppresses the rotation. Worse, near the grips that must be used to actually pull the device, the rotation field decays exponentially into the array, so the centre of rotation is not where the design assumed it to be.

This seminar follows how these problems are resolved in a hybrid auxetic metamaterial platform reported in *Nature Communications* in 2026 (Kim, Kim *et al.*, **17**, 7389). The methodology has two steps. First, the frame geometry is optimised with the grips written into the objective function — minimising *S*(*W*,*L*,*D*), which yields a distinct optimum near *W*/*D* ≈ 35. Second, the frame is bonded to a continuous elastomeric substrate **only at coordinates that satisfy u = *k*r**. The frame thus becomes an exoskeleton that imposes the boundary conditions for isotropic expansion rather than a substrate carrying the devices, and the continuous membrane interpolates smoothly in between. To quantify the local distortion that an outer Poisson's ratio cannot see, the authors borrow the correlation coefficient from statistics to define the **isotropic expansion factor *ρ***, and propose *ρ* ≥ 0.99 as the threshold of perceptual indistinguishability. The claims are validated by FEA, digital image correlation, spray-printed patterns, and a passive-matrix LED display operated under 20% strain, with a route to beyond 1000 PPI argued through strain shielding rather than structural downscaling.

The talk closes with a critical appraisal. Because *ρ* is a global correlation coefficient, it is insensitive to precisely the high-spatial-frequency local distortion it was introduced to detect; its value depends on the aspect ratio of the sampled region; the supporting theory rationalises rather than predicts the observed optimum; and the control experiment without the auxetic frame is absent. As the field moves from stretchable to expandable displays and screen size becomes something we actively manipulate, the question left open is what the right yardstick for distortion should be.

**Keywords:** stretchable display, auxetic metamaterial, negative Poisson's ratio, isotropic expansion, image distortion, selective bonding, kirigami, display form factor

### Short version for announcements (~120 words)

Display form factors have moved from curved to flexible, foldable, and rollable — but stretchable is the first in which the physical area of the screen actually grows. That is the moment image distortion becomes a first-order problem. This seminar explains why auxetic structures with a negative Poisson's ratio are a necessary but not sufficient condition: a global *ν* = −1 does not mean every interior point expands isotropically. It then covers the hybrid auxetic platform reported in *Nature Communications* (2026), which resolves the deadlock by bonding a rigid auxetic frame to a continuous elastomer only at selected anchor points, and introduces the isotropic expansion factor *ρ* as a field-level distortion metric. The talk ends with a critical look at the limits of that metric and of the design methodology.

---

## 참고 / Reference

Su-Bon Kim, Junho Kim, Sejin Kim, Dongho Choi, Chang-Yeon Gu, Taek-Soo Kim, Hanul Moon & Seunghyup Yoo,
*Hybrid auxetic metamaterial platforms enabling multiscale isotropic expansion for distortion-free stretchable displays*,
**Nature Communications 17, 7389 (2026)**. doi:10.1038/s41467-026-74141-6

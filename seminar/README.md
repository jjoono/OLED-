# 스트레처블 디스플레이 세미나 자료

`stretchable_display_seminar.pptx` — 52장, 발표 시간 약 35분 + Q&A. 모든 슬라이드에 발표자 노트 포함.

## 다루는 논문

Su-Bon Kim, Junho Kim, Sejin Kim, Dongho Choi, Chang-Yeon Gu, Taek-Soo Kim, Hanul Moon & Seunghyup Yoo,
*Hybrid auxetic metamaterial platforms enabling multiscale isotropic expansion for distortion-free stretchable displays*,
**Nature Communications 17, 7389 (2026)**. doi:10.1038/s41467-026-74141-6 (KAIST School of EE · Dong-A University)

## 구성

| 부 | 주제 | 슬라이드 |
|---|---|---|
| 01 | 디스플레이 폼팩터는 왜 계속 변하는가 — Flat → Curved → Flexible → Foldable → Rollable → Stretchable | 3–7 |
| 02 | 스트레처블은 결국 '디스플레이'다 — 화면이 커지는 순간 생기는 이미지 왜곡 | 8–13 |
| 03 | Auxetic: 음의 푸아송비라는 답, 그리고 그 계보 (2022 KIMM 최초 auxetic meta-display 포함) | 14–18 |
| 04 | 진짜 문제: Global NPR ≠ Local isotropic expansion — 회전 메커니즘·이중 구속·그립 경계층 | 19–23 |
| 05 | 이 논문의 설계 방법론 — S(W,L,D) 기하 최적화 → 선택적 접합점 설계 | 24–35 |
| 06 | 새 지표: 등방 팽창 계수 ρ — 통계 상관계수에서 빌려온 잣대 | 36–39 |
| 07 | 실험 검증 — 제작·계면 신뢰성·FEA vs DIC·스프레이 패턴·PM LED 디스플레이 | 40–45 |
| 08 | 확장성과 의의 — downsizing의 벽, strain shielding, >1000 PPI, 선행 연구 비교, 토의 | 46–52 |

## 스토리라인 요약

1. 커브드부터 롤러블까지는 **형태**만 바뀌었고 화면 면적과 픽셀 간격은 불변이었다. 스트레처블은 폼팩터 역사상 처음으로 **면적 자체가 커지는** 단계다.
2. 그래서 이미지 왜곡이 새로운 1차 문제가 된다. 엘라스토머는 ν ≈ +0.5 이므로 15% 인장 시 종횡비가 1.24배 일그러진다. 디지털 핑거 줌(u = k·r)과 달리, 그냥 늘리는 것은 확대가 아니다.
3. Auxetic(ν = −1)이 답처럼 보이지만 **필요조건일 뿐이다**. 2022년 KIMM의 첫 auxetic meta-display는 개념적으로 훌륭했으나 void-rich kirigami 시트가 곧 소자 기판이라 공정 호환성·해상도에 제약이 있었다.
4. 더 근본적으로, rotating-square의 메커니즘은 **회전**이라 국소적으로 등방 스케일링이 아니다. 전면 접합하면 회전이 기판에 그대로 복사되고, 반대로 기판은 회전을 구속한다. 게다가 **잡고 당기는 그립** 근처에서 회전장이 지수적으로 감쇠해 회전축이 중앙만이 아니게 된다.
5. 이 논문은 그립까지 목적함수에 넣어 기하를 먼저 최적화하고(W/D ≈ 35), 그 위에서 **u = k·r 을 만족하는 지점에만 선택적으로 접합**한다. 순서를 뒤집으면 유효 연결점이 아예 존재하지 않는다.
6. 그리고 바깥쪽 푸아송비로는 판별할 수 없는 국부 왜곡을 재기 위해 통계의 상관계수에서 **등방 팽창 계수 ρ**를 정의했다. ρ ≥ 0.99가 지각 임계다.

## 자료 출처

논문 Figure는 원저작물(CC BY-NC-ND 4.0)에서 인용했으며, 슬라이드에 '개념도 (본 발표용 작성)'으로 표기된 그림·그래프는 발표용으로 새로 제작한 것이다.

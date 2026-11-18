# 두 원고의 타겟 저널 제안

- 작성일: 2026-11-18
- 근거: 각 저널의 범위와 일반적인 채택 경향을 바탕으로 한 판단이다. 최근 게재 동향이나 편집 방침을 조회해 확인한 것은 아니다.

| 원고 | 도전적 | 현실적 (주력) | 안전망 |
|---|---|---|---|
| **A. Reciprocity와 angular mixing** (`manuscript_draft.md`) | Optica | Optics Express | JOSA B, Optics Letters(축약판) |
| **B. Round-trip loss** (Nature 원고) | Nature, Nature Photonics | Nature Communications, Light: Science & Applications | Science Advances, Advanced Optical Materials, ACS Photonics |

## A. Reciprocity와 angular mixing 원고

**현실적: Optics Express**
- 이미 OE 형식으로 준비되어 있다(`oe_submission/`).
- 이론 유도, 축소 모델, 수치 검증, 소규모 실험으로 이어지는 구성은 OE가 선호하는 "광학 모델 + 엄밀한 검증" 논문이다.
- 검토 보고서의 주요 문제(이상적 투과 혼합의 실현 가능성, 초록·결론의 과장)를 고치면 채택 가능성이 높다.

**도전적: Optica**
- "OLED 광추출 모델"이 아니라 일반 원리로 포장해야 한다. 즉 "상반성이 반사 혼합에 분산만큼의 벌점을 부과하고, statistical ray optics는 투과 혼합에만 맞는다"는 메시지를 LED와 태양전지 light trapping으로 확장한다. $n$=3.5 slab 결과가 그 근거다.
- 약점은 두 가지다. 실험(4.5절)이 핵심 주장을 검증하지 않고, 실제 성능 이득이 "최대 약 2 pp"라는 부정적 결과다. 편집자 단계에서 반려될 위험이 있다.
- **Physical Review Applied**도 대안이다. 이론 비중이 크고 소자 응용이 있는 논문을 잘 받는다.

**전략**
- B 원고가 A의 출발점이다. A는 [14]로 B를 인용하고 η 식을 B에서 가져온다.
- 따라서 B가 먼저 출판되거나 최소한 arXiv에 올라가 있어야 A를 심사받기 수월하다.
- B가 높은 저널에 실리면, A는 OE에서 "후속 정밀 이론"으로 자리 잡기 좋다.

## B. Round-trip loss 원고

**도전적: Nature, 이어서 Nature Photonics**
- 메시지 자체(50–60% 한계 원인 규명, 77% 실증, 90%로 가는 길)는 Nature급이다.
- 다만 현재 형태로는 고전할 가능성이 크다. 활성 면적 대비 반사체 면적이 60배이고, 76.7%가 소자 1개의 결과이며, 새로움의 범위가 [30] Levell, [6]과 겹친다.
- Nature에서 반려되면 Nature Photonics로 이관하는 것이 자연스럽다.

**현실적: Nature Communications 또는 Light: Science & Applications**
- 같은 그룹의 직전 연구가 Nat. Commun.에 실렸다([6] 2018, [23] 2025). 비교 문헌의 수준과 독자층이 맞는다.
- 대면적 실증 없이도 "메커니즘 규명 + 77% + 설계 규칙"이면 충분히 경쟁력이 있다.

## 등급을 올리려면

- **B를 Nature 또는 Nature Photonics에 내려면:** 활성 면적과 반사체 면적이 비슷한 소자에서 70%를 넘기는 결과(예: Ag 음극과 저손실 반사체를 소자 안에 통합)가 사실상 필수다. 조건마다 소자 여러 개의 통계와 적분구 교차검증도 필요하다. 이것이 없으면 Nat. Commun.을 1순위로 두는 편이 시간 면에서 유리하다.
- **A를 Optica에 내려면:** 투과 혼합 기준선에 실제로 접근하는 구조를 하나 시뮬레이션이나 실험으로 보여야 한다("이렇게 하면 이득을 얻는다"는 긍정적 결과). 태양전지 light trapping 쪽 함의도 한 문단 이상 다뤄야 한다.

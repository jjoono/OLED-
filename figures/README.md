# Figures

모두 순수 벡터 SVG이고 SVG 필터를 쓰지 않으므로, Illustrator/Inkscape에서
래스터화 없이 그대로 편집됩니다.

## `tandem_blue_pep_oled.svg`

Tandem blue PEP-OLED 소자 구조 모식도 (아이소메트릭 분해도).
순수 벡터 SVG이므로 해상도 제한이 없고, 논문/발표자료에 그대로 삽입할 수 있습니다.

- 캔버스: 440 x 520 (viewBox 동일) — 어떤 배율로도 깨지지 않습니다.
- 층 구성 (아래 → 위): ITO/Ag/ITO anode · Polaritonic HTL · Blue EML · CGL ·
  Blue EML · Polaritonic ETL · Ag/Al/Liq cathode
- 발광은 두 EML 유닛 바깥쪽 계면(HTL 위, ETL 아래)에서 나오는 것으로 표현했습니다.

### 편집 방법

| 프로그램 | 비고 |
|---|---|
| **Inkscape** (무료, Win/Mac/Linux) | SVG가 네이티브 포맷이라 손실 없이 열리고 저장됩니다. 가장 권장 |
| **Adobe Illustrator** | 열 때 "Convert to CSS / Use system fonts" 대신 **Cascading Style Sheets 유지** 옵션을 쓰면 그라디언트가 보존됩니다 |
| **Affinity Designer / Figma** | 그라디언트·블러 대부분 그대로 들어옵니다 |
| **웹 브라우저** | 보기 전용 (파일을 드래그해서 열기) |
| **텍스트 에디터** | SVG는 XML이라 색상 코드나 라벨을 직접 고칠 수 있습니다 |

SVG 안의 주요 그룹에 id를 달아두었으므로 레이어 패널에서 바로 찾을 수 있습니다:
`bottom-unit`, `glow-lower`, `middle-unit`, `glow-upper`, `top-unit`, `labels`,
`ground-shadow`.

### 스크립트로 다시 만들기

`make_tandem_figure.py` 가 이 SVG를 생성합니다. 층 두께·색·간격을 바꾸려면
스크립트 상단의 투영 벡터(`EX/EY`, `DX/DY`)와 `layers` 리스트, `<defs>`의
그라디언트만 고치면 됩니다.

```bash
python3 figures/make_tandem_figure.py
```

---

## `trans_scale_simulation.svg`

Trans-scale 광학 시뮬레이션 개념도 — 컷어웨이(cut-away) 아이소메트릭.
소자 앞모서리를 사각으로 잘라내어 단면을 드러내고, 그 단면 위에 물리를 그렸습니다.

- 캔버스: 1020 x 565
- 층 구성 (위 → 아래): Metal cathode · OLED · (ITO anode) · Glass ·
  Outcoupling structure (아랫면에 마이크로렌즈 어레이)
- 단면 위 물리 표현
  - **파동광학 영역**: 쌍극자 방사 로브 2개 + 쌍극자 모멘트 화살표 + 퍼져나가는 파면
  - **기하광학 영역**: 유리를 가로지르는 광선, 원거리장 세기 분포(반구 극좌표 + 적색 로브),
    산란 입자와 산란 경로, 마이크로렌즈를 통한 적색 추출광
- 오른쪽: 두 영역을 묶는 브레이스와 노란 점선 패널

레이어 패널에서 찾을 수 있는 id: `layer-Cath`, `layer-Oled`, `layer-Ito`,
`layer-Glass`, `layer-Out`.

### 구조 바꾸기

`make_transscale_figure.py` 상단에서 조절합니다.

| 변수 | 의미 |
|---|---|
| `EX, EY` / `DX, DY` | 아이소메트릭 투영 벡터 (보는 각도) |
| `U0, V0` | 잘라낸 모서리의 위치 (0~1). 값을 키우면 단면이 좁아집니다 |
| `LAYERS` | `(이름, 윗면 높이, 두께, 평면 크기 배율)` — 층 추가·삭제·두께 조절 |

```bash
python3 figures/make_transscale_figure.py
```

---

## `outcoupling_roundtrip.svg`

그림 1(a) — 외부 광추출 구조가 부착된 OLED에서 빛이 소자와 공기 사이를
여러 번 왕복하는 과정. 2차원 단면이며, 왼쪽/오른쪽 두 패널을 비교합니다.

- 캔버스: 2120 x 442
- 원본과 달리 `⋮` 로 끊지 않고, **금속 음극부터 마이크로렌즈까지 한 장에** 그렸습니다.
- 왼쪽: 일반적인 소자. 금속 ohmic 흡수·TCO 흡수가 커서 몇 번의 왕복 만에 빛이 소멸
- 오른쪽: 제안하는 디자인 룰로 설계한 소자. 흡수가 억제되어 세기가 유지됨

### 광선 기하 (스크립트에서 검증됨)

| 항목 | 값 |
|---|---|
| 유리 내부 전파각 | 45° (유리 n=1.5 의 임계각 41.8° 초과 → 평탄면이면 전반사로 갇힘) |
| 반 왕복당 수평 이동 | 120 = 렌즈 피치의 정확히 2배 |
| 결과 | 네 번의 출광 시도가 모두 렌즈 **정중앙**에 입사 |

화살표의 선 굵기는 광출력의 0.75 제곱에 비례합니다. 한 패널 안에서 굵기가
줄어드는 정도가 곧 왕복당 손실입니다.

| | 금속 반사율 | TCO 투과율(1회) | 누적 추출광 |
|---|---|---|---|
| 일반 소자 | 0.72 | 0.90 | 49 % |
| 디자인 룰 | 0.98 | 0.995 | 74 % |

`make_roundtrip_figure.py` 상단의 `panel(...)` 호출에서 이 세 값
(`eta`, `r_met`, `t_tco`)만 바꾸면 화살표 굵기와 손실 표시 크기가 전부
자동으로 다시 계산됩니다.

```bash
python3 figures/make_roundtrip_figure.py
```

---

세 스크립트 모두 의존성 없이 파이썬 표준 라이브러리만 사용합니다.

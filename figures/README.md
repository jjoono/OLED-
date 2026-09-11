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

두 스크립트 모두 의존성 없이 파이썬 표준 라이브러리만 사용합니다.

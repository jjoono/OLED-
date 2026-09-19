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

### 색상 — Nature Publishing Group 계열

`make_roundtrip_figure.py` 상단 팔레트 블록 한 곳에서 관리하며, PPTX 빌더가
이 상수들을 그대로 import 하므로 두 포맷의 색이 어긋날 수 없습니다.

| 역할 | 색 | 비고 |
|---|---|---|
| 광선 | `#00a087` | NPG teal-green |
| 흡수 손실 | `#e64b35` | NPG coral red |
| 발광점 | `#e8901f` | gold |
| 유리/렌즈 | `#eaf1f6` / `#6e93ae` | 구조는 전부 채도를 낮춤 |
| TCO | `#cbe1ee` / `#6e93ae` | |
| 유기층 | `#fdf2e2` / `#c99a55` | |
| 금속 음극 | `#c3c7cb` / `#868c92` | |

채도가 높은 색은 **빛(teal)과 손실(coral) 둘뿐**이고 나머지 구조는 모두
저채도 틴트입니다. 선 두께도 Nature 인쇄 기준(0.4–0.6 pt)에 맞췄습니다.

### PowerPoint 버전 — `outcoupling_roundtrip.pptx`

**그룹이 하나도 없고 이미지도 하나도 없는** 132개의 네이티브 PPT 도형입니다.
열자마자 개별 클릭·편집이 됩니다.

| 그림 요소 | PPT 도형 | 편집 방법 |
|---|---|---|
| 마이크로렌즈 15개 | 타원(Oval) | 개별 선택해 크기·개수 변경 |
| 유리/TCO/유기층/금속 | 직사각형 | 높이를 끌어 두께 조절 |
| 광선·출광 화살표 | 직선 화살표 커넥터 | 끝점을 끌어 각도 변경, 선 두께로 세기 표현 |
| 발광점·흡수 버스트 | 별 도형 | 꼭짓점 수·크기 조절 |
| 흡수 물결 | 자유형 | 점 편집 가능 |
| 모든 글자 | 텍스트 상자 | 바로 타이핑 |

- 슬라이드 크기는 그림 비율(13.33 × 2.90 in)에 맞춰 두었습니다.
  `디자인 > 슬라이드 크기`에서 바꾸거나, 전체 선택 후 다른 덱에 붙여넣으면 됩니다.
- 도형마다 이름이 있어 `홈 > 선택 > 선택 창`에서 바로 찾을 수 있습니다
  (`L lens 03`, `R ray up 2`, `L ohmic loss 1L` …).

```bash
python3 figures/make_roundtrip_pptx.py     # needs python-pptx
```

> 나머지 두 3D 그림은 그라디언트를 많이 쓰므로, PowerPoint에서
> `삽입 > 그림`으로 SVG를 넣은 뒤 `그래픽 형식 > 도형으로 변환` → `그룹 해제`
> 를 하면 그라디언트를 유지한 채 네이티브 도형으로 바뀝니다.

---

## `emitting_area.svg`

저손실 설계 소자와 일반 소자의 **발광 면적 · 밝기** 비교. 아이소메트릭 스택 2개.

- 캔버스: 1200 x 492
- 왼쪽 = 일반 소자, 오른쪽 = 디자인 룰 적용 (Fig.1(a)와 좌우 순서를 맞춤)

### 발광 영역 크기를 눈대중으로 정하지 않았습니다

Fig.1(a)와 **같은 왕복 손실 모델**에서 두 수치를 유도해 그대로 씁니다.

```
rho = (1-eta) * R_metal * T_TCO^2     한 왕복당 살아남는 비율
Lam = 1 / ln(1/rho)                   측방 확산 길이
tot = eta / (1-rho)                   최종적으로 빠져나오는 총 광량
```

| | rho | 확산 길이 | 총 추출광 |
|---|---|---|---|
| 일반 소자 | 0.408 | 1.12 step | 51 % |
| 디자인 룰 | 0.679 | 2.58 step | 94 % |

→ **발광 영역 반지름 = 확산 길이 비(2.32×)**,
**발광 밝기 = 총 추출광 비(1.84×)** 로 대응시켰습니다.
(Fig.1(a)의 49 % / 74 % 는 왕복 4회까지만 센 값, 위 51 % / 94 % 는 무한 합계입니다.)

발광은 표면 위의 **원**이므로, 등방성 확산이 아이소메트릭 평면에서 타원으로
보이도록 radial gradient에 `gradientTransform` 으로 같은 투영 기저를 실었습니다.
손으로 타원을 그린 것이 아닙니다.

```bash
python3 figures/make_emitting_area_figure.py
```

팔레트는 `make_roundtrip_figure.py` 에서 import 하므로 Fig.1(a)와 색이 어긋나지
않습니다.

### Blender 버전 — `emitting_area_render.blend`

같은 비교를 실제 3D 렌더로 만든 Blender 4.0 파일입니다. 회색 스튜디오 바닥 위에
소자 스택 두 개, 각 상면에 발광 영역과 그 위로 퍼지는 빛 안개(볼륨) 콘.

- **발광 영역의 반지름 = 확산 길이 비(2.32×)**, **발광·안개 세기 = 총 추출광 비(1.84×)**
  — 2D 그림과 같은 손실 모델에서 나온 값을 그대로 씁니다.
- 발광은 `Emitting area (...)` 평면의 재질 노드에서 `exp(-r/Λ)` 로 감쇠합니다.
  Λ 값(`Divide` 노드)과 세기(`Multiply` 노드)만 고치면 됩니다.
- 빛 안개는 `Light haze (...)` 콘의 Principled Volume. 높이에 따라 옅어지고
  콘 벽 쪽으로 부드럽게 사라집니다. 밀도·발광은 재질의 마지막 두 `Multiply` 노드.
- **발광 색** — 스크립트 상단 `GLOW` (선형 RGB, 기본 `(0.06, 1.0, 0.30)` = 쨍한 초록).
  안개가 산란시키는 색은 `HAZE_SCATTER` 로 따로 두어, 흰 조명이 안개에 반사될 때
  너무 진해지지 않게 했습니다. 색을 바꾸면 두 상수만 고치면 됩니다.
- **마이크로렌즈 어레이 (실제 지오메트리)** — 광추출 필름 위에 구면 캡을
  **육각(삼각) 격자**로 깔았습니다. 행 간격 `pitch*sqrt(3)/2`, 행마다 반 피치씩
  어긋나는 실제 MLA 배열입니다. 소자당 201개, 하나의 메시로 합쳐져 있습니다
  (`Micro-lens array (...)`).

  | 상수 | 뜻 |
  |---|---|
  | `LENS_PITCH` | 렌즈 중심 간격 (0.125 → 가로 16개) |
  | `LENS_FILL` | 렌즈 반지름 / 반 피치. 1.0 이면 서로 맞닿음 |
  | `LENS_SINK` | 캡이 필름에 잠기는 깊이 (반지름 대비). 반구가 아니라 캡으로 보이게 함 |
  | `LENS_SEGS`, `LENS_RINGS` | 돔 하나의 분할 수 |

  렌즈 재질만 `shadow_through=False` 입니다. 다른 투과 층들은 그림자 광선을
  통과시키지만, 렌즈는 서로 그림자를 드리워야 형상이 읽히기 때문입니다.
  투과율이 높고 그림자도 없으면 흰 필름 위의 흰 돔이라 **아예 보이지 않습니다**.
- **발광면 위치** — 발광 평면은 렌즈 어레이 **아래**(필름 상면)에 있습니다.
  빛이 렌즈를 통과해 나오므로, 렌즈 격자가 초록빛 위로 비쳐 보입니다.
- 유리 층은 그림자 광선을 통과시키는 Light Path 트릭이 걸려 있어, Cycles에서
  굴절 유리 아래가 검게 나오는 문제가 없습니다.
- 컬렉션: `Device conventional`, `Device design rule`, `Studio`, `Cameras`
  (카메라 3개: `Cam both`, `Cam conventional`, `Cam design rule`).

### 투명 배경 (기본값)

`--bg transparent` 가 기본이라 렌더는 **알파 채널이 있는 PNG**로 나옵니다.
바닥은 지워지지 않고 shadow catcher 로 바뀌므로, 배경은 비면서 접지 그림자는
남습니다. 회색 스튜디오 배경이 필요하면 `--bg studio`.

한 가지 주의할 점이 있어 컴포지터가 자동으로 붙습니다. Cycles 에서 옅은 발광
볼륨은 알파가 거의 0 이라, 투명 배경으로 렌더하면 **빛 안개가 RGB 에는 있는데
합성하면 사라집니다.** 그래서 컴포지터에서

- 알파 = `max(원래 알파, 채널 최댓값)` 로 다시 만들고,
- 원래 알파가 0 이던 곳(= 안개)만 색을 un-premultiply 합니다.

휘도가 아니라 **채널 최댓값**을 쓰는 이유가 있습니다. 채도 높은 색을 휘도로
나누면 강한 채널이 1 을 넘겨 클리핑되면서 색이 흰색으로 돌아갑니다. 채널
최댓값으로 나누면 색조가 그대로 보존됩니다.

덕분에 밝은 배경에도 어두운 배경에도 초록빛이 그대로 얹힙니다.
이 처리가 싫으면 `Compositing` 노드 트리를 끄면 됩니다.

같은 노드 트리 맨 앞의 `Hue Saturation Value` 가 채도를 1.30 배 올립니다.
이 Blender 빌드는 OCIO 룩(`AgX - Punchy` 등)을 노출하지 않아서 직접 넣었습니다.
AgX 는 밝은 부분의 채도를 깎기 때문에, 발광 세기(`EMIT_STRENGTH`)를 올려
코어를 날리면 오히려 색이 흰색으로 빠집니다 — 색을 진하게 하고 싶으면
세기가 아니라 `GLOW` 와 이 채도 값을 조절하세요.

### 렌더 품질

| 항목 | 값 | 이유 |
|---|---|---|
| 해상도 | 2400 x 1200 (클로즈업 1200 x 1200) | Nature 2단 폭 180 mm @ 300 dpi |
| 샘플 | 320, adaptive threshold 0.006 | OIDN 이 뒤를 받치므로 과하게 올릴 필요 없음 |
| 바운스 | total 16 / transmission 16 / glossy 8 | 렌즈-유리-유리로 투과가 겹침 |
| 코스틱 | **끔** | 투과 캡 400개에서 렌더 시간이 3배가 되는데 눈에 띄는 차이가 없었음 |
| 픽셀 필터 | 1.20 (기본 1.5) | 조금 더 또렷하게 |
| 뷰 트랜스폼 | AgX + 컴포지터 채도 1.30 | |

`--samples N`, `--res PCT` 로 임시 오버라이드할 수 있어 빠른 확인에 씁니다.

### 미감 관련 설정

- **모서리 베벨** — 모든 층 박스에 폭 0.0035 / 4 세그먼트 베벨이 걸려 있습니다.
  완벽하게 날카로운 모서리는 CG 티가 가장 많이 나는 요소입니다.
- **환경광은 옆으로 기운 그라디언트**입니다. 천정이 밝으면 위를 향한 렌즈 캡들이
  전부 같은 양의 빛을 받아 **배열의 음영이 통째로 사라집니다.** 실제로 이 문제로
  렌즈가 안 보였고, 환경광을 어둡게(0.11) 하고 좌우 방향으로 바꿔 해결했습니다.
- **렌즈 캡이 직접 발광합니다** (`lens_material`). 빛은 렌즈를 통과해 나오므로
  물리적으로도 맞고, 덕분에 캡을 거의 불투명하게(투과 0.08) 두어 음영을 살리면서도
  초록이 나옵니다. 발광 평면과 같은 `radial_falloff` 노드를 공유하므로 어긋나지
  않습니다. 비율은 `LENS_EMIT_FRAC`.
- **키 라이트는 크고 부드럽게**(size 7.5). 작고 단단한 광원은 렌즈 어레이의 그림자를
  바닥에 빗살처럼 찍어냅니다.
- 카메라는 85 mm 로 멀리서 — 거의 정사영에 가까워 기술 도면 같은 깔끔함이 납니다.

렌더는 Cycles(CPU). 빌드에 OpenImageDenoise가 있으면 켜지고, 없으면
`denoise_renders.py` 로 외부 디노이즈합니다(PyPI `oidn`, `OpenEXR` 패키지;
빌더가 알베도·노멀 패스를 EXR로 같이 저장해 둡니다).

```bash
blender -b -P figures/build_emitting_area_blend.py -- --out figures                  # 씬 + 렌더 3장
blender -b -P figures/build_emitting_area_blend.py -- --out figures --no-render      # 씬만
blender -b -P figures/build_emitting_area_blend.py -- --out figures --bg studio      # 회색 배경
blender -b -P figures/build_emitting_area_blend.py -- --out . --quality test --cam conventional
python3 figures/denoise_renders.py figures/emitting_area_render*.png
```

`--cam` 은 이름 일부만 줘도 됩니다(`both`, `conventional`, `design`).
디노이저는 알파 채널도 같이 정리합니다 — 컴포지터가 글로우의 노이즈를
알파 쪽으로 옮겨 놓기 때문입니다.

렌더 결과: `emitting_area_render.png` (두 소자, 2400 x 1200), `..._conventional.png`,
`..._designrule.png`. 라벨·인셋 도식은 렌더에 넣지 않았으니 Illustrator/PowerPoint에서
얹으시면 됩니다.

---

2D 스크립트들은 파이썬 표준 라이브러리만 사용합니다
(PPTX 빌더는 `python-pptx`, Blender 빌더는 Blender 4.0 내장 파이썬, 디노이저는 `oidn`+`OpenEXR`).

# 비대칭 3D Freeform MLA 최적화 (BO 기반 trans-scale)

hemisphere/회전대칭 렌즈로는 **원리적으로 불가능한 방향성 비대칭 발광**을, 모든 렌즈가
동일한 비대칭 freeform 형상을 갖고 같은 방향으로 정렬된 MLA로 실현하기 위한 설계 코드.

## 지오메트리 방식 — LightTools 네이티브 Freeform (.ent 직접 생성)

LightTools **FreeformEntity(.ent)** 라이브러리 파일은 표면 격자점을 X,Y,Z 텍스트로
담는다(`restorePoints` 블록). 따라서:

- MATLAB이 검증된 템플릿 `.ent`에서 **FrontSurface 격자의 Z만 바꿔 `.ent`를 직접 작성**
- 배열 모델 텍스처의 `LibraryElement` `Filename`을 그 `.ent`로 지정 (기존 코드가
  `swept_XXX.1.ent`를 물리던 자리와 동일)
- **STL/mesh/import/SaveLibrary/SweptEntity 인스턴스 전부 불필요**

> 이 방식으로 온 이유: 이 LightTools 설치본은 STL 임포트가 없고 CAD 임포트(SAT/STEP/
> IGES)는 비활성이라, 외부 지오메트리를 넣을 길이 `.ent` 직접 작성뿐이었다. 다행히
> FreeformEntity가 격자 sag를 네이티브로 지원해 정확히 들어맞았다.

## 파일 구성

| 파일 | 역할 |
|------|------|
| `BO_Freeform3D_asym.m` | 메인 드라이버. **이 파일을 Run.** BO로 freeform 격자 Z + 마이크로캐비티 최적화. 목적함수 = 목표 방향 조향 EQE(`EQE_cone`). |
| `generate_freeform_ent.m` | 템플릿 `.ent`에서 FrontSurface 격자 Z만 치환해 새 `.ent` 작성 (X,Y·구조·RearSurface 평면 그대로). |
| `freeform_grid_info.m` | 템플릿 격자의 (X,Y)와 내부점 인덱스 파싱. |
| `freeform_template.ent` | **사용자가 GUI에서 만든 "유효 solid" freeform 렌즈.** Front/Rear 규칙 5×5, 두께 1mm, 경계 Z=0, NURBS Off. |
| `RenewLightTools.m` | 배열 모델 LightTools 1개 재시작/연결 (경로·사용자명 여기서 수정). |

## 파라미터화 & 형상

```
z1..z9      : FrontSurface 내부 3x3 격자 높이 (경계 16점은 Z=0 고정)   ← 형상 DOF
dETL, dHTL  : OLED 마이크로캐비티 두께 (나노 스케일)
stretchZ    : 텍스처 z 스트레치
```
= 9 형상 DOF + 3 = **12 DOF** (GP-bayesopt에 적합).

- 렌즈 = plano-freeform: **FrontSurface = 자유형 곡면, RearSurface = 평면**(기판에 파묻힌 바닥).
- 경계 Z=0 → bump가 기판에 매끄럽게 맞물림. 내부 3×3만 자유 → 비대칭/조향은 내부
  높이 분포로 만든다(예: +x쪽을 높이면 발광이 +x로 조향).
- **중심 두께 1mm > 내부 Z 상한 0.8** → 두 표면이 안 겹쳐 항상 유효 solid.

### [중요] NURBS는 Off로 유지
템플릿은 `restoreNURBS: "No"`. **NURBS를 켜면 LightTools가 U,V를 바꿔 면을
재구성**해 격자가 흔들린다. 생성기는 템플릿의 NURBS Off를 그대로 두므로 격자가 안정적
5×5로 고정된다(`restoreSmoothResample: "Yes"`가 매끄러움은 유지). 격자를 키우려면
템플릿 자체를 GUI에서 원하는 U×V로 다시 만들어 저장하면, 코드가 자동으로 그 크기에
맞춘다(`freeform_grid_info`가 N·내부점을 파싱).

## 목적함수: 특정 (θ밴드 × φ구간) EQE 최대화

`J = EQE_region` = 목표 **극각 밴드 `θ∈[TGT_TH_LO,TGT_TH_HI]` × 방위각 구간
`φ∈[TGT_PHI_C ± TGT_PHI_HALF]`** 로 방출된 EQE. (대비 페널티 없음 — hemisphere는
별도 대조군으로 비교.) 영역이 φ로 국한돼 있어 EQE_region 최대화 자체가 방위각 집중을
보상하고, hemisphere는 φ 전체로 퍼져 같은 영역에 조금만 넣으므로 비교에서 낮게 나온다.

- 리포트: `EQE_total`, `EQE_region`, `directionality = EQE_region/EQE_total`.
- 검증(합성 원거리장, θ∈[20,40] φ∈[0±45]): 비대칭 tilt=0.35 vs hemisphere=0.07 (~5×).
- **φ_c = 렌즈 정렬 방향**. 배열을 Δφ 회전 → 발광도 Δφ 회전(방향 제어). θ밴드 중심을
  바꿔 달성 가능한 최대 tilt를 특성화. hemisphere를 같은 조건에서 돌려 대조군으로.

## 실행 전 확인 (@@VERIFY)

Windows + LightTools 2023.03 + MATLAB(COM), Global/Statistics Toolbox 필요.

1. **템플릿 배치**: `freeform_template.ent`(유효 solid freeform 렌즈)를
   `RenewLightTools.m`/`BO_Freeform3D_asym.m`의 `FF_BASE` 경로에 둘 것.
2. **텍스처 LibraryElement**: 배열 모델(Lens_size_effect...lts)의 텍스처 unit-cell이
   이 freeform `.ent`를 받도록 설정돼 있어야 함(기존 swept `.ent` 자리 대체).
3. **far-field mesh 차원** `MESH_NLONG/MESH_NLAT`을 모델 실제 셀 수에 맞추고,
   **Symmetry = No Symmetry** 로 둘 것(안 그러면 방위각 비대칭이 평균화되어 안 보임).
4. 경로/파일명/사용자명(`jhkim`), `.mat` 데이터, CPS 관련 함수가 경로에 있을 것.

## 방향 제어성 증명 그림 (권장)
`target_phi`를 0:45:315로 스윕하며 각각 독립 최적화 → 발광 피크가 그 방위각을 따라
도는지 보이면 "동일 렌즈의 회전정렬만으로 발광 방향 제어" 주장을 직접 입증.

## 계승 (v4)
탐색/검증 ray 분리, 크래시 시 NaN(GP 오염 방지), warm-start, 수렴 적응예산,
patternsearch 정련, 고정밀 반복검증 mean±std.

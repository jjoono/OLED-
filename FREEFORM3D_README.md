# 비대칭 3D Freeform MLA 최적화 (BO 기반 trans-scale)

hemisphere/회전대칭 렌즈로는 **원리적으로 불가능한 방향성 비대칭 발광**을, 모든 렌즈가
동일한 비대칭 freeform 형상을 갖고 같은 방향으로 정렬된 MLA로 실현하기 위한 설계 코드.

## 파일 구성

| 파일 | 역할 |
|------|------|
| `freeform_height.m` | 비대칭 렌즈 높이장 `z(ρ,φ)=H·P(ρ)·[1+Σ ρ^m(c_m cosmφ+s_m sinmφ)]`. **c1,s1(1차 harmonic)이 정점 tilt = 방향성 핵심 DOF**. 회전대칭이면 harmonic 항이 0(= hemisphere가 갇히는 상태). |
| `generate_freeform_mesh.m` | 높이장 → **watertight 삼각 메쉬** → ASCII STL. plano-convex(평평한 바닥) 닫힌 solid. LightTools import용. |
| `BO_Freeform3D_asym.m` | 메인 드라이버. Bayesian Optimization으로 freeform + 마이크로캐비티 동시 최적화. 목적함수 = **목표 방향으로 조향 추출된 EQE (`EQE_cone`)**. |

## 파라미터화 (13-dim)

```
p1..p5      : 반경 base 프로파일 제어높이 (정점=1, rim=0 고정)
H           : 렌즈 높이/aspect
c1, s1      : 1차 harmonic  (정점 tilt = 방위각 비대칭·방향 제어)   ← 논문 core
c2, s2      : 2차 harmonic  (rim shaping)
dETL, dHTL  : OLED 마이크로캐비티 두께 (나노 스케일, EQE에 직접)
stretchZ    : 텍스처 z 스트레치
```

## 목적함수: 왜 hemisphere를 배제하는가

`EQE_cone` = 목표 방향 `(θ_t, φ_t)` 중심 반각 `cone_half` 원뿔로 방출된 EQE.
검증 결과(합성 원거리장):

- (θ=30°,φ=0°)로 tilt된 빔 → 목표 (30,0) 분율 **0.478** vs 반대편 (30,180) **0.006**
- 대칭 빔(θ=0°) → (30,0)=(30,180)=**0.170 동일**

즉 방위각 대칭 발광(hemisphere가 만드는 것)은 어떤 방향을 목표로 해도 대칭이라
`EQE_cone`을 올릴 수 없다. 오직 비대칭 freeform만 특정 방위각으로 조향 가능.

## Trans-scale 흐름 (1회 평가)

```
BO 점 x → freeform_height/generate_freeform_mesh → STL
      → [LightTools swept 인스턴스] import + RepairEntities + SaveLibrary → .ent(ACIS)
      → [LightTools 배열 인스턴스] 텍스처 unit-cell 파일 교체 + StretchZ
   나노: CPS_for_Isub(OLED 스택) → 각도별 I_sub, EQE_sub, 하단 반사율 → .coa 코팅
      → 파장별 소스 apodizer 주입 → BeginAllSimulations
      → far-field INTENSITY_MESH를 (θ×φ) 2D로 읽어 cone_power_fraction → EQE_cone
```

## 실행 전 확인 필요 (`@@VERIFY` 표시)

이 코드는 **Windows + LightTools 2023.03 + MATLAB(COM)** 환경에서 실행된다. LightTools
DB access 이름은 모델/버전에 따라 다를 수 있으므로 아래를 GUI/매크로 기록으로 확인:

1. **STL import 명령** (`updateFreeformGeometry`의 `ImportCADFile`): File > Import를
   매크로 기록해 정확한 LTCmd 문자열/임포트된 엔티티 이름 확보.
2. **far-field mesh 차원** (`read_intensity_grid`의 `nLat=90, nLong=72`): 모델의
   far-field `INTENSITY_MESH` 실제 셀 수로 맞출 것. **longitude 셀 수 > 1 이고
   Symmetry = No Symmetry** 여야 방위각 비대칭이 관측된다 (모델 기본값이 No Symmetry).
3. **CellValue_UI 인자 순서** (longitude, latitude): 기존 회전대칭 코드는 단일
   longitude 슬라이스만 읽었으므로 2D 확장 시 순서를 실측 확인.
4. 경로/파일명(`C:\Users\jhkim\...`)은 사용자 머신 기준. `.mat`(nk_JH33 등)과
   `SweptEntity.2.lts`, `Lens_size_effect_..._v1.1.lts`가 해당 경로에 있어야 함.

## 방향 제어성 증명 그림 (권장)

`target_phi`를 0:45:315로 스윕하며 각각 독립 최적화 → 각 φ에서 발광 피크가 실제로
그 방위각을 따라 도는지 보이면 "동일 렌즈의 회전정렬만으로 발광 방향 제어" 주장을
직접 입증. `target_theta` 스윕은 달성 가능한 최대 tilt 각도 경계를 보여준다.

## 계승 (기존 v4에서)

탐색/검증 ray 분리, 크래시 시 NaN 반환(GP 오염 방지), warm-start 상위 N점,
수렴 판정 적응예산, patternsearch 국소정련, 고정밀 반복검증 mean±std.

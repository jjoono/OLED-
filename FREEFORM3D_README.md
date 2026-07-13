# 비대칭 3D Freeform MLA 최적화 (BO 기반 trans-scale)

hemisphere/회전대칭 렌즈로는 **원리적으로 불가능한 방향성 비대칭 발광**을, 모든 렌즈가
동일한 비대칭 freeform 형상을 갖고 같은 방향으로 정렬된 MLA로 실현하기 위한 설계 코드.

## 파일 구성

| 파일 | 역할 |
|------|------|
| `freeform_grid_height.m` | **임의 3D** 렌즈 높이장 `z(x,y)= clamp(interp2_spline(C)) · (1-ρ²)^p`. NgxNg 제어높이 격자 `C` → 국소 bump/dimple 포함 임의 비대칭 형상. **현재 사용 파라미터화**. |
| `generate_freeform_mesh.m` | 높이장 함수핸들 → **watertight 삼각 메쉬** → ASCII STL. plano-convex(평평한 바닥) 닫힌 solid. 파라미터화 무관(핸들만 받음). |
| `BO_Freeform3D_asym.m` | 메인 드라이버. Bayesian Optimization으로 freeform 격자 + 마이크로캐비티 동시 최적화. 목적함수 = **목표 방향으로 조향 추출된 EQE (`EQE_cone`)**. **← 이 파일을 Run**. |
| `RenewLightTools.m` | LightTools 2개 인스턴스 재시작/연결 (메인·테스트 공용). 경로/사용자명 여기서 수정. |
| `test_freeform_geom.m` | **BO 전에 먼저 실행.** STL import 명령과 far-field mesh 차원(아래 @@VERIFY)만 빠르게 점검. |
| `freeform_height.m` | (대안) 회전대칭+harmonic 파라미터화 `z=H·P(ρ)·[1+Σρ^m(c_m cosmφ+s_m sinmφ)]`. 저차 "약한 tilt"용. 메쉬가 핸들 기반이라 교체 가능하나, 현재 메인은 격자 방식 사용. |

## 파라미터화 (19-dim; 형상 16 + 마이크로캐비티 3)

```
g11..g44    : 4x4 B-spline 제어높이 격자 [mm], 각 [0, 0.8]   ← 임의 3D 형상 (16 DOF)
dETL, dHTL  : OLED 마이크로캐비티 두께 (나노 스케일, EQE에 직접)
stretchZ    : 텍스처 z 스트레치
```

- 높이장 `z(x,y)`는 정의상 **single-valued** → 자기교차 불가 → 기하학적으로 항상 유효.
- rim(ρ=1)은 윈도우 `(1-ρ²)^p (p≥1)`로 0 고정 → 메쉬 watertight + 유한 rim 기울기.
- 형상 유효성(`isValidFreeform`): 최소 정점높이 `≥0.08mm`, 최대 기울기 `≤3.0`
  (제조/레이트레이싱 안정). 랜덤 시드 수용률 ~100%(검증됨).
- **DOF/optimizer**: 19 DOF는 GP-BO(`bayesopt`) 상한(~20) 내. `NGRID`를 5(→28 DOF)
  이상으로 키우려면 `surrogateopt`/CMA-ES 등으로 optimizer 교체 필요.
  (`BO_Freeform3D_asym.m` 상단 `NGRID`, `GRID_HMAX`, `GRID_WINPOW`로 조절)

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

## 실행 전 확인 필요 (`@@VERIFY`) — 먼저 `test_freeform_geom.m` 로 점검

이 코드는 **Windows + LightTools 2023.03 + MATLAB(COM)** 환경에서 실행된다. LightTools
명령/DB 이름은 버전마다 달라 아래 두 가지를 **본 BO 를 돌리기 전에** 확정해야 한다.
`test_freeform_geom.m` 이 이 두 가지만 빠르게 반복 점검하는 스크립트다.

### 1. STL import 명령 (`LT_IMPORT_CMD`) — 매크로 레코더로 확보 (필수)
초기 버전이 쓴 `ImportCADFile` 은 **존재하지 않는 명령**이라 import 가 실패했다
(콘솔: `Unknown variable "ImportCADFile"`, 이어서 RepairEntities/SaveLibrary 도 선택
지오메트리가 없어 실패). 올바른 명령은 버전마다 다르므로 직접 뽑아야 한다:

1. LightTools 에서 **매크로 레코더**(Console/Macro) 시작
2. **File > Import** 로 `.stl` 파일 하나를 수동 임포트
3. 레코더 정지 → 기록된 명령을 복사
4. 그 명령의 파일경로 자리를 `%s` 로 바꿔 **`BO_Freeform3D_asym.m` 상단
   `LT_IMPORT_CMD`** 에 넣고, 임포트된 solid 이름(레코더의 `DefaultSelect` 줄)을
   `LT_IMPORT_SOLIDNAME` 에 넣는다.
5. `test_freeform_geom.m` 실행 → `[OK] .ent 생성됨` 이 나오면 통과.

> `LTCmd` 는 잘못된 명령이어도 MATLAB 예외를 안 던지고 콘솔에만 오류를 찍는다. 그래서
> 본 코드는 **.ent 파일이 실제로 생겼는지**로 성공을 판정하고, 실패 시 그 평가를
> 조용히 넘기지 않고 NaN 처리한다.

### 2. far-field mesh 차원 (`MESH_NLONG`, `MESH_NLAT`) — 모델과 일치 (필수)
초기 버전이 `nLong=36` 으로 읽어 `잘못된 인덱스(31,·) ... CellValue UI` 에러가 났다.
= 그 메쉬의 longitude 축은 **30 셀**뿐인데 넘어선 것. 에러의 index `N` 에서 `N-1` 이
실제 longitude 셀 수다. 현재 기본값을 **longitude 30 × latitude 90** 으로 맞춰뒀고,
GUI 의 far-field `INTENSITY_MESH` 셀 수와 대조해 `BO_Freeform3D_asym.m` 상단
`MESH_NLONG/MESH_NLAT` (및 각도 범위 `MESH_LAT_MAX` 등)에 반영할 것.
`CellValue_UI(iLong, iLat)` 규약(첫 인자=longitude, 둘째=latitude)은 기존 코드와 동일.

> **[가장 중요] 방위각 비대칭이 실제로 관측되려면 이 far-field 메쉬의
> Symmetry = "No Symmetry" 여야 한다.** 회전대칭/미러 대칭이 켜져 있으면 렌즈를
> 비대칭으로 만들어도 메쉬가 방위각을 평균화해 EQE_cone 이 조향에 반응하지 않는다.
> (모델 기본 receiver default 는 No Symmetry 지만, 해당 mesh 개체에서 재확인할 것.)

### 3. 기타
- 경로/파일명(`C:\Users\jhkim\...`), 사용자명(`jhkim`, `RenewLightTools.m`)은 사용자
  머신 기준. `.mat`(nk_JH33 등), 두 `.lts`, CPS 관련 함수가 경로에 있어야 함.
- import 파이프라인이 통과하고 mesh 차원이 맞으면, 초기 로그의 매 평가 `Error`(NaN)는
  사라진다 (그 NaN 들은 import 실패 + mesh 인덱스 초과의 복합 결과였다).

## 방향 제어성 증명 그림 (권장)

`target_phi`를 0:45:315로 스윕하며 각각 독립 최적화 → 각 φ에서 발광 피크가 실제로
그 방위각을 따라 도는지 보이면 "동일 렌즈의 회전정렬만으로 발광 방향 제어" 주장을
직접 입증. `target_theta` 스윕은 달성 가능한 최대 tilt 각도 경계를 보여준다.

## 계승 (기존 v4에서)

탐색/검증 ray 분리, 크래시 시 NaN 반환(GP 오염 방지), warm-start 상위 N점,
수렴 판정 적응예산, patternsearch 국소정련, 고정밀 반복검증 mean±std.

% ============================================================
%  BO_Freeform3D_asym.m
%  비대칭(non-rotationally-symmetric) freeform MLA 최적화 — 방향성 발광
%
%  [논문 core] 모든 개별 렌즈가 "동일한 비대칭 형상"을 갖고 같은 방향으로 정렬된
%  MLA 는 방향성 비대칭 발광을 만든다. hemisphere/회전대칭 렌즈로는 원리적으로
%  불가능(azimuthal 대칭이 기하학적으로 강제됨). 이 코드는 Bayesian Optimization
%  으로 그 비대칭 freeform 프로파일을 trans-scale(나노 CPS + 마이크로 레이트레이싱)
%  시뮬레이션 하에 설계한다.
%
%  --- 기존 회전대칭 코드(PSO_test_..._rOLED_sweep_v4.m) 대비 변경점 ---
%   (1) 형상 표현: 7점 회전스윕 스플라인 -> 임의 3D freeform 높이장 z(x,y).
%       NgxNg B-spline 제어높이 격자(interp2 spline) + rim 윈도우 (1-rho^2)^p.
%       (freeform_grid_height.m). Ng=4 -> 16 형상 DOF: 국소 bump/dimple 포함
%       임의 3D 비대칭 형상 표현. 높이장이라 항상 single-valued -> 기하 유효.
%   (2) 지오메트리 주입: SetSweptProfilePoints(회전스윕) -> generate_freeform_mesh
%       로 watertight STL 생성 -> LightTools import -> SaveLibrary 로 .ent 변환
%       -> 기존 텍스처 unit-cell 파이프라인 재사용.
%   (3) 목적함수: 단일 방위각 슬라이스 -> 2D 원거리장(theta,phi) 전체를 읽어
%       "목표 방향(theta_t, phi_t) 원뿔로 추출된 EQE" (EQE_cone) 최대화.
%       추출효율 + 방향성 조향을 동시에 보상한다.
%   (4) 제약: isValidFreeform (최소높이 + 최대기울기; 제조/레이트레이싱 안정).
%
%   [DOF/optimizer] 형상 16 + 마이크로캐비티 3 = 19 DOF. GP-BO(bayesopt) 상한(~20)
%   내에서 동작. Ng 를 5(28 DOF) 이상으로 키우면 bayesopt 대신 surrogateopt/CMA-ES
%   등으로 교체 필요(README 참고).
%
%   [v4 로부터 계승] 탐색/검증 ray 분리, 크래시 시 NaN 반환, warm-start 상위 N점,
%   수렴 판정 적응예산, patternsearch 국소정련, 고정밀 반복검증 mean±std.
%
%   [!] ▶ ON-MACHINE 확인 필요 지점은 코드에서 '@@VERIFY' 로 표시했다.
%       (LightTools DB access 이름은 모델/버전에 따라 다를 수 있으므로 GUI 에서
%        Data Access Name 을 확인 후 맞춰야 한다.)
% ============================================================
clear;

%% ===== STL import 명령 (@@VERIFY: 매크로 기록으로 확보) — 최상단에서 사전 점검 =====
% LTCmd 로 SweptEntity 인스턴스에 freeform STL 을 임포트하는 명령. 버전마다 이름이
% 달라 추측이 무의미하므로, LightTools 매크로 레코더로 File>Import 를 한 번 수동
% 실행해 정확한 명령을 얻어 아래에 넣는다. 명령 안의 파일경로 자리에 '%s' 를 둘 것
% (sprintf 로 STL 경로가 치환된다). 예(형식 예시, 실제 명령 아님): 'ImportFile "%s"'
% import 후 생성되는 solid 이름을 LT_IMPORT_SOLIDNAME 에(레코더 DefaultSelect 줄에서 확인).
% (먼저 test_freeform_geom.m 로 이 두 값을 확정한 뒤 여기에 복사할 것.)
global LT_IMPORT_CMD LT_IMPORT_SOLIDNAME
LT_IMPORT_CMD       = '';              % <- 매크로 기록으로 얻은 명령 문자열 ('%s'=STL경로)
LT_IMPORT_SOLIDNAME = 'freeform_lens'; % <- import 후 solid 이름

% LightTools 를 재시작하기 전에(=사용자의 열린 세션을 죽이기 전에) 미리 중단.
if isempty(LT_IMPORT_CMD)
    error(['STL import 명령(LT_IMPORT_CMD)이 비어 있습니다. 먼저 test_freeform_geom.m 로 ' ...
        'File>Import 매크로를 기록해 정확한 명령을 확보한 뒤, 이 스크립트 상단 ' ...
        'LT_IMPORT_CMD/LT_IMPORT_SOLIDNAME 에 넣고 다시 실행하세요. (README 의 @@VERIFY 절)']);
end

%% ===== LightTools 연결 (v4 계승) =====
global ID_swept ID_LT ltml ltloc count eval_count restart_interval ray_nums_current
global target_theta target_phi cone_half   % 방향성 목적함수 타겟(전역: objFcn 접근)
RenewLightTools();
try
    ltml.LTCmd(ltml.GetLTAPI(ID_LT), 'Message "Check Connection"');
catch
    ltml  = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
end
count = 1;
restart_interval = 20;
lt = ltloc.GetLTAPI(ID_swept);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

%% ===== 정확도/시간 트레이드오프 (v4 계승) =====
RAY_SEARCH   = 20000;    % 탐색용 ray
RAY_FINAL    = 100000;   % 검증용 ray
N_FINAL_REP  = 3;        % 최종 반복 평가 -> mean±std
N_WARM       = 3;        % warm-start 상위점
N_FRESH      = 7;        % fresh valid 시드
POLISH_EVALS = 15;       % patternsearch 정련 예산
EXPLORATION_RATIO = 0.7;

INIT_EVAL_FIRST   = 60;
CONV_BLOCK        = 10;
CONV_TOL          = 0.002;
CONV_PATIENCE     = 2;
MAX_EVAL_PER_TGT  = 150;

%% ===== 방향성 발광 타겟 =====
% 목적함수 = 목표 방향 (target_theta, target_phi) 를 중심으로 반각 cone_half
% 원뿔 안으로 추출된 EQE (EQE_cone). phi 를 스윕하면 "같은 렌즈를 회전정렬만 해도
% 발광 방향이 따라 도는지"를, theta 를 스윕하면 "얼마나 크게 기울일 수 있는지"를
% 보여주는 논문 그림이 된다. 기본은 정면에서 30도 기울인 단일 방향 최적화.
target_theta = 30;    % [deg] 목표 극각(빔을 정면에서 이만큼 기울임)
target_phi   = 0;     % [deg] 목표 방위각(비대칭 정렬 방향)
cone_half    = 20;    % [deg] 원뿔 반각(집광 목표 범위)

%% ===== far-field INTENSITY_MESH 격자 사양 (@@VERIFY: 모델과 반드시 일치) =====
% [중요] 에러 "잘못된 인덱스(31,·) ... CellValue UI" 는 여기 값이 모델 메쉬의 실제
%  셀 수보다 크면 발생한다. 에러가 index 31 에서 났으므로 첫 축(longitude)은 30셀.
%  아래 값을 모델의 far-field 메쉬(기존 코드가 CellValue_UI 로 읽던 pos 3)의 실제
%  longitude/latitude 셀 수와 각도 범위에 정확히 맞출 것. 방위각 비대칭을 보려면
%  longitude 셀 수 > 1 이고 메쉬 Symmetry = No Symmetry 여야 한다.
global MESH_NLONG MESH_NLAT MESH_LONG_MIN MESH_LONG_MAX MESH_LAT_MIN MESH_LAT_MAX MESH_POS
MESH_POS      = 3;     % INTENSITY_MESH 리스트에서 far-field 메쉬 위치 (기존 코드와 동일)
MESH_NLONG    = 30;    % @@VERIFY: longitude(방위각) 셀 수  (에러 신호상 30)
MESH_NLAT     = 90;    % @@VERIFY: latitude(극각) 셀 수      (기존 코드가 90까지 읽음)
MESH_LONG_MIN = 0;     MESH_LONG_MAX = 360;   % 방위각 범위 [deg]
MESH_LAT_MIN  = 0;     MESH_LAT_MAX  = 90;    % 극각 범위 [deg] (전방 반구; 모델에 맞게)


% (선택) 방위각 스윕으로 방향 제어성 증명하려면 아래 리스트 사용:
%   phi_sweep = 0:45:315;  % 각 phi 마다 아래 최적화를 독립 실행
% 여기서는 단일 타겟 1회 최적화를 수행(스윕은 for 로 감싸면 됨).

%% ===== 최적화 변수 : 임의 3D freeform = NgxNg B-spline 제어격자 + 마이크로캐비티 =====
%   g11..gNN : Ng x Ng 제어높이 [mm] (임의 3D 형상; 높이장이라 항상 single-valued)
%   dETL,dHTL: 마이크로캐비티 두께 (나노 스케일, EQE 에 직접 영향)
%   stretchZ : 텍스처 z 스트레치
%
%   Ng=4 -> 16 형상 DOF + 3 = 19 DOF. (GP-BO 상한 ~20 내. Ng=5 로 올리면 28 DOF 가
%   되어 bayesopt 로는 버거우니 그때는 optimizer 교체 필요 - README 참고.)
global NGRID GRID_WINPOW GRID_HMAX
NGRID       = 4;      % 격자 한 변 (형상 DOF = NGRID^2)
GRID_WINPOW = 1.0;    % rim 윈도우 (1-rho^2)^winPow 지수 (>=1 : rim 기울기 유한)
GRID_HMAX   = 0.80;   % 제어높이 상한 [mm] (footprint 반경 1mm -> aspect 상한 ~0.8)

nGridDOF = NGRID^2;
gridNames = cell(1, nGridDOF);
for i = 1:NGRID
    for j = 1:NGRID
        gridNames{(i-1)*NGRID + j} = sprintf('g%d%d', i, j);
    end
end
varNames = [gridNames, {'dETL','dHTL','stretchZ'}];
lb = [zeros(1,nGridDOF)              , 10   10   0.1];
ub = [GRID_HMAX*ones(1,nGridDOF)     , 150  150  3.0];

optVars = optimizableVariable.empty(0, numel(lb));
for i = 1:numel(lb)
    optVars(i) = optimizableVariable(varNames{i}, [lb(i), ub(i)]);
end
fprintf('형상 DOF = %d (%dx%d 격자) + 3 = %d 총 DOF\n', nGridDOF, NGRID, NGRID, numel(lb));

%% ===== 초기 시드 + BO =====
eval_count = 0;
RenewLightTools();
lt = ltloc.GetLTAPI(ID_swept);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

initX = array2table(genValidFreeform(20, lb, ub), 'VariableNames', varNames);
ray_nums_current = RAY_SEARCH;

fprintf('\n######## Freeform BO 시작: target (theta=%d, phi=%d, cone=%d) ########\n', ...
    target_theta, target_phi, cone_half);

results = bayesopt(@bo_objective, optVars, ...
    'MaxObjectiveEvaluations', INIT_EVAL_FIRST, ...
    'XConstraintFcn',          @bo_xconstraint, ...
    'IsObjectiveDeterministic', false, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'ExplorationRatio',        EXPLORATION_RATIO, ...
    'Verbose', 1, 'PlotFcn', {}, 'InitialX', initX);

% --- 수렴 판정 적응예산 (v4 계승) ---
noImpCount = 0;
while noImpCount < CONV_PATIENCE && results.NumObjectiveEvaluations < MAX_EVAL_PER_TGT
    prevBestEst = -results.MinEstimatedObjective;
    addEval = min(CONV_BLOCK, MAX_EVAL_PER_TGT - results.NumObjectiveEvaluations);
    results = resume(results, 'MaxObjectiveEvaluations', addEval);
    newBestEst = -results.MinEstimatedObjective;
    relImp = (newBestEst - prevBestEst) / max(abs(newBestEst), eps);
    if relImp < CONV_TOL, noImpCount = noImpCount + 1; else, noImpCount = 0; end
    fprintf('[Converge] evals=%3d | bestEst EQE_cone=%.5f | relImp=%+.4f | noImp %d/%d\n', ...
        results.NumObjectiveEvaluations, newBestEst, relImp, noImpCount, CONV_PATIENCE);
end

%% ===== 국소 정련 (patternsearch) =====
bestX_BO = results.XAtMinEstimatedObjective;
x0 = table2array(bestX_BO);
fprintf('--- Local polish (patternsearch, %d evals) ---\n', POLISH_EVALS);
psOpts = optimoptions('patternsearch', ...
    'MaxFunctionEvaluations', POLISH_EVALS, 'InitialMeshSize', 0.1, ...
    'MeshTolerance', 1e-3, 'Cache', 'on', 'Display', 'iter');
try
    xPol = patternsearch(@polish_objective, x0, [],[],[],[], lb, ub, [], psOpts);
catch perr
    fprintf('[Warn] patternsearch 실패(%s). BO 결과만 사용.\n', perr.message);
    xPol = x0;
end

%% ===== 최종 고정밀 검증 =====
if isequal(xPol, x0), candX = {x0}; else, candX = {x0, xPol}; end
ray_nums_current = RAY_FINAL;
candMean = -inf(1, numel(candX));  candStd = zeros(1, numel(candX));
for c = 1:numel(candX)
    if ~isValidFreeform(candX{c}), continue; end
    e = nan(1, N_FINAL_REP);
    for rrep = 1:N_FINAL_REP, e(rrep) = simulate_dirEQE(candX{c}); end
    candMean(c) = mean(e, 'omitnan');  candStd(c) = std(e, 'omitnan');
    fprintf('  candidate %d : EQE_cone = %.5f ± %.5f (N=%d, %d rays)\n', ...
        c, candMean(c), candStd(c), N_FINAL_REP, RAY_FINAL);
end
[bestEQE, ci] = max(candMean);
if ~isfinite(bestEQE), bestEQE = -results.MinEstimatedObjective; ci = 1; end
bestX = array2table(candX{ci}, 'VariableNames', varNames);

%% ===== 결과 저장 + 형상/원거리장 리포트 =====
save('BO_Freeform_result.mat', 'bestX', 'bestEQE', 'results', ...
    'target_theta', 'target_phi', 'cone_half', 'varNames', 'lb', 'ub');
fprintf('\n######## Done: best EQE_cone = %.5f ± %.5f ########\n', bestEQE, candStd(ci));
disp(bestX);

% 최적 형상 미리보기(비대칭 tilt 방향 확인)
report_best_shape(table2array(bestX), varNames);


%% =====================================================================
%%  1회 평가 공용 래퍼 (주기 재시작 + 크래시 -> NaN) (v4 계승)
%% =====================================================================
function eqe = simulate_dirEQE(pt)
global ID_swept ltml ltloc eval_count restart_interval
eval_count = eval_count + 1;
if mod(eval_count, restart_interval) == 0
    fprintf('\n[Refresh] 시뮬 %d회. LightTools 재시작...\n', eval_count);
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept); ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    pause(2);
end
try
    eqe = objFcn_directionalEQE(pt).EQE_cone;
    if eqe == 0, eqe = NaN; end   % 기하오류 등: 실패로 처리(0 아님)
catch err
    fprintf('\n[Error] eval %d LightTools 충돌: %s\n', eval_count, err.message);
    eqe = NaN;
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept); ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
end
end

%% ===== bayesopt 목적함수 (EQE_cone 최대화 -> -EQE_cone 최소화) =====
function obj = bo_objective(Xtbl)
obj = -simulate_dirEQE(table2array(Xtbl));
end

%% ===== patternsearch 정련용 =====
function f = polish_objective(x)
if ~isValidFreeform(x), f = 0; return; end
e = simulate_dirEQE(x);
if isnan(e), e = 0; end
f = -e;
end

%% ===== bayesopt 제약: valid freeform 만 통과 =====
function tf = bo_xconstraint(Xtbl)
pts = table2array(Xtbl);
n = size(pts,1);  tf = false(n,1);
for k = 1:n, tf(k) = isValidFreeform(pts(k,:)); end
end

%% ===== warm-start 상위 K점 (v4 계승) =====
function T = topKPoints(results, K)
X = results.XTrace; f = results.ObjectiveTrace;
ok = isfinite(f); X = X(ok,:); f = f(ok);
[~, idx] = sort(f, 'ascend');  K = min(K, numel(idx));
T = X(idx(1:K), :);
end

%% ===== 무작위 valid freeform 시드 (rejection sampling; 격자는 수용률 ~100%) =====
function P = genValidFreeform(K, lb, ub)
dim = numel(lb);  P = zeros(K, dim);
for i = 1:K
    ok = false; tries = 0;
    while ~ok
        p = lb + rand(1, dim) .* (ub - lb);
        if isValidFreeform(p), ok = true; P(i,:) = p; end
        tries = tries + 1;
        if tries > 5000, error('genValidFreeform:seed', 'valid 시드 생성 실패 - 제약이 너무 빡빡함'); end
    end
end
end


%% =====================================================================
%%  Objective: freeform 지오메트리 주입 + 나노 CPS + 2D 방향성 EQE
%%  (나노 스케일 CPS/코팅 블록은 검증된 objFcn_angularEQE 에서 이식, 형상/각도
%%   읽기만 신규)
%% =====================================================================
function output = objFcn_directionalEQE(point)
global ID_LT ID_swept ltml ltloc count ray_nums_current
global target_theta target_phi cone_half MESH_POS

lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

% --- 파라미터 언팩 ---
global NGRID
C        = reshape(point(1:NGRID^2), [NGRID, NGRID]);  % 제어높이 격자
dETL     = point(NGRID^2 + 1);
dHTL     = point(NGRID^2 + 2);
stretchZ = point(NGRID^2 + 3);

% --- 고정 설정(기존과 동일 스케일) ---
d_sub = 1.3;  r_OLED = 1;  r_pat = 1;   % r_pat: 텍스처 패턴 크기(필요시 상수/변수화)
Lensheight = 0.01;
wavelength_start = 580;  wavelength_end = 590;  n = 10;
if isempty(ray_nums_current), ray_nums = 50000; else, ray_nums = ray_nums_current; end

% --- 시뮬 설정 (v4 계승) ---
List = ltml.LTDbList(lt,'lens_manager[1]','SIMULATIONS');
Key  = ltml.LTListByName(lt,List,'ForwardAll');
ltml.LTDbSet(lt,Key,'MaxProgress',ray_nums);
List = ltml.LTDbList(lt,'lens_manager[1]','CUBE_PRIMITIVE');
Key  = ltml.LTListByName(lt,List,'Substrate');
ltml.LTDbSet(lt,Key,'Height',d_sub);  ltml.LTDbSet(lt,Key,'Y',d_sub/2);
SRList = ltml.LTDbList(lt,'lens_manager[1]','CUBE_PRIMITIVE');
SRKey  = ltml.LTListAtPos(lt,SRList,2);
ltml.LTDbSet(lt,SRKey,'Y',d_sub+Lensheight/2);
List = ltml.LTDbList(lt,'lens_manager[1]','TEXTURE_ZONE_EXTENT');
Key  = ltml.LTListByName(lt,List,'zone');
ltml.LTDbSet(lt,Key,'Geometry_1',r_pat);  ltml.LTDbSet(lt,Key,'Geometry_2',r_pat);
List = ltml.LTDbList(lt,'lens_manager[1]','DISK_SOURCE');
Key  = ltml.LTListByName(lt,List,'DiskSource_18');
ltml.LTDbSet(lt,Key,'Radius',r_OLED);

% --- (신규) 임의 3D freeform 지오메트리 주입 ---
% 형상 유효성(최소높이/기울기)은 objective 진입 전 isValidFreeform 이 보장하지만
% 방어적으로 재확인.
if ~isValidFreeform(point)
    output = zero_output(); return;
end
entPath = updateFreeformGeometry(C, stretchZ);
if isempty(entPath)
    output = zero_output(); return;   % 지오메트리 생성/import 실패
end

% --- 나노 스케일 CPS + 하단 반사율 + 코팅 (v4 objFcn 과 동일) ---
load('nk_JH33.mat');  load('Photopic_400_800.mat');  load('CIE_1931.mat');  load('R_pd.mat');
wavelength = (wavelength_start:wavelength_end).';
wavelength_num = length(wavelength);
emission_spectrum = spectrum.l_I_Irdmppyph2tmd(wavelength_start-399:wavelength_end-399,:);
eta_rad = 0.98;  horizontal_dipole_ratio = 0.865;
bottom_air_refractive_index = ones(wavelength_num,1);

no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
layer_num = size(no_bar,2);
sin089 = sind(0:89);  cos089 = cosd(0:89);
no_bar = no_bar(wavelength_start-399:wavelength_end-399,:);
ne_bar = ne_bar(wavelength_start-399:wavelength_end-399,:);
thickness = [100 dETL 25 10 dHTL 150];
EML_position = 4;  z0 = 12.5;  u_data_num = 499;  max_u = 3;

CPS_result = CPS_for_Isub(no_bar,ne_bar,thickness,emission_spectrum,eta_rad, ...
    horizontal_dipole_ratio,bottom_air_refractive_index,EML_position,z0,u_data_num,max_u,wavelength);
EQE_sub_CPS = CPS_result.EQE_sub;

% 하단 반사율
TMF_p = TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_s = TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
R_p = abs(TMF_p.r_p).^2;  R_s = abs(TMF_s.r_s).^2;
Reflectance = (R_p + R_s)/2;

% .coa 코팅 파일 write (v4 버그수정: 읽기 전 fclose)
fileID = fopen(sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count), 'w');
fprintf(fileID,'%s\n%s%d\n%s\n%s\n%s\n%s\n ','DFAT Version 1.0','DATANAME: R_Bottom_',count,'ABSORBING: YES','INDEX: 1.51','DATAITEMS: TAVG RAVG');
for i = wavelength_start:wavelength_end
    fprintf(fileID,'%s  %d\n','wv',i);
    for j = 0:89
        fprintf(fileID,'%s  %d  %d  %.3f\n','AOI',j,0,Reflectance(i-wavelength_start+1,j+1));
    end
end
fclose(fileID);
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings]" LoadFileName="' sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count) '"']);
List = ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');  Key = ltml.LTListByName(lt,List,'R_Al');
List = ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');  Key = ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName',sprintf('R_Bottom_%d', count));

% 파장별 가중/소스 세팅 (v4 계승)
I_white = 0.5*(CPS_result.I_sub_s + CPS_result.I_sub_p);
P_white = I_white .* repmat(sin089, wavelength_num, 1);
weight_factor = sum(P_white, 2);

% --- (신규) 2D 원거리장 방향성 EQE ---
K = (wavelength_num-1)/n + 1;
Power_output   = zeros(1, wavelength_num);
% 목표 원뿔로 들어가는 파워 분율을 파장 샘플마다 저장
coneFrac       = zeros(1, wavelength_num);

for wv = 1:n:wavelength_num
    % 각도별 소스 세팅 (v4 계승)
    fileID = fopen('C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','w');
    fprintf(fileID,'%s  %d  %d  %d  %d  %d  %d','SPHEREMESH:',1,90,0,0,360,90);
    writematrix(flip(I_white(wv,:).'),'C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','Delimiter','tab','WriteMode','append');
    fclose(fileID);
    SRList = ltml.LTDbList(lt,'Lens_manager[1]','DISK_SOURCE');  SRKey = ltml.LTListAtPos(lt,SRList,1);
    ltml.LTDbSet(lt,SRKey,'Radiant_Power', weight_factor(wv));
    SRList = ltml.LTDbList(lt,'Lens_manager[1]','Spectral_region');  SRKey = ltml.LTListAtPos(lt,SRList,2);
    ltml.LTDbSet(lt,SRKey,'Spectral_Definition','Monochromatic');
    ltml.LTDbSet(lt,SRKey,'Single_Wavelength', wv+wavelength_start-1);
    List = ltml.LTDbList(lt,'lens_manager[1]','DIRECTION_GRID_APODIZER');  Key = ltml.LTListAtPos(lt,List,1);
    ltml.LTDbSet(lt,Key,'LoadFileName','C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt');

    % 시뮬레이션
    ltml.LTBegin(lt);  ltml.LTCmd(lt,'\V3D BeginAllSimulations');  ltml.LTEnd(lt);

    % 전체 파워 (mesh 1)
    List = ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');  Key = ltml.LTListAtPos(lt,List,1);
    Power_output(wv) = ltml.LTDbGet(lt,Key,'TotalPower');

    % --- 2D 원거리장 읽기: 방위각(longitude) x 극각(latitude) 그리드 ---
    % @@VERIFY: 아래 mesh index/차원/DB이름은 모델의 far-field INTENSITY_MESH 설정에
    %  맞춰야 한다. 모델 기본: longitude 0..360, latitude 0..180, No Symmetry.
    %  방위각 비대칭이 나타나려면 (a) longitude 셀 수 > 1, (b) Symmetry=None 필수.
    Key = ltml.LTListAtPos(lt,List,MESH_POS);   % far-field 방향성 mesh
    [Igrid, thC, phC] = read_intensity_grid(ltml, lt, Key);  % Igrid: [nLat x nLong]
    coneFrac(wv) = cone_power_fraction(Igrid, thC, phC, target_theta, target_phi, cone_half);
end

% 파장 샘플 축약 (v4 계승)
weight_factor_2 = zeros(K,1);  Power_output_2 = zeros(K,1);
EQE_sub_matrix_2 = zeros(K,1);  coneFrac_2 = zeros(K,1);
for k = 1:K
    idx = n*(k-1) + 1;
    weight_factor_2(k)  = weight_factor(idx);
    Power_output_2(k)   = Power_output(idx);
    EQE_sub_matrix_2(k) = CPS_result.EQE_sub_matrix(idx);
    coneFrac_2(k)       = coneFrac(idx);
end
EQE_wv_matrix = Power_output_2 ./ weight_factor_2;
EQE_sub_matrix_2 = EQE_sub_matrix_2 / sum(EQE_sub_matrix_2) * EQE_sub_CPS;

% 총 추출 EQE 및 목표 원뿔로 조향된 EQE
EQE_total = sum(EQE_wv_matrix .* EQE_sub_matrix_2);
EQE_cone  = sum(EQE_wv_matrix .* EQE_sub_matrix_2 .* coneFrac_2);   % <- 목적함수

output = struct();
output.EQE_total = EQE_total;
output.EQE_cone  = EQE_cone;
output.coneFrac  = EQE_cone / max(EQE_total, eps);   % 목표 원뿔로 조향된 비율(방향 선택성)

% 코팅 정리 (v4 계승)
List = ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');  Key = ltml.LTListByName(lt,List,'R_Al');
List = ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');  Key = ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName','R_temp');
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings].COATING[' sprintf('R_Bottom_%d', count) ']" Delete= \Q']);
fclose('all');
end

%% ===== 실패 시 0 출력 =====
function output = zero_output()
output = struct('EQE_total',0,'EQE_cone',0,'coneFrac',0);
end


%% =====================================================================
%%  (신규) 비대칭 freeform -> STL -> LightTools import -> .ent
%% =====================================================================
function entPath = updateFreeformGeometry(C, stretchZ)
% 임의 3D freeform 렌즈를 STL 로 만들어 SweptEntity LightTools 인스턴스에 import
% 하고, 기존 파이프라인이 기대하는 .ent(ACIS) 로 SaveLibrary 하여 텍스처 unit-cell
% 파일을 교체한다. STL->.ent 변환(ACIS kernel)은 LightTools 가 담당하므로 MATLAB
% 이 직접 ACIS 를 쓸 필요가 없다.
global ID_swept ID_LT ltml ltloc GRID_WINPOW
entPath = '';
lt  = ltloc.GetLTAPI(ID_swept);   % 렌즈 마스터 인스턴스
lt2 = ltloc.GetLTAPI(ID_LT);      % 배열 시뮬 인스턴스

% 랜덤 파일명 (rng 오염 방지 위해 datetime 기반)
tagc = char(datetime('now','Format','yyyyMMddHHmmssSSS'));
base = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\';
stlPath = [base 'freeform_' tagc '.stl'];
entPath_local = [base 'swept_' tagc '.ent'];
entPath_mod   = [base 'swept_' tagc '.1.ent'];   % Repair 후 텍스처가 참조하는 이름

% 1) STL 생성 (footprint 반경 = 렌즈 기존 스케일 1 mm)
%    파라미터화 무관 mesh: 높이장 핸들만 넘긴다.
winPow = GRID_WINPOW;
hfun   = @(X,Y) freeform_grid_height(X, Y, C, winPow);
mopts  = struct('nr',40,'nt',120,'Rfoot',1,'solidName','freeform_lens');
try
    minfo = generate_freeform_mesh(hfun, stlPath, mopts); %#ok<NASGU>
catch me
    fprintf('[Geom] 메쉬 생성 실패: %s\n', me.message);  return;
end

% 2) LightTools 로 import + Repair + .ent 저장
%    [중요] LTCmd 는 잘못된 명령이어도 MATLAB 예외를 던지지 않고 LightTools 콘솔에만
%    "Unknown variable ..." 를 찍는다. 따라서 성공 여부는 .ent 파일이 실제로
%    생성됐는지로 판정한다(아래). import 명령 문자열은 버전마다 다르므로 전역
%    LT_IMPORT_CMD 로 분리했다 -> 매크로 기록으로 얻은 정확한 명령을 넣을 것.
global LT_IMPORT_CMD LT_IMPORT_SOLIDNAME
if isempty(LT_IMPORT_CMD)
    error(['[Geom] STL import 명령이 설정되지 않았습니다. LightTools 매크로 레코더로 ' ...
        'File>Import 로 .stl 하나를 수동 임포트해 정확한 명령을 얻은 뒤, 스크립트 상단 ' ...
        'LT_IMPORT_CMD 에 넣으세요. (자세한 방법: FREEFORM3D_README.md 의 STL import 절)']);
end
if exist(entPath_local, 'file'), delete(entPath_local); end   % 이전 잔재 제거(판정 신뢰성)
try
    ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    % LT_IMPORT_CMD 안의 '%s' 자리에 STL 경로가 들어간다 (sprintf).
    ltml.LTCmd(lt, sprintf(LT_IMPORT_CMD, stlPath));
    ltml.LTCmd(lt, sprintf('DefaultSelect "%s"', LT_IMPORT_SOLIDNAME));  % 임포트된 solid 선택
    ltml.LTCmd(lt, 'RepairEntities');
    ltml.LTCmd(lt, sprintf('SaveLibrary XYZ 0,0,0 "%s"', entPath_local));
catch me
    fprintf('[Geom] LightTools import/저장 예외: %s\n', me.message);  return;
end
pause(0.3);
if ~exist(entPath_local, 'file')
    fprintf(['[Geom] .ent 저장 실패 -> import 명령이 안 먹혔을 가능성. ' ...
        'LightTools 콘솔의 오류 메시지와 LT_IMPORT_CMD/LT_IMPORT_SOLIDNAME 를 확인하세요.\n']);
    return;   % 조용히 진행하지 않고 이 평가를 실패(NaN) 처리
end

% 2b) 임포트된 solid 정리 — 기존 코드는 SetSweptProfilePoints 로 "한 개"를 수정했지만
%     지금은 매 평가마다 새 solid 를 import 하므로, 저장 후 삭제하지 않으면 마스터
%     모델에 계속 쌓여 이름 충돌/메모리 문제가 생긴다. 선택 후 삭제로 클린 상태 유지.
%     @@VERIFY: 삭제 명령/선택자는 매크로 기록으로 확인(레코더에서 solid 삭제 시 기록됨).
try
    ltml.LTCmd(lt, sprintf('DefaultSelect "%s"', LT_IMPORT_SOLIDNAME));
    ltml.LTCmd(lt, 'Delete');
catch
    % 삭제 실패는 치명적이지 않으나, 반복 누적을 막으려면 restart_interval 을 짧게.
end

% 3) 배열 모델의 텍스처 unit-cell 파일 + StretchZ 갱신 (v4 계승)
try
    List = ltml.LTDbList(lt2,'LENS_MANAGER[1]','LIBRARY_ELEMENT_UNIT_CELL');
    Key  = ltml.LTListByName(lt2,List,'LibraryElement');
    ltml.LTDbSet(lt2,Key,'Filename', entPath_mod);
    List = ltml.LTDbList(lt2,'LENS_MANAGER[1]','TEXTURE_PARAMETER');
    Key  = ltml.LTListByName(lt2,List,'StretchZ');
    ltml.LTDbSet(lt2,Key,'Value', stretchZ);
catch me
    fprintf('[Geom] 텍스처 갱신 실패: %s\n', me.message);  return;
end

entPath = entPath_mod;
end


%% =====================================================================
%%  (신규) 2D 인텐시티 그리드 읽기 + 원뿔 파워 분율
%% =====================================================================
function [Igrid, thetaC, phiC] = read_intensity_grid(ltml, lt, Key)
% far-field INTENSITY_MESH 를 [nLat x nLong] 격자로 읽는다.
% 차원/각도범위는 상단 전역(MESH_*)에서 가져오며 모델 메쉬와 정확히 일치해야 한다.
% CellValue_UI(iLong, iLat) 규약(기존 코드가 CellValue_UI(1,lat) 로 읽던 것과 동일:
% 첫 인자=longitude, 둘째 인자=latitude).
% [비용 주의] 셀별 LTDbGet 는 COM 호출이므로 nLong*nLat 회/파장/평가로 급증한다.
%  30x90=2700회. 느리면 상단에서 해상도를 줄이거나(모델 메쉬도 함께 축소),
%  LightTools 메쉬 데이터 일괄 export 커맨드로 대체(매크로 기록으로 확인).
global MESH_NLONG MESH_NLAT MESH_LONG_MIN MESH_LONG_MAX MESH_LAT_MIN MESH_LAT_MAX
nLong = MESH_NLONG;  nLat = MESH_NLAT;

Igrid = zeros(nLat, nLong);
for iL = 1:nLong
    for iT = 1:nLat
        Igrid(iT, iL) = ltml.LTDbGet(lt, Key, 'CellValue_UI', iL, iT);
    end
end
% 셀 중심 각도 [deg]
thetaC = MESH_LAT_MIN  + (MESH_LAT_MAX -MESH_LAT_MIN) /nLat  * ((1:nLat)  - 0.5);  % 1xnLat
phiC   = MESH_LONG_MIN + (MESH_LONG_MAX-MESH_LONG_MIN)/nLong * ((1:nLong) - 0.5);  % 1xnLong
end

function frac = cone_power_fraction(Igrid, thetaC, phiC, th_t, ph_t, halfAng)
% 목표 방향(th_t, ph_t) 중심 반각 halfAng 원뿔로 방출되는 파워 분율.
%   dP ∝ I(theta,phi) * sin(theta) dtheta dphi
%   원뿔 판정: 목표 방향 단위벡터와 셀 방향 단위벡터의 각거리 <= halfAng
[TH, PH] = ndgrid(thetaC, phiC);          % [nLat x nLong], deg
sinth = sind(TH);
W = Igrid .* sinth;                        % 파워 가중(상수는 분율에서 상쇄)
Wtot = sum(W(:));
if Wtot <= 0, frac = 0; return; end

% 목표 방향 벡터
d_t = ang2vec(th_t, ph_t);
% 셀 방향과의 각거리
V = ang2vec(TH, PH);                        % 3 x (nLat*nLong) 형태로 계산
cosd_sep = V(1,:).*d_t(1) + V(2,:).*d_t(2) + V(3,:).*d_t(3);
cosd_sep = reshape(cosd_sep, size(TH));
inCone = (cosd_sep >= cosd(halfAng));

frac = sum(W(inCone)) / Wtot;
end

function V = ang2vec(theta, phi)
% 구면각(theta=극각, phi=방위각, deg) -> 단위벡터. 배열 입력 시 3 x N.
th = theta(:).'; ph = phi(:).';
V = [sind(th).*cosd(ph); sind(th).*sind(ph); cosd(th)];
if isscalar(theta), V = V(:); end
end


%% =====================================================================
%%  (신규) freeform 형상 유효성
%% =====================================================================
function TF = isValidFreeform(X)
% X: 1xD (또는 NxD). 격자 freeform 파라미터가 물리적/제조적으로 타당한지 검사.
%   높이장 z(x,y) 는 정의상 항상 single-valued 이므로 자기교차는 불가능하다.
%   실제 제약은 (1) 형상이 유의미(너무 납작하지 않음), (2) 표면 기울기가 제조/
%   레이트레이싱 가능 범위 이내.  (수용률 ~100% 설계 - README/프로토타입 검증)
global NGRID GRID_WINPOW
persistent XX YY
if isempty(XX)
    t = linspace(-1, 1, 61);   % 유효성 평가용 저해상 격자(속도)
    [XX, YY] = meshgrid(t, t);
end
SLOPE_MAX = 3.0;    % 최대 |grad z| (제조/레이트레이싱 안정)
Z_MIN     = 0.08;   % 최소 정점 높이 [mm] (degenerate 방지)

numRows = size(X,1);
TF = true(numRows,1);
nGrid = NGRID^2;
t = linspace(-1, 1, 61);
for k = 1:numRows
    C = reshape(X(k, 1:nGrid), [NGRID, NGRID]);
    Z = freeform_grid_height(XX, YY, C, GRID_WINPOW);
    Zf = Z;  Zf(~isfinite(Zf)) = 0;

    bad = false;
    if max(Zf(:)) < Z_MIN, bad = true; end
    if ~bad
        [gx, gy] = gradient(Zf, t, t);
        if max(hypot(gx(:), gy(:))) > SLOPE_MAX, bad = true; end
    end
    TF(k) = ~bad;
end
end


%% =====================================================================
%%  형상 리포트 (최적 결과 시각화; LightTools 불필요, 순수 MATLAB)
%% =====================================================================
function report_best_shape(x, varNames) %#ok<INUSD>
global NGRID GRID_WINPOW
C = reshape(x(1:NGRID^2), [NGRID, NGRID]);
t = linspace(-1,1,200);  [Xg,Yg] = meshgrid(t,t);
Z = freeform_grid_height(Xg, Yg, C, GRID_WINPOW);
m = Z; m(~isfinite(m)) = 0;
cx = sum(Xg(:).*m(:))/sum(m(:));  cy = sum(Yg(:).*m(:))/sum(m(:));
figure('Name','Best 3D freeform (grid)');
subplot(1,2,1); surf(Xg,Yg,Z,'EdgeColor','none'); axis tight; view(35,30);
title(sprintf('Best freeform  z_{max}=%.3f mm', max(m(:)))); xlabel x; ylabel y; zlabel z;
subplot(1,2,2); contourf(Xg,Yg,Z,25,'LineColor','none'); axis equal tight; hold on;
plot(cx,cy,'r+','MarkerSize',14,'LineWidth',2); plot(0,0,'w.','MarkerSize',8);
title(sprintf('height centroid=(%+.3f,%+.3f)  [조향 \\rightarrow +x=phi 0]', cx, cy));
fprintf('[Shape] height centroid (비대칭/조향 방향 지표) = (%+.4f, %+.4f)\n', cx, cy);
end


%% RenewLightTools 는 별도 파일(RenewLightTools.m)로 분리됨 (test 스크립트와 공용).

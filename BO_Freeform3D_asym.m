% ============================================================
%  BO_Freeform3D_asym.m
%  비대칭(non-rotationally-symmetric) freeform MLA 최적화 — 방향성 발광
%  [지오메트리 방식: LightTools 네이티브 Freeform LensElement (.ent 직접 생성)]
%
%  [논문 core] 모든 개별 렌즈가 "동일한 비대칭 형상"을 갖고 같은 방향으로 정렬된
%  MLA 는 방향성 비대칭 발광을 만든다. hemisphere/회전대칭 렌즈로는 원리적으로
%  불가능. 이 코드는 Bayesian Optimization 으로 그 비대칭 freeform 을 trans-scale
%  (나노 CPS + 마이크로 레이트레이싱) 시뮬 하에 설계한다.
%
%  --- 지오메트리 파이프라인 (STL/import/SaveLibrary 전부 제거) ---
%   LightTools FreeformEntity(.ent) 는 표면 격자점을 X,Y,Z 텍스트로 담는다.
%   따라서 MATLAB 이 검증된 템플릿(.ent)에서 FrontSurface 격자의 Z 만 바꿔 .ent 를
%   직접 쓰고, 배열 모델 텍스처의 LibraryElement Filename 을 그 .ent 로 지정한다
%   (기존 코드가 swept_XXX.1.ent 를 물리던 자리와 동일). import 불필요.
%
%   * 템플릿은 사용자가 GUI 에서 만든 "유효 solid" freeform 렌즈:
%     Front/Rear 모두 규칙 5x5 격자, 두께 1mm, 경계 Z=0, NURBS Off(+SmoothResample).
%   * BO 변수 = FrontSurface 내부점 Z (5x5 -> 내부 3x3 = 9개). 경계 16점은 Z=0 고정
%     -> bump 가 기판에 매끄럽게 맞물림. + 마이크로캐비티 dETL,dHTL,stretchZ.
%   * [NURBS 주의] 템플릿은 NURBS Off. 켜면 LightTools 가 U,V 를 바꿔 면을 재구성해
%     격자가 흔들리므로, 생성기는 템플릿의 restoreNURBS:"No" 를 그대로 유지한다.
%
%  --- 목적함수 ---
%   2D 원거리장(theta,phi)을 읽어 "목표 방향(theta_t,phi_t) 원뿔로 추출된 EQE"
%   (EQE_cone) 최대화. 추출효율 + 방향성 조향 동시 보상. 대칭 발광은 어떤 방향으로도
%   대칭이라 EQE_cone 을 못 올림 -> 비대칭 freeform 만 조향 가능(= 논문 판별 지표).
%
%  [v4 계승] 탐색/검증 ray 분리, 크래시 시 NaN, warm-start, 수렴 적응예산,
%  patternsearch 정련, 고정밀 반복검증 mean±std.
%
%  [!] @@VERIFY: 경로/파일명, far-field mesh 차원(MESH_*), 텍스처 LibraryElement
%      세팅은 사용자 모델에 맞춰 확인. far-field mesh 는 Symmetry=No Symmetry 필수.
% ============================================================
clear;

%% ===== 경로/템플릿 설정 (@@VERIFY: 사용자 머신 경로) =====
global FF_TEMPLATE FF_XY FF_INNER FF_N FF_BASE
FF_BASE     = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\';
FF_TEMPLATE = [FF_BASE 'freeform_template.ent'];   % 검증된 유효 solid 템플릿을 이 경로에 둘 것
[FF_XY, FF_INNER, FF_N] = freeform_grid_info(FF_TEMPLATE);
nInner = numel(FF_INNER);
fprintf('Freeform 격자: 총 %d 점, 내부(자유) %d 점 -> 형상 DOF=%d\n', FF_N, nInner, nInner);

%% ===== LightTools 연결 (배열 모델 1개) =====
global ID_LT ltml ltloc count eval_count restart_interval ray_nums_current
global target_theta target_phi cone_half MESH_POS
global MESH_NLONG MESH_NLAT MESH_LONG_MIN MESH_LONG_MAX MESH_LAT_MIN MESH_LAT_MAX
RenewLightTools();
try
    ltml.LTCmd(ltml.GetLTAPI(ID_LT), 'Message "Check Connection"');
catch
    ltml  = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
end
count = 1;
restart_interval = 20;
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

%% ===== 정확도/시간 트레이드오프 (v4 계승) =====
RAY_SEARCH   = 20000;    RAY_FINAL = 100000;   N_FINAL_REP = 3;
N_WARM       = 3;        POLISH_EVALS = 15;    EXPLORATION_RATIO = 0.7;
INIT_EVAL_FIRST = 50;    CONV_BLOCK = 10;      CONV_TOL = 0.002;
CONV_PATIENCE   = 2;     MAX_EVAL_PER_TGT = 120;

%% ===== 방향성 발광 타겟 =====
% 목적함수 J = EQE_cone_fwd(theta_t,phi_t) - W_CONTRAST * EQE_cone_opp(theta_t,phi_t+180)
%   - EQE_cone_fwd : 목표 원뿔로 추출된 EQE (HUD eyebox 로 가는 유용 효율)
%   - EQE_cone_opp : 반대 방위각 같은 극각 원뿔로 새는 EQE (손실/글레어)
%   hemisphere(방위각 대칭)는 fwd=opp -> 대비항이 이득을 상쇄 -> 구조적으로 불리.
%   비대칭 freeform 만 fwd>>opp 로 J 를 키운다(= 논문 판별 지표).
global W_CONTRAST
target_theta = 30;    % [deg] 목표 극각(빔을 정면에서 이만큼 기울임)
target_phi   = 0;     % [deg] 목표 방위각(= 렌즈 정렬 방향; 배열 회전으로 임의 제어)
cone_half    = 20;    % [deg] 원뿔 반각(eyebox 크기)
W_CONTRAST   = 0.5;   % 반대편 누설 페널티 가중(0=순수 EQE_cone, 1=순수 대비차)

%% ===== far-field INTENSITY_MESH 격자 사양 (@@VERIFY: 모델과 일치) =====
% "잘못된 인덱스(N,·) CellValue UI" 에러 시 N-1 이 실제 longitude 셀 수.
% 방위각 비대칭을 보려면 longitude 셀>1 이고 mesh Symmetry=No Symmetry 여야 함.
MESH_POS      = 3;    MESH_NLONG = 30;   MESH_NLAT = 90;
MESH_LONG_MIN = 0;    MESH_LONG_MAX = 360;
MESH_LAT_MIN  = 0;    MESH_LAT_MAX  = 90;

%% ===== 최적화 변수 : FrontSurface 내부 Z + 마이크로캐비티 =====
%   z1..z{nInner} : FrontSurface 내부 격자점 높이 [정규 단위]. 경계는 0 고정.
%   dETL,dHTL     : OLED 마이크로캐비티 두께 (나노 스케일)
%   stretchZ      : 텍스처 z 스트레치
ZLO = 0.0;  ZHI = 0.8;   % 내부점 Z 범위 (두께 1mm 보다 작게 -> solid 유효 유지)
zNames = arrayfun(@(k) sprintf('z%d',k), 1:nInner, 'UniformOutput', false);
varNames = [zNames, {'dETL','dHTL','stretchZ'}];
lb = [ZLO*ones(1,nInner),  10  10  0.1];
ub = [ZHI*ones(1,nInner),  150 150 3.0];

optVars = optimizableVariable.empty(0, numel(lb));
for i = 1:numel(lb)
    optVars(i) = optimizableVariable(varNames{i}, [lb(i), ub(i)]);
end
fprintf('총 DOF = %d (형상 %d + 마이크로캐비티 3)\n', numel(lb), nInner);

%% ===== 초기 시드 + BO =====
eval_count = 0;
initX = array2table(lb + rand(20, numel(lb)) .* (ub - lb), 'VariableNames', varNames);
ray_nums_current = RAY_SEARCH;

fprintf('\n######## Freeform BO 시작: target (theta=%d, phi=%d, cone=%d) ########\n', ...
    target_theta, target_phi, cone_half);

results = bayesopt(@bo_objective, optVars, ...
    'MaxObjectiveEvaluations', INIT_EVAL_FIRST, ...
    'IsObjectiveDeterministic', false, ...
    'AcquisitionFunctionName', 'expected-improvement-plus', ...
    'ExplorationRatio', EXPLORATION_RATIO, ...
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
x0 = table2array(results.XAtMinEstimatedObjective);
fprintf('--- Local polish (patternsearch, %d evals) ---\n', POLISH_EVALS);
psOpts = optimoptions('patternsearch', 'MaxFunctionEvaluations', POLISH_EVALS, ...
    'InitialMeshSize', 0.05, 'MeshTolerance', 1e-3, 'Cache', 'on', 'Display', 'iter');
try
    xPol = patternsearch(@polish_objective, x0, [],[],[],[], lb, ub, [], psOpts);
catch perr
    fprintf('[Warn] patternsearch 실패(%s). BO 결과만 사용.\n', perr.message);  xPol = x0;
end

%% ===== 최종 고정밀 검증 =====
if isequal(xPol, x0), candX = {x0}; else, candX = {x0, xPol}; end
ray_nums_current = RAY_FINAL;
candMean = -inf(1, numel(candX));  candStd = zeros(1, numel(candX));
for c = 1:numel(candX)
    e = nan(1, N_FINAL_REP);
    for rrep = 1:N_FINAL_REP, e(rrep) = simulate_dirEQE(candX{c}); end
    candMean(c) = mean(e, 'omitnan');  candStd(c) = std(e, 'omitnan');
    fprintf('  candidate %d : EQE_cone = %.5f ± %.5f (%d rays)\n', c, candMean(c), candStd(c), RAY_FINAL);
end
[bestEQE, ci] = max(candMean);
if ~isfinite(bestEQE), bestEQE = -results.MinEstimatedObjective; ci = 1; end
bestX = array2table(candX{ci}, 'VariableNames', varNames);

%% ===== 결과 저장 + 형상 리포트 =====
% 최적 형상의 방향성 분해(목표 원뿔 vs 반대편) 재평가 (고정밀)
ray_nums_current = RAY_FINAL;
try
    bd = objFcn_directionalEQE(table2array(bestX));
catch
    bd = struct('EQE_total',NaN,'EQE_cone_fwd',NaN,'EQE_cone_opp',NaN,'asym',NaN);
end

save('BO_Freeform_result.mat', 'bestX', 'bestEQE', 'results', 'bd', ...
    'target_theta', 'target_phi', 'cone_half', 'W_CONTRAST', 'varNames', 'lb', 'ub', 'FF_XY', 'FF_INNER');
fprintf('\n######## Done ########\n');
fprintf('  목적함수 J (fwd - %.2f*opp) = %.5f ± %.5f\n', W_CONTRAST, bestEQE, candStd(ci));
fprintf('  EQE_total        = %.5f\n', bd.EQE_total);
fprintf('  EQE_cone_fwd (θ=%d,φ=%d)   = %.5f  <- HUD 방향 유용 효율\n', target_theta, target_phi, bd.EQE_cone_fwd);
fprintf('  EQE_cone_opp (θ=%d,φ=%d) = %.5f  <- 반대편 누설\n', target_theta, mod(target_phi+180,360), bd.EQE_cone_opp);
fprintf('  방위각 비대칭 asym = (fwd-opp)/(fwd+opp) = %.3f   [hemisphere→0]\n', bd.asym);
disp(bestX);
report_best_shape(table2array(bestX));


%% =====================================================================
%%  1회 평가 래퍼 (주기 재시작 + 크래시 -> NaN)
%% =====================================================================
function eqe = simulate_dirEQE(pt)
global ID_LT ltml ltloc eval_count restart_interval
eval_count = eval_count + 1;
if mod(eval_count, restart_interval) == 0
    fprintf('\n[Refresh] 시뮬 %d회. LightTools 재시작...\n', eval_count);
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_LT);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);  pause(2);
end
try
    eqe = objFcn_directionalEQE(pt).EQE_cone;   % 목적함수 J (음수/0 도 유효값)
catch err
    fprintf('\n[Error] eval %d LightTools 충돌: %s\n', eval_count, err.message);
    eqe = NaN;
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_LT);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
end
end

function obj = bo_objective(Xtbl)
obj = -simulate_dirEQE(table2array(Xtbl));
end

function f = polish_objective(x)
e = simulate_dirEQE(x);
if isnan(e), e = 0; end
f = -e;
end


%% =====================================================================
%%  Objective: freeform .ent 주입 + 나노 CPS + 2D 방향성 EQE
%% =====================================================================
function output = objFcn_directionalEQE(point)
global ID_LT ltml ltloc count ray_nums_current
global target_theta target_phi cone_half MESH_POS FF_N FF_INNER W_CONTRAST

lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

% --- 파라미터 언팩 ---
nInner   = numel(FF_INNER);
zInner   = point(1:nInner);
dETL     = point(nInner + 1);
dHTL     = point(nInner + 2);
stretchZ = point(nInner + 3);

% --- 고정 설정(기존과 동일 스케일) ---
d_sub = 1.3;  r_OLED = 1;  r_pat = 1;
Lensheight = 0.01;
wavelength_start = 580;  wavelength_end = 590;  n = 10;
if isempty(ray_nums_current), ray_nums = 50000; else, ray_nums = ray_nums_current; end

% --- (신규) freeform .ent 직접 생성 + 텍스처 unit-cell 에 물리기 ---
entPath = updateFreeformGeometry_ent(zInner, stretchZ);
if isempty(entPath)
    output = zero_output(); return;
end

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
sin089 = sind(0:89);  cos089 = cosd(0:89); %#ok<NASGU>
no_bar = no_bar(wavelength_start-399:wavelength_end-399,:);
ne_bar = ne_bar(wavelength_start-399:wavelength_end-399,:);
thickness = [100 dETL 25 10 dHTL 150];
EML_position = 4;  z0 = 12.5;  u_data_num = 499;  max_u = 3;

CPS_result = CPS_for_Isub(no_bar,ne_bar,thickness,emission_spectrum,eta_rad, ...
    horizontal_dipole_ratio,bottom_air_refractive_index,EML_position,z0,u_data_num,max_u,wavelength);
EQE_sub_CPS = CPS_result.EQE_sub;

TMF_p = TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_s = TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
R_p = abs(TMF_p.r_p).^2;  R_s = abs(TMF_s.r_s).^2;
Reflectance = (R_p + R_s)/2;

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

I_white = 0.5*(CPS_result.I_sub_s + CPS_result.I_sub_p);
P_white = I_white .* repmat(sin089, wavelength_num, 1);
weight_factor = sum(P_white, 2);

% --- 2D 원거리장 방향성 EQE (목표 원뿔 + 반대편 원뿔) ---
K = (wavelength_num-1)/n + 1;
Power_output = zeros(1, wavelength_num);
coneFrac_fwd = zeros(1, wavelength_num);   % 목표 방향 (theta_t, phi_t)
coneFrac_opp = zeros(1, wavelength_num);   % 반대 방위각 (theta_t, phi_t+180)
phi_opp = mod(target_phi + 180, 360);

for wv = 1:n:wavelength_num
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

    ltml.LTBegin(lt);  ltml.LTCmd(lt,'\V3D BeginAllSimulations');  ltml.LTEnd(lt);

    List = ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');  Key = ltml.LTListAtPos(lt,List,1);
    Power_output(wv) = ltml.LTDbGet(lt,Key,'TotalPower');

    Key = ltml.LTListAtPos(lt,List,MESH_POS);   % far-field 방향성 mesh
    [Igrid, thC, phC] = read_intensity_grid(ltml, lt, Key);
    coneFrac_fwd(wv) = cone_power_fraction(Igrid, thC, phC, target_theta, target_phi, cone_half);
    coneFrac_opp(wv) = cone_power_fraction(Igrid, thC, phC, target_theta, phi_opp,   cone_half);
end

weight_factor_2 = zeros(K,1);  Power_output_2 = zeros(K,1);
EQE_sub_matrix_2 = zeros(K,1);  cFwd_2 = zeros(K,1);  cOpp_2 = zeros(K,1);
for k = 1:K
    idx = n*(k-1) + 1;
    weight_factor_2(k)  = weight_factor(idx);
    Power_output_2(k)   = Power_output(idx);
    EQE_sub_matrix_2(k) = CPS_result.EQE_sub_matrix(idx);
    cFwd_2(k)           = coneFrac_fwd(idx);
    cOpp_2(k)           = coneFrac_opp(idx);
end
EQE_wv_matrix = Power_output_2 ./ weight_factor_2;
EQE_sub_matrix_2 = EQE_sub_matrix_2 / sum(EQE_sub_matrix_2) * EQE_sub_CPS;

wq = EQE_wv_matrix .* EQE_sub_matrix_2;   % 파장별 EQE 가중
EQE_total    = sum(wq);
EQE_cone_fwd = sum(wq .* cFwd_2);         % 목표 원뿔 EQE
EQE_cone_opp = sum(wq .* cOpp_2);         % 반대편 원뿔 EQE(누설)

% 목적함수: 목표 원뿔 EQE - W*반대편 누설 (hemisphere 는 fwd=opp -> 불리)
J = EQE_cone_fwd - W_CONTRAST * EQE_cone_opp;
asym = (EQE_cone_fwd - EQE_cone_opp) / max(EQE_cone_fwd + EQE_cone_opp, eps);

output = struct('EQE_total', EQE_total, 'EQE_cone', J, ...
    'EQE_cone_fwd', EQE_cone_fwd, 'EQE_cone_opp', EQE_cone_opp, 'asym', asym);

% 코팅 정리 (v4 계승)
List = ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');  Key = ltml.LTListByName(lt,List,'R_Al');
List = ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');  Key = ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName','R_temp');
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings].COATING[' sprintf('R_Bottom_%d', count) ']" Delete= \Q']);
fclose('all');
end

function output = zero_output()
% 지오메트리 생성/주입 실패: 목적함수를 NaN 으로(= bayesopt 오류점, GP 비오염).
output = struct('EQE_total',0,'EQE_cone',NaN,'EQE_cone_fwd',0,'EQE_cone_opp',0,'asym',0);
end


%% =====================================================================
%%  (신규) freeform .ent 직접 생성 + 텍스처 unit-cell Filename 지정
%% =====================================================================
function entPath = updateFreeformGeometry_ent(zInner, stretchZ)
% FrontSurface 내부점 Z 를 zInner 로, 경계는 0 으로 채운 25-벡터를 만들고,
% 템플릿 .ent 를 복제해 Z 만 치환한 .ent 를 쓴 뒤, 배열 모델 텍스처의
% LibraryElement Filename 을 그 .ent 로 지정한다. (import/SaveLibrary 불필요)
global ID_LT ltml ltloc FF_TEMPLATE FF_XY FF_INNER FF_N FF_BASE eval_count
entPath = '';
lt = ltloc.GetLTAPI(ID_LT);

% 전체 격자 Z: 경계 0, 내부 = zInner
zFull = zeros(FF_N, 1);
zFull(FF_INNER) = zInner(:);

% 고유 파일명(rng 오염 없이): 시간+eval 카운터
tag = sprintf('%s_%d', datestr(now,'yyyymmddHHMMSSFFF'), eval_count); %#ok<TNOW1,DATST>
entPath_out = [FF_BASE 'ff_' tag '.1.ent'];

try
    generate_freeform_ent(zFull, FF_TEMPLATE, entPath_out);
catch me
    fprintf('[Geom] .ent 생성 실패: %s\n', me.message);  return;
end
if ~exist(entPath_out, 'file')
    fprintf('[Geom] .ent 파일이 생성되지 않음.\n');  return;
end

% 텍스처 unit-cell 파일 + StretchZ 갱신 (v4 계승)
% @@VERIFY: 배열 모델 텍스처의 LibraryElement 가 이 freeform .ent 를 받도록 설정돼야 함.
try
    List = ltml.LTDbList(lt,'LENS_MANAGER[1]','LIBRARY_ELEMENT_UNIT_CELL');
    Key  = ltml.LTListByName(lt,List,'LibraryElement');
    ltml.LTDbSet(lt,Key,'Filename', entPath_out);
    List = ltml.LTDbList(lt,'LENS_MANAGER[1]','TEXTURE_PARAMETER');
    Key  = ltml.LTListByName(lt,List,'StretchZ');
    ltml.LTDbSet(lt,Key,'Value', stretchZ);
catch me
    fprintf('[Geom] 텍스처 갱신 실패: %s\n', me.message);  return;
end
entPath = entPath_out;
end


%% =====================================================================
%%  2D 인텐시티 그리드 읽기 + 원뿔 파워 분율
%% =====================================================================
function [Igrid, thetaC, phiC] = read_intensity_grid(ltml, lt, Key)
% far-field INTENSITY_MESH 를 [nLat x nLong] 로 읽는다. 차원/각도범위는 상단 전역.
% CellValue_UI(iLong, iLat): 첫 인자=longitude, 둘째=latitude (기존 코드와 동일).
global MESH_NLONG MESH_NLAT MESH_LONG_MIN MESH_LONG_MAX MESH_LAT_MIN MESH_LAT_MAX
nLong = MESH_NLONG;  nLat = MESH_NLAT;
Igrid = zeros(nLat, nLong);
for iL = 1:nLong
    for iT = 1:nLat
        v = ltml.LTDbGet(lt, Key, 'CellValue_UI', iL, iT);
        if isempty(v) || ~isfinite(v), v = 0; end   % 잘못된 인덱스 방어
        Igrid(iT, iL) = v;
    end
end
thetaC = MESH_LAT_MIN  + (MESH_LAT_MAX -MESH_LAT_MIN) /nLat  * ((1:nLat)  - 0.5);
phiC   = MESH_LONG_MIN + (MESH_LONG_MAX-MESH_LONG_MIN)/nLong * ((1:nLong) - 0.5);
end

function frac = cone_power_fraction(Igrid, thetaC, phiC, th_t, ph_t, halfAng)
% 목표 방향(th_t,ph_t) 중심 반각 halfAng 원뿔로 방출되는 파워 분율.
%   dP ∝ I(theta,phi) sin(theta).  원뿔 판정: 방향 단위벡터 각거리 <= halfAng.
[TH, PH] = ndgrid(thetaC, phiC);
W = Igrid .* sind(TH);
Wtot = sum(W(:));
if Wtot <= 0, frac = 0; return; end
d_t = ang2vec(th_t, ph_t);
V = ang2vec(TH, PH);
csep = reshape(V(1,:).*d_t(1) + V(2,:).*d_t(2) + V(3,:).*d_t(3), size(TH));
frac = sum(W(csep >= cosd(halfAng))) / Wtot;
end

function V = ang2vec(theta, phi)
th = theta(:).';  ph = phi(:).';
V = [sind(th).*cosd(ph); sind(th).*sind(ph); cosd(th)];
if isscalar(theta), V = V(:); end
end


%% =====================================================================
%%  형상 리포트 (격자 산점 -> 표면 보간; 순수 MATLAB)
%% =====================================================================
function report_best_shape(x)
global FF_XY FF_INNER FF_N
nInner = numel(FF_INNER);
zFull = zeros(FF_N,1);  zFull(FF_INNER) = x(1:nInner);
% 산점 보간으로 표면 시각화
F = scatteredInterpolant(FF_XY(:,1), FF_XY(:,2), zFull, 'natural', 'none');
g = linspace(min(FF_XY(:,1)), max(FF_XY(:,1)), 120);
[Xg,Yg] = meshgrid(g,g);
Zg = F(Xg,Yg);
m = Zg; m(~isfinite(m)) = 0;
cx = sum(Xg(:).*m(:))/sum(m(:));  cy = sum(Yg(:).*m(:))/sum(m(:));
figure('Name','Best asymmetric freeform (.ent grid)');
subplot(1,2,1); surf(Xg,Yg,Zg,'EdgeColor','none'); hold on;
plot3(FF_XY(:,1),FF_XY(:,2),zFull,'k.','MarkerSize',12);
axis tight; view(35,30); title(sprintf('Best freeform  z_{max}=%.3f', max(zFull))); xlabel x; ylabel y; zlabel z;
subplot(1,2,2); contourf(Xg,Yg,Zg,25,'LineColor','none'); axis equal tight; hold on;
plot(cx,cy,'r+','MarkerSize',14,'LineWidth',2); plot(0,0,'w.','MarkerSize',8);
title(sprintf('height centroid=(%+.3f,%+.3f)  [조향 \\rightarrow +x=phi 0]', cx, cy));
fprintf('[Shape] height centroid (조향 방향 지표) = (%+.4f, %+.4f)\n', cx, cy);
end

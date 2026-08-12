% ============================================================
%  Symmetric (rotationally-symmetric) swept-entity MLA + OLED co-optimization
%  -- GLOBAL optimizer (surrogateopt, feasibility-aware) + multi-start + polish
%
%  이 파일은 PSO_test_260401_..._auto_BO_rOLED_sweep_v4.m 를
%  "전역 최적화가 제대로 안 되는" 문제를 고쳐 다시 쓴 것이다.
%
%  [무엇이 문제였나]
%   - r_pat_list = 25  (단일 값)  ->  r_pat sweep / warm-start-from-prev-rpat /
%     단조성 자기진단(monotonicity) 로직이 전부 dead code (n_pat=1 이므로 s>1 미도달).
%     즉 실질은 "r_pat=25 고정, 13차원 EQE_total 단일 전역최적화" 문제인데,
%     전역 탐색이 bayesopt 단일 실행(13 DOF)이라 좁은 feasible 영역(~0.3%)에서
%     local optimum 에 잘 갇힌다.
%
%  [무엇을 고쳤나]  (물리/기하/스택/목적함수는 그대로 두고 "탐색 전략"만 교체)
%   (1) bayesopt(13D 단일) -> surrogateopt (RBF 대리모형, 표본효율 + 고차원에 강함)
%       + 여러 번 재시작(multi-start)으로 서로 다른 basin 을 탐색해 전역성 확보.
%       (프로젝트의 BO_asym_v3~v6 에서 이미 검증된 surrogateopt 로 통일)
%   (2) 좁은 feasible 영역 대응: 기하 제약 isValidPoints 를 surrogateopt 의
%       "비선형 제약(Ineq)" 으로 직접 연결 -> 대리모형이 feasible 쪽으로 샘플링을
%       스스로 편향. 게다가 infeasible 점은 값비싼 LightTools 시뮬 없이 즉시 반환
%       (기하검사만 수행)하므로 예산 낭비가 없다.
%   (3) 각 start 마다 valid 시드(rejection sampling) 다수 공급 + 직전 best 앵커
%       -> 대리모형이 feasible manifold 를 처음부터 넓게 커버.
%   (4) 각 start 종료 후 patternsearch 국소정련, 그리고 후보들(surrogate best,
%       polish, 전역 best 앵커)을 고정밀 ray 로 N_FINAL_REP 회 재평가해 승자 채택.
%       -> ray-tracing 노이즈로 인한 잘못된 최적점 선택 방지 (mean±std 기록).
%
%  [바뀌지 않은 것]  r_pat=25 고정, 대칭 7점 스플라인 swept-entity 형상,
%    ITO 스택(dAg 미사용), 목적함수 = EQE_total, CPS/TMM 물리, LightTools I/O.
%    -> freeform 렌즈(비대칭) 대비 "공정하게 공동최적화된 대칭/반구형 MLA" 기준선.
% ============================================================
clear;
%% For LightTools Connection
global ID_swept ID_LT ltml ltloc count r_pat eval_count restart_interval ray_nums_current
RenewLightTools();
try
    ltml.LTCmd(ltml.GetLTAPI(ID_LT), 'Message "Check Connection"');
catch
    ltml = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
end
count = 1;
restart_interval = 20;   % 시뮬 N회마다 LightTools 재시작
lt = ltloc.GetLTAPI(ID_swept);
ltx= getltpointer(ID_swept);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

%% ===== 고정 설정 =====
r_pat = 25;              % (고정) 대칭 MLA 패턴 크기. objFcn 이 global 로 사용.

%% ===== 정확도/시간 트레이드오프 =====
RAY_SEARCH   = 20000;    % 전역탐색/정련용 ray 수 (속도 우선, 노이즈 ~1/sqrt(N))
RAY_FINAL    = 100000;   % 최종 검증용 ray 수 (보고값 정밀도 우선)
N_FINAL_REP  = 3;        % 최종 후보 반복 평가 횟수 -> mean±std

%% ===== 전역 최적화(surrogateopt + multi-start) 예산 =====
NUM_STARTS       = 3;    % 서로 다른 시드로 surrogateopt 재시작 횟수 (전역성 ↑)
EVALS_PER_START  = 140;  % start 당 surrogateopt 최대 함수평가 (infeasible 은 저가)
MIN_SURR_POINTS  = 30;   % 초기 대리모형 표본 수 (>= 시드 수와 균형)
N_SEED_VALID     = 40;   % start 당 valid 초기 시드 개수 (feasible manifold 커버)
POLISH_EVALS     = 20;   % 각 start 후 patternsearch 국소정련 예산

%% Optimization Variables (13-dim: x2..x6, y2..y6, dETL, dHTL, stretchZ)
varNames = {'x2','x3','x4','x5','x6', 'y2','y3','y4','y5','y6', 'dETL','dHTL','stretchZ'};
lb = [0, 0, 0, 0, 0, 0,   0,   0,   0,   0,   10, 10, 0.1];
ub = [1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 1.5, 1.5, 150,150, 3];
nvar = numel(lb);

psOpts = optimoptions('patternsearch', ...
    'MaxFunctionEvaluations', POLISH_EVALS, ...
    'InitialMeshSize', 0.1, ...
    'MeshTolerance', 1e-3, ...
    'Cache', 'on', ...
    'Display', 'iter');

%% ===== Multi-start 전역 최적화 =====
gBestX    = [];          % 전역 best 설계변수 (row vector)
gBestEQE  = -inf;        % 전역 best 검증 EQE_total (고정밀 mean)
gBestStd  = NaN;
start_log = struct('surrEQE',{},'polishEQE',{},'bestEQE',{},'bestStd',{},'evals',{});

for st = 1:NUM_STARTS
    fprintf('\n############ Global start %d/%d (surrogateopt) ############\n', st, NUM_STARTS);

    % start 시작 시 LightTools 클린 상태 보장
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);
    ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    eval_count = 0;

    % --- valid 초기 시드(+직전 전역 best 앵커) ---
    seedMat = genValidPoints(N_SEED_VALID, lb, ub);
    if ~isempty(gBestX)
        seedMat = [gBestX; seedMat];   % 앵커: 이전 basin 을 잃지 않도록
    end
    initPts = struct('X', seedMat);

    % --- (1) 전역 탐색: surrogateopt (저정밀 ray) ---
    ray_nums_current = RAY_SEARCH;
    sopts = optimoptions('surrogateopt', ...
        'MaxFunctionEvaluations', EVALS_PER_START, ...
        'MinSurrogatePoints',     MIN_SURR_POINTS, ...
        'InitialPoints',          initPts, ...
        'UseParallel',            false, ...
        'PlotFcn', [], 'Display', 'iter');
    % surrogate_objconstr: 구조체(Fval, Ineq) 반환 -> surrogateopt 가
    %   feasibility-aware 로 동작. infeasible 은 시뮬 없이 즉시 반환.
    [xS, fS, ~, outS] = surrogateopt(@surrogate_objconstr, lb, ub, sopts);

    surrEQE = NaN;
    if ~isempty(xS) && isValidPoints(xS(:).')
        surrEQE = -fS;
        xS = xS(:).';
    else
        fprintf('[Warn] start %d: surrogateopt 이 feasible 해를 반환하지 못함.\n', st);
        xS = [];
    end

    % --- (2) 국소 정련: patternsearch (surrogate best 에서 출발) ---
    xPol = []; polishEQE = NaN;
    x0 = xS;
    if isempty(x0) && ~isempty(gBestX), x0 = gBestX; end
    if ~isempty(x0)
        fprintf('--- Local polish (patternsearch, %d evals) ---\n', POLISH_EVALS);
        ray_nums_current = RAY_SEARCH;
        try
            xPol = patternsearch(@polish_objective, x0, [],[],[],[], lb, ub, [], psOpts);
            xPol = xPol(:).';
            polishEQE = -polish_objective(xPol);
        catch perr
            fprintf('[Warn] patternsearch 실패(%s). surrogate 결과만 사용.\n', perr.message);
            xPol = [];
        end
    end

    % --- (3) 최종 검증: 후보들을 고정밀 ray 로 반복 평가 후 승자 채택 ---
    candX = {};
    if ~isempty(xS),    candX{end+1} = xS;    end
    if ~isempty(xPol) && (isempty(xS) || ~isequal(xPol, xS)), candX{end+1} = xPol; end
    if ~isempty(gBestX) && ~any(cellfun(@(c) isequal(c, gBestX), candX))
        candX{end+1} = gBestX;   % 전역 best 앵커도 현 조건에서 재검증 -> 하한선
    end

    ray_nums_current = RAY_FINAL;
    candMean = -inf(1, numel(candX));
    candStd  = zeros(1, numel(candX));
    for c = 1:numel(candX)
        if ~isValidPoints(candX{c}), continue; end
        e = nan(1, N_FINAL_REP);
        for rrep = 1:N_FINAL_REP
            e(rrep) = simulate_EQE(candX{c});
        end
        candMean(c) = mean(e, 'omitnan');
        candStd(c)  = std(e, 'omitnan');
        fprintf('  start %d cand %d: EQE_total = %.5f ± %.5f (N=%d, %d rays)\n', ...
            st, c, candMean(c), candStd(c), N_FINAL_REP, RAY_FINAL);
    end
    [bestEQE, ci] = max(candMean);

    if isfinite(bestEQE) && bestEQE > gBestEQE
        gBestEQE = bestEQE;
        gBestStd = candStd(ci);
        gBestX   = candX{ci};
        fprintf('  [Global] start %d 에서 전역 best 갱신: EQE_total = %.5f ± %.5f\n', ...
            st, gBestEQE, gBestStd);
    else
        fprintf('  [Global] start %d: 전역 best 갱신 없음 (현 best = %.5f).\n', st, gBestEQE);
    end

    start_log(st).surrEQE   = surrEQE;
    start_log(st).polishEQE = polishEQE;
    start_log(st).bestEQE   = bestEQE;
    start_log(st).bestStd   = candStd(min(ci,numel(candStd)));
    start_log(st).evals     = outS.funccount;

    % 체크포인트 저장 (중간에 죽어도 최선값 보존)
    save('symmetric_MLA_global_result.mat', 'gBestX', 'gBestEQE', 'gBestStd', ...
        'varNames', 'r_pat', 'start_log', 'lb', 'ub');
    fprintf('############ start %d done | 전역 best EQE_total = %.5f ± %.5f ############\n', ...
        st, gBestEQE, gBestStd);
end

%% ===== 결과 요약 =====
disp('=== Symmetric MLA global optimization finished ===');
if isempty(gBestX)
    warning('feasible 해를 찾지 못했습니다. 시드/제약을 점검하세요.');
else
    bestT = array2table(gBestX, 'VariableNames', varNames);
    fprintf('\n######## Best symmetric MLA (r_pat=%g) ########\n', r_pat);
    fprintf('  EQE_total = %.5f ± %.5f  (고정밀 %d rays, N=%d)\n', ...
        gBestEQE, gBestStd, RAY_FINAL, N_FINAL_REP);
    disp('  best design variables:'); disp(bestT);
    save('symmetric_MLA_global_result.mat', 'gBestX', 'gBestEQE', 'gBestStd', ...
        'varNames', 'r_pat', 'start_log', 'lb', 'ub', 'bestT');
end

% start 별 수렴 추이 (전역성 점검용)
figure('Name','symmetric MLA multi-start','Color','w');
sEQE = arrayfun(@(s) s.bestEQE, start_log);
bar(1:NUM_STARTS, sEQE); grid on;
xlabel('global start #'); ylabel('verified best EQE\_total');
title(sprintf('Symmetric MLA co-optimization (r\\_pat=%g): multi-start best', r_pat));


%% ===== LightTools 1회 평가 공용 래퍼 =====
% 주기적 재시작 + 크래시 처리. 크래시/기하오류는 NaN 반환.
function eqe = simulate_EQE(pt)
global ID_swept ltml ltloc eval_count restart_interval

eval_count = eval_count + 1;
if mod(eval_count, restart_interval) == 0
    fprintf('\n[Refresh] 시뮬레이션 %d회 수행. LightTools를 재시작합니다...\n', eval_count);
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);
    ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    pause(2);
end

try
    eqe = objFcn_angularEQE(pt).EQE_total;
    if eqe == 0
        eqe = NaN;   % 파셋 불일치 등 기하 오류: 값 0이 아니라 "평가 실패"로 처리
    end
catch err
    fprintf('\n[Error] eval %d 평가 중 LightTools 충돌: %s\n', eval_count, err.message);
    eqe = NaN;
    fprintf('LightTools를 긴급 재시작합니다...\n');
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);
    ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
end
end

%% ===== surrogateopt 목적+제약 결합 함수 =====
% 반환 구조체:
%   .Ineq <= 0  이면 feasible (기하 valid), > 0 이면 infeasible
%   .Fval       최소화 대상 (= -EQE_total). infeasible 이면 시뮬 없이 상수.
% -> surrogateopt 가 제약 대리모형으로 feasible 영역 위주 샘플링.
%    infeasible 점은 값비싼 LightTools 호출 없이 즉시 반환(예산 절약).
function out = surrogate_objconstr(x)
x = x(:).';
if ~isValidPoints(x)
    out.Ineq = 1;      % infeasible (기하 위반)
    out.Fval = 1;      % 임의 상수(신뢰 안 함) - 시뮬 미수행
    return;
end
e = simulate_EQE(x);
if ~isfinite(e)
    out.Ineq = 1;      % 시뮬 실패 지점도 회피하도록 infeasible 처리
    out.Fval = 1;
else
    out.Ineq = -1;     % feasible
    out.Fval = -e;     % EQE 최대화 == -EQE 최소화
end
end

%% ===== patternsearch 정련용 목적함수 =====
function f = polish_objective(x)
x = x(:).';
if ~isValidPoints(x)
    f = 0;          % 무효 형상: 시뮬 없이 벌점 (valid면 -EQE < 0 이므로 항상 열등)
    return;
end
e = simulate_EQE(x);
if isnan(e), e = 0; end
f = -e;
end

%% ===== 무작위 valid 시드 생성 (rejection sampling) =====
function P = genValidPoints(K, lb, ub)
dim = numel(lb);
P = zeros(K, dim);
for i = 1:K
    ok = false;
    while ~ok
        p = lb + rand(1, dim) .* (ub - lb);
        if isValidPoints(p)
            ok = true;
            P(i, :) = p;
        end
    end
end
end


%% Objective Function (ray 수 가변 + .coa fclose + rng 비오염 파일명)
function output = objFcn_angularEQE(point)
global ID_LT ID_swept ltml ltloc count r_pat ray_nums_current
% Define segment length and other necessary parameters
lt = ltloc.GetLTAPI(ID_LT);  % lenssizeeffect
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

d_sub = 1.3;
r_OLED=1;
x_pattern=r_pat;
y_pattern=r_pat;
Lensheight=0.01;
wavelength_start=450;
wavelength_end=750;
n=10; % step size for wavelength

% 탐색/검증 단계에 따라 ray 수 가변
if isempty(ray_nums_current)
    ray_nums = 50000;
else
    ray_nums = ray_nums_current;
end

List=ltml.LTDbList(lt,'lens_manager[1]','SIMULATIONS');
Key=ltml.LTListByName(lt,List,'ForwardAll');
ltml.LTDbSet(lt,Key,'MaxProgress',ray_nums);
List=ltml.LTDbList(lt,'lens_manager[1]','CUBE_PRIMITIVE');
Key=ltml.LTListByName(lt,List,'Substrate');
ltml.LTDbSet(lt,Key,'Height',d_sub);
ltml.LTDbSet(lt,Key,'Y',d_sub/2);
SRList=ltml.LTDbList(lt,'lens_manager[1]','CUBE_PRIMITIVE');
SRKey=ltml.LTListAtPos(lt,SRList,2);
ltml.LTDbSet(lt,SRKey,'Y',d_sub+Lensheight/2);
List=ltml.LTDbList(lt,'lens_manager[1]','TEXTURE_ZONE_EXTENT');
Key=ltml.LTListByName(lt,List,'zone');
ltml.LTDbSet(lt,Key,'Geometry_1',x_pattern);
ltml.LTDbSet(lt,Key,'Geometry_2',y_pattern);
List=ltml.LTDbList(lt,'lens_manager[1]','DISK_SOURCE');
Key=ltml.LTListByName(lt,List,'DiskSource_18');
ltml.LTDbSet(lt,Key,'Radius',r_OLED);

% passing input points
x2 = point(1);  x3 = point(2);  x4 = point(3);  x5 = point(4);  x6 = point(5);
y2 = point(6);  y3 = point(7);  y4 = point(8);  y5 = point(9);  y6 = point(10);
dETL = point(11); dHTL = point(12);
% dAg = point(13);
stretchZ=point(13);

% Create spline control points
xy = zeros(7,2);
xy(1,:) = [0, 1];
xy(7,:) = [1, 0];
xy(2,:) = [x2, y2];
xy(3,:) = [x3, y3];
xy(4,:) = [x4, y4];
xy(5,:) = [x5, y5];
xy(6,:) = [x6, y6];

lt = ltloc.GetLTAPI(ID_swept); % swept entity
ltx= getltpointer(ID_swept);  % swept entity
lt2 = ltloc.GetLTAPI(ID_LT); % LT simulation

Curve="LENS_MANAGER[1].COMPONENTS[Components].SWEPT_SOLID[SweptEntity].SWEPT_PRIMITIVE[SweptPrimitive].SWEPT_PROFILE[SweptProfile].FITTED_CURVE[SweptSurface_1]";
ltx.SetSweptProfilePoints(Curve,xy,7); % 7*2 double
ltx.DbSet(Curve,'StartSlopeMode',"Auto");
ltx.DbSet(Curve,'EndSlopeMode',"Auto");

List=ltml.LTDbList(lt,'LENS_MANAGER[1]','FITTED_CURVE');
Key=ltml.LTListByName(lt,List,'SweptSurface_1');

ltml.LTDbSet(lt, Key,'NumFacets',100);
x_values = zeros(101,1);

for a=1:101
    x_values(a)=ltml.LTDbGet(lt,Key,'YFacetsAt',a);
end
max_length = max(x_values);

if max_length > 1
    xy = xy / max_length;
end

ltx.SetSweptProfilePoints(Curve,xy,7); % 7*2 double
ltx.DbSet(Curve,'StartSlopeMode',"Auto");
ltx.DbSet(Curve,'EndSlopeMode',"Auto");

xy_l = zeros(7,2); % x,y coordinates in LightTools

for j=1:7
    xy_l(j,1) = ltml.LTDbGet(lt, Key, 'YAt',j);
    xy_l(j,2) = ltml.LTDbGet(lt, Key, 'ZAt',j);
end

tol = 1e-4;  % 필요시 조정
if max(abs(xy(:) - xy_l(:))) > tol
    output = struct();
    output.EQE_0_20 = 0;
    output.EQE_20_40 = 0;
    output.EQE_40_60 = 0;
    output.EQE_60_80 = 0;
    output.EQE_total = 0;
    return;
end


% File name and path configuration
rng('shuffle')
strLength = 10;
charSet = ['a':'z' 'A':'Z' '0':'9'];
numChars = length(charSet);
randIndices = randi(numChars, 1, strLength);
index = charSet(randIndices);

pathname = '"C:\Users\jhkim\Desktop\Green_CE_Calculation\swept_';
pathname_unrepaired = '"C:\Users\jhkim\Desktop\Green_CE_Calculation\unrepaired\swept_unrepaired_';

totalpath = [pathname index '.ent"'];
totalpath_unrepaired = [pathname_unrepaired index '.ent"'];

ltml.LTCmd(lt, 'DefaultSelect "SweptEntity.tag_1"');
ltml.LTCmd(lt, sprintf('SaveLibrary XYZ 0,0,0 %s ', totalpath_unrepaired));
ltml.LTCmd(lt, 'DefaultSelect "SweptEntity.tag_1"');
ltml.LTCmd(lt, 'RepairEntities');
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
ltml.LTSetOption(lt2, "ShowFileDialogBox", 0);
ltml.LTCmd(lt, 'DefaultSelect "SweptEntity.tag_1"');
ltml.LTCmd(lt, sprintf('SaveLibrary XYZ 0,0,0 %s ', totalpath));
ltml.LTCmd(lt, 'Undo');
ltml.LTCmd(lt, 'Undo');

totalpathmod = [pathname index '.1.ent"'];

List = ltml.LTDbList(lt2, 'LENS_MANAGER[1]', 'LIBRARY_ELEMENT_UNIT_CELL');
Key = ltml.LTListByName(lt2, List, 'LibraryElement');
ltml.LTDbSet(lt2, Key, 'Filename', totalpathmod);

List = ltml.LTDbList(lt2, 'LENS_MANAGER[1]', 'TEXTURE_PARAMETER');
Key = ltml.LTListByName(lt2, List, 'StretchZ');
ltml.LTDbSet(lt2, Key, 'Value', stretchZ);

%% Define layer (CPS)
load('nk_JH33.mat');
load('Photopic_400_800.mat');
load('CIE_1931.mat');
load('R_pd.mat');
wavelength=(wavelength_start:wavelength_end).';

wavelength_num=length(wavelength);
emission_spectrum=spectrum.l_I_Irdmppyph2tmd(wavelength_start-399:wavelength_end-399,:);
eta_rad=0.98;
horizontal_dipole_ratio=0.865;
bottom_air_refractive_index=ones(wavelength_num,1);

no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
layer_num=size(no_bar,2);
sin089=sind(0:89);
cos089=cosd(0:89);
no_bar=no_bar(wavelength_start-399:wavelength_end-399,:);
ne_bar=ne_bar(wavelength_start-399:wavelength_end-399,:);
thickness=[100 dETL 25 10 dHTL 150];

EML_position=4; % count from left side (+air)
z0=12.5;
u_data_num=499;
max_u=3;

CPS_result=CPS_for_Isub(no_bar,ne_bar,thickness,emission_spectrum,eta_rad,horizontal_dipole_ratio,bottom_air_refractive_index,EML_position,z0,u_data_num,max_u,wavelength);
EQE_air_CPS=CPS_result.EQE_air;
EQE_sub_CPS=CPS_result.EQE_sub;

%% bottom reflectance
TMF_OLED_bottom_p=TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_OLED_bottom_s=TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);

R_p_bottom=abs(TMF_OLED_bottom_p.r_p).^2;
T_p_bottom=no_bar(:,1)./no_bar(:,layer_num)*(1./cos089).*sqrt(1-(ne_bar(:,layer_num)./ne_bar(:,1)*sin089).^2).*abs(TMF_OLED_bottom_p.t_p).^2;

R_s_bottom=abs(TMF_OLED_bottom_s.r_s).^2;
T_s_bottom=no_bar(:,1)./no_bar(:,layer_num)*(1./cos089).*sqrt(1-(no_bar(:,layer_num)./no_bar(:,1)*sin089).^2).*abs(TMF_OLED_bottom_s.t_s).^2;

for i=1:wavelength_num
    T_p_bottom(i,ceil(asind(ne_bar(i,1)/ne_bar(i,layer_num)))+1:end)=0;
    T_s_bottom(i,ceil(asind(no_bar(i,1)/no_bar(i,layer_num)))+1:end)=0;
end

Transmittance=(T_p_bottom+T_s_bottom)/2;
Reflectance=(R_p_bottom+R_s_bottom)/2;

%% Coating (.mat to .coa)
lt = ltloc.GetLTAPI(ID_LT); % LT simulation
fileID = fopen(sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count), 'w');
fprintf(fileID,'%s\n%s%d\n%s\n%s\n%s\n%s\n ','DFAT Version 1.0', 'DATANAME: R_Bottom_',count, 'ABSORBING: YES', 'INDEX: 1.51', 'DATAITEMS: TAVG RAVG');
for i=wavelength_start:wavelength_end
    fprintf(fileID,'%s  %d\n','wv',i);
    for j=0:89
        fprintf(fileID,'%s  %d  %d  %.3f\n', 'AOI',j, 0, Reflectance(i-wavelength_start+1,j+1));
    end
end
fclose(fileID);  % LightTools가 읽기 전에 버퍼 플러시 + 파일 잠금 해제

ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings]" LoadFileName="' sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count) '"']);

List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');
Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');
Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName',sprintf('R_Bottom_%d', count));

%%
I_white=0.5*(CPS_result.I_sub_s+CPS_result.I_sub_p); % s랑 p 따로 구분하지 않음 일단
sin089=sind(0:89);
P_white=I_white.*repmat(sin089,wavelength_num,1);
weight_factor=sum(P_white,2); % I_white : I_sub의 파장별 intensity 301x90행렬
I_white_ang=sum(P_white);
%     weight_factor(1,1)=weight_factor(2,1);

wavelength_num=length(wavelength);

I_air_1_2=zeros(90,(wavelength_num+n-1)/n);
Luminance=cell((wavelength_num+n-1)/n,1);
Ray_wv=zeros(1,(wavelength_num+n-1)/n);
Cell_flux= zeros((wavelength_num+n-1)/n,9);
for wv=1:n:wavelength_num
    fileID = fopen('C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','w');
    fprintf(fileID,'%s  %d  %d  %d  %d  %d  %d','SPHEREMESH:',1, 90, 0, 0, 360, 90);
    writematrix(flip(I_white(wv,:).'),'C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','Delimiter','tab','WriteMode','append');
    fclose(fileID);
    SRList=ltml.LTDbList(lt, 'Lens_manager[1]','DISK_SOURCE');
    SRKey=ltml.LTListAtPos(lt,SRList,1);
    ltml.LTDbSet(lt,SRKey,'Radiant_Power', weight_factor(wv)); % 파장에 따른 파워를 다르게 설정, 그 안에서 각도별 파워는 grid에서 조정
    for k=1:1  % 예전에 광원 많았을때는 k=1:광원수 였었음
        SRList=ltml.LTDbList(lt, 'Lens_manager[1]','Spectral_region');
        SRKey=ltml.LTListAtPos(lt,SRList,k+1);
        ltml.LTDbSet(lt,SRKey,'Spectral_Definition', 'Monochromatic');
        ltml.LTDbSet(lt,SRKey,'Single_Wavelength', wv+wavelength_start-1);
        List=ltml.LTDbList(lt,'lens_manager[1]','DIRECTION_GRID_APODIZER');
        Key=ltml.LTListAtPos(lt,List,k);
        pathname='C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\'; % have to change pathname
        ltml.LTDbSet(lt,Key,'LoadFileName',[pathname sprintf('AI_temp.txt')]);
    end
    %% 시뮬레이션 및 후처리
    ltml.LTBegin(lt);
    ltml.LTCmd(lt,'\V3D BeginAllSimulations');
    ltml.LTEnd(lt);
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,1);
    Power_output(wv)=ltml.LTDbGet(lt,Key,'TotalPower');  % [W]
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,2);
    Power_output_30(wv)=ltml.LTDbGet(lt,Key,'TotalPower');  % [W]
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,3);
    for j=1:90
        I_air_1_JH(91-j,:)=ltml.LTDbGet(lt,Key,'CellValue_UI',1,91-j);
    end
    I_air_1_2(:,(wv+n-1)/n)=smooth(I_air_1_JH);
    %     I_air_1_2(:,(wv+n-1)/n)=I_air_1_JH;
end

K = (wavelength_num-1)/n + 1;

weight_factor_2  = zeros(K,1);
Power_output_2   = zeros(K,1);
EQE_sub_matrix_2 = zeros(K,1);

for k = 1:K
    idx = n*(k-1) + 1;

    weight_factor_2(k)  = weight_factor(idx);
    Power_output_2(k)   = Power_output(idx);
    EQE_sub_matrix_2(k) = CPS_result.EQE_sub_matrix(idx);
end

EQE_wv_matrix = Power_output_2 ./ weight_factor_2;  % (Kx1)

% 3) Normalize CPS spectral EQE_sub distribution to match EQE_sub_CPS
EQE_sub_matrix_2 = EQE_sub_matrix_2 / sum(EQE_sub_matrix_2) * EQE_sub_CPS;  % (Kx1)

% 4) Total EQE after optics
EQE_total = sum(EQE_wv_matrix .* EQE_sub_matrix_2);

% 5) Angular EQEs using LT angular intensity distribution per sampled wavelength
EQE_0_20   = 0;
EQE_20_40  = 0;
EQE_40_60  = 0;
EQE_60_80  = 0;

sin_col = sin089(:);  % 90x1 for elementwise multiply

for k = 1:K
    % Per-wavelength contribution to total EQE
    contrib_k = EQE_wv_matrix(k) * EQE_sub_matrix_2(k);

    % Angular radiant intensity vs theta for this wavelength sample
    I_theta = I_air_1_2(:,k);  % 90x1, theta = 0..89 deg

    % Convert to proportional angular power weights (constants cancel in fractions)
    W_theta = I_theta .* sin_col;  % 90x1, proportional to dP/dtheta integrated over azimuth
    W_tot   = sum(W_theta);

    % Fractions in bins (using [a,b) convention)
    f_0_20   = sum(W_theta(1:20))   / W_tot;  % 0..19 deg
    f_20_40  = sum(W_theta(21:40))  / W_tot;  % 20..39 deg
    f_40_60  = sum(W_theta(41:60))  / W_tot;  % 40..59 deg
    f_60_80  = sum(W_theta(61:80))  / W_tot;  % 60..79 deg

    % Accumulate angular EQEs
    EQE_0_20   = EQE_0_20   + contrib_k * f_0_20;
    EQE_20_40  = EQE_20_40  + contrib_k * f_20_40;
    EQE_40_60  = EQE_40_60  + contrib_k * f_40_60;
    EQE_60_80  = EQE_60_80  + contrib_k * f_60_80;
end

output = struct();
output.EQE_0_20 = EQE_0_20;
output.EQE_20_40 = EQE_20_40;
output.EQE_40_60 = EQE_40_60;
output.EQE_60_80 = EQE_60_80;
output.EQE_total = EQE_total;

List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');
Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');
Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName','R_temp');
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings].COATING[' sprintf('R_Bottom_%d', count) ']" Delete= \Q']);
fclose('all');

end

%% Spline Constraints Function (기존 코드 그대로 유지)
function TF = isValidPoints(X)
% X: N x 12+ matrix (numeric). 앞 10개 열(x2~x6, y2~y6)만 사용.
numRows = size(X,1);
numPts  = 7;
TF = true(numRows,1);

for k = 1:numRows
    x = [0, X(k,1:5), 1];    % x2~x6
    y = [1, X(k,6:10), 0];   % y2~y6

    violates = false;

    % (1) Intersection
    for i = 1:numPts - 1
        for j = i + 2:numPts - 1
            if i == 1 && j == numPts - 1
                continue;
            end
            if checkIntersection([x(i), y(i)], [x(i+1), y(i+1)], ...
                    [x(j), y(j)], [x(j+1), y(j+1)])
                violates = true;
                break;
            end
        end
        if violates, break; end
    end

    % (2) Collinearity
    if ~violates
        for i = 1:numPts - 2
            if isCollinear([x(i), y(i)], [x(i+1), y(i+1)], [x(i+2), y(i+2)])
                violates = true;
                break;
            end
        end
    end

    % (3) Spacing
    if ~violates
        minD = 0.05; maxD = 1.0;
        d = hypot(diff(x), diff(y));
        if any(d < minD | d > maxD)
            violates = true;
        end
    end

    % (4) Angle
    if ~violates
        maxAng = 2 * pi / 3;
        for i = 2:numPts - 1
            v1 = [x(i), y(i)] - [x(i-1), y(i-1)];
            v2 = [x(i+1), y(i+1)] - [x(i), y(i)];
            ang = atan2(norm(cross([v1,0], [v2,0])), dot(v1, v2));
            if ang > maxAng
                violates = true;
                break;
            end
        end
    end

    TF(k) = ~violates;
end

% === Helper Functions ===
    function isCol = isCollinear(p1, p2, p3)
        area = 0.5 * det([p1 1; p2 1; p3 1]);
        isCol = abs(area) < 1e-5;
    end

    function intersects = checkIntersection(p1, p2, p3, p4)
        function o = orientation(p, q, r)
            o = (q(2) - p(2)) * (r(1) - q(1)) - (q(1) - p(1)) * (r(2) - q(2));
        end
        o1 = orientation(p1, p2, p3);
        o2 = orientation(p1, p2, p4);
        o3 = orientation(p3, p4, p1);
        o4 = orientation(p3, p4, p2);
        intersects = (o1 * o2 < 0) && (o3 * o4 < 0);
    end
end


function RenewLightTools()
global ID_LT ID_swept ltml ltloc lt
lt_exe_path = 'C:\Program Files\Optical Research Associates\LightTools 2023.03\lt.exe';
model_file_path_swept = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\SweptEntity.2.lts';
model_file_path_LT = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\Lens_size_effect_for_PSO_bump_modified_v1.1.lts';
% =========================================================================

fprintf('--- Restarting LightTools ---\n');

% 1. 기존 LightTools 강제 종료
target_user = 'jhkim';
kill_cmd = sprintf('taskkill /F /FI "USERNAME eq %s" /IM lt.exe', target_user);
[~, ~] = system(kill_cmd);
pause(2);

% 2. 시스템 명령어로 .lts 파일 직접 실행
cmd = sprintf('"%s" "%s" &', lt_exe_path, model_file_path_swept);
status = system(cmd);
% 2. LightTools 재실행 및 연결
try
    ltml = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
catch
    error('LightTools 재시작 실패. 라이선스나 설치 상태를 확인하세요.');
end

find_cmd = sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);

[status, cmdout] = system(find_cmd);
if status == 0 && contains(cmdout, 'lt.exe')
    tokens = regexp(cmdout, '"(\d+)"', 'tokens');
    if ~isempty(tokens)
        pid_str = tokens{1}{1};
        ID_swept = str2double(pid_str);
        fprintf('PID found for user %s: %d\n', target_user, ID_swept);
    else
        error('프로세스는 찾았으나 PID 추출 실패. 정규식 확인 필요.');
    end
else
    error('사용자 %s 로 실행된 LightTools(lt.exe)를 찾을 수 없습니다.', target_user);
end
cmd = sprintf('"%s" "%s" &', lt_exe_path, model_file_path_LT);

status = system(cmd);
% 2. LightTools 재실행 및 연결
find_cmd = sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);

[status, cmdout] = system(find_cmd);
if status == 0 && contains(cmdout, 'lt.exe')
    tokens = regexp(cmdout, '"(\d+)"', 'tokens');
    if ~isempty(tokens)
        pid_str = tokens{3}{1};
        ID_LT = str2double(pid_str);
        fprintf('PID found for user %s: %d\n', target_user, ID_LT);
    else
        error('프로세스는 찾았으나 PID 추출 실패. 정규식 확인 필요.');
    end
else
    error('사용자 %s 로 실행된 LightTools(lt.exe)를 찾을 수 없습니다.', target_user);
end
pause(5);
end

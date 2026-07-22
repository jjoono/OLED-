% ============================================================
%  Freeform swept-lens + OLED co-optimization  --  FAST global optimizer
%  목적함수: EQE_total  (반구형 MLA 55% 에 대항하는 freeform 60% 타겟)
%
%  [무엇을 대체하나]  업로드된 PSO 코드(PSO_test_260205_stretch_z_auto.m)
%   - PSO는 초기화 200 + 20 iter × 200 = 최대 4,200회 평가.
%     full-wavelength LightTools(453–753, n=10 → evaluation당 각도시뮬 31회)면
%     현실적으로 며칠~주 단위. 표본효율이 매우 나쁨(비싼 블랙박스에 부적합).
%
%  [내 선택]  surrogateopt (RBF 대리모형) + multi-fidelity + multi-start + polish
%   - 이유(성능): 13 DOF·비싼 블랙박스·기하제약 → GP(bayesopt)보다 RBF surrogate가
%     적합하고 표본효율이 PSO보다 1~2 order 좋음.
%   - 이유(공정성): 대칭 baseline(BO_symmetric_MLA_global_v5.m)과 최적화기/예산/
%     변수 박스를 동일하게 맞춰 "freeform한테만 좋은 옵티마이저 줬다"는 반박을 차단.
%   - [핵심 속도] multi-fidelity: 탐색은 저파장해상도(n=30→K=11)+저ray, 최종 검증만
%     full(n=10→K=31)+고ray. 렌즈 전달효율 EQE_wv 는 파장에 매끄러워 형상 순위가
%     보존됨(EQE_sub 는 CPS full 해상도라 renormalize로 편향도 작음).
%
%  [원본 대비 버그 수정]
%   (1) line 302 `~isequal(xy,xy_l)` : float 완전일치라 사실상 항상 0 반환.
%       → tolerance 비교(max(abs(Δ))>tol)로 교체.
%   (2) .coa 파일을 LightTools가 읽기 전에 fclose 누락 → 미플러시 값 흔들림.
%       → coa 작성 직후 fclose 추가.
%
%  [바뀌지 않은 것]  swept 7점 스플라인 형상, ITO 스택, x_pattern=15, d_sub=1.295,
%    CPS/TMM 물리, 파장창 453–753. (freeform 형상 자유도는 대칭 baseline과 동일 박스)
% ============================================================
clear;
%% For LightTools Connection
global ID_swept ID_LT ltml ltloc count eval_count restart_interval ...
       ray_nums_current wave_n_current
RenewLightTools();
try
    ltml.LTCmd(ltml.GetLTAPI(ID_LT), 'Message "Check Connection"');
catch
    ltml = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
end
count = 1;
restart_interval = 20;
lt = ltloc.GetLTAPI(ID_swept);
ltx= getltpointer(ID_swept);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

%% ===== Multi-fidelity 설정 (속도의 핵심) =====
% 탐색(저정밀): 파장 성기게 + ray 적게   |   검증(고정밀): full 파장 + ray 많이
WAVE_N_SEARCH = 30;      % 파장 step (453:30:753 → K=11) : 탐색용
WAVE_N_FINAL  = 10;      % 파장 step (453:10:753 → K=31) : full 검증용
RAY_SEARCH    = 10000;   % 탐색 ray 수
RAY_FINAL     = 50000;   % 검증 ray 수
N_FINAL_REP   = 3;       % 최종 후보 반복 평가 → mean±std
% [주의] (wavelength_num-1)=300 이 WAVE_N 으로 나누어떨어져야 함(30,10 모두 OK).

%% ===== 전역 최적화(surrogateopt + multi-start) 예산 =====
NUM_STARTS      = 3;     % 서로 다른 시드로 재시작 (전역성)
EVALS_PER_START = 140;   % start 당 surrogateopt 최대 함수평가 (infeasible 은 저가)
MIN_SURR_POINTS = 30;
N_SEED_VALID    = 40;    % start 당 valid 초기 시드
POLISH_EVALS    = 20;    % patternsearch 국소정련

%% Optimization Variables (13-dim: x2..x6, y2..y6, dETL, dHTL, stretchZ)
% [정리] 원본 lb/ub 는 14-dim 이었으나 objFcn 이 실제 쓰는 건 13개뿐이고
%   stretchZ 범위도 [5,30]으로 잘못 들어가 있었음 → 대칭 baseline 과 동일한
%   13-dim 박스(stretchZ ∈ [0.1,3])로 정정.
varNames = {'x2','x3','x4','x5','x6', 'y2','y3','y4','y5','y6', 'dETL','dHTL','stretchZ'};
lb = [0, 0, 0, 0, 0, 0,   0,   0,   0,   0,   10, 10, 0.1];
ub = [1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 1.5, 1.5, 150,150, 3];
nvar = numel(lb);

psOpts = optimoptions('patternsearch', ...
    'MaxFunctionEvaluations', POLISH_EVALS, ...
    'InitialMeshSize', 0.1, 'MeshTolerance', 1e-3, ...
    'Cache', 'on', 'Display', 'iter');

%% ===== Multi-start 전역 최적화 =====
gBestX   = [];  gBestEQE = -inf;  gBestStd = NaN;
start_log = struct('surrEQE',{},'polishEQE',{},'bestEQE',{},'bestStd',{},'evals',{});

for st = 1:NUM_STARTS
    fprintf('\n############ Global start %d/%d (surrogateopt) ############\n', st, NUM_STARTS);
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);
    ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    eval_count = 0;

    seedMat = genValidPoints(N_SEED_VALID, lb, ub);
    if ~isempty(gBestX), seedMat = [gBestX; seedMat]; end
    initPts = struct('X', seedMat);

    % --- (1) 전역 탐색: surrogateopt (저정밀 fidelity) ---
    ray_nums_current = RAY_SEARCH;  wave_n_current = WAVE_N_SEARCH;
    sopts = optimoptions('surrogateopt', ...
        'MaxFunctionEvaluations', EVALS_PER_START, ...
        'MinSurrogatePoints',     MIN_SURR_POINTS, ...
        'InitialPoints',          initPts, ...
        'UseParallel', false, 'PlotFcn', [], 'Display', 'iter');
    [xS, fS, ~, outS] = surrogateopt(@surrogate_objconstr, lb, ub, sopts);

    surrEQE = NaN;
    if ~isempty(xS) && isValidPoints(xS(:).')
        surrEQE = -fS;  xS = xS(:).';
    else
        fprintf('[Warn] start %d: surrogateopt feasible 해 미반환.\n', st);  xS = [];
    end

    % --- (2) 국소 정련: patternsearch (저정밀 fidelity) ---
    xPol = []; polishEQE = NaN;  x0 = xS;
    if isempty(x0) && ~isempty(gBestX), x0 = gBestX; end
    if ~isempty(x0)
        fprintf('--- Local polish (patternsearch, %d evals) ---\n', POLISH_EVALS);
        ray_nums_current = RAY_SEARCH;  wave_n_current = WAVE_N_SEARCH;
        try
            xPol = patternsearch(@polish_objective, x0, [],[],[],[], lb, ub, [], psOpts);
            xPol = xPol(:).';  polishEQE = -polish_objective(xPol);
        catch perr
            fprintf('[Warn] patternsearch 실패(%s).\n', perr.message);  xPol = [];
        end
    end

    % --- (3) 최종 검증: full-fidelity 반복 평가 후 승자 채택 ---
    candX = {};
    if ~isempty(xS),  candX{end+1} = xS;  end
    if ~isempty(xPol) && (isempty(xS) || ~isequal(xPol, xS)), candX{end+1} = xPol; end
    if ~isempty(gBestX) && ~any(cellfun(@(c) isequal(c, gBestX), candX))
        candX{end+1} = gBestX;
    end

    ray_nums_current = RAY_FINAL;  wave_n_current = WAVE_N_FINAL;
    candMean = -inf(1, numel(candX));  candStd = zeros(1, numel(candX));
    for c = 1:numel(candX)
        if ~isValidPoints(candX{c}), continue; end
        e = nan(1, N_FINAL_REP);
        for rrep = 1:N_FINAL_REP, e(rrep) = simulate_EQE(candX{c}); end
        candMean(c) = mean(e, 'omitnan');  candStd(c) = std(e, 'omitnan');
        fprintf('  start %d cand %d: EQE_total = %.5f ± %.5f (N=%d, %d rays, full-λ)\n', ...
            st, c, candMean(c), candStd(c), N_FINAL_REP, RAY_FINAL);
    end
    [bestEQE, ci] = max(candMean);

    if isfinite(bestEQE) && bestEQE > gBestEQE
        gBestEQE = bestEQE;  gBestStd = candStd(ci);  gBestX = candX{ci};
        fprintf('  [Global] start %d 에서 전역 best 갱신: EQE_total = %.5f ± %.5f\n', ...
            st, gBestEQE, gBestStd);
    else
        fprintf('  [Global] start %d: 갱신 없음 (현 best = %.5f).\n', st, gBestEQE);
    end

    start_log(st).surrEQE=surrEQE; start_log(st).polishEQE=polishEQE;
    start_log(st).bestEQE=bestEQE; start_log(st).bestStd=candStd(min(ci,numel(candStd)));
    start_log(st).evals=outS.funccount;

    save('freeform_EQEtotal_result.mat', 'gBestX','gBestEQE','gBestStd', ...
        'varNames','start_log','lb','ub');
    fprintf('############ start %d done | 전역 best EQE_total = %.5f ± %.5f ############\n', ...
        st, gBestEQE, gBestStd);
end

%% ===== 결과 요약 =====
disp('=== Freeform EQE_total global optimization finished ===');
if isempty(gBestX)
    warning('feasible 해를 찾지 못했습니다. 시드/제약을 점검하세요.');
else
    bestT = array2table(gBestX, 'VariableNames', varNames);
    fprintf('\n######## Best freeform lens (EQE_total) ########\n');
    fprintf('  EQE_total = %.5f ± %.5f  (full-λ %d rays, N=%d)\n', ...
        gBestEQE, gBestStd, RAY_FINAL, N_FINAL_REP);
    disp('  best design variables:'); disp(bestT);
    save('freeform_EQEtotal_result.mat', 'gBestX','gBestEQE','gBestStd', ...
        'varNames','start_log','lb','ub','bestT');
end

figure('Name','freeform EQE_total multi-start','Color','w');
sEQE = arrayfun(@(s) s.bestEQE, start_log);
bar(1:NUM_STARTS, sEQE); grid on;
xlabel('global start #'); ylabel('verified best EQE\_total');
title('Freeform lens co-optimization: multi-start best (full-\lambda verified)');


%% ===== LightTools 1회 평가 공용 래퍼 =====
function eqe = simulate_EQE(pt)
global ID_swept ltml ltloc eval_count restart_interval
eval_count = eval_count + 1;
if mod(eval_count, restart_interval) == 0
    fprintf('\n[Refresh] 시뮬 %d회. LightTools 재시작...\n', eval_count);
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    pause(2);
end
try
    eqe = objFcn_EQEtotal(pt).EQE_total;
    if eqe == 0, eqe = NaN; end   % 파셋 불일치 등 기하 오류: 평가 실패로 처리
catch err
    fprintf('\n[Error] eval %d LightTools 충돌: %s\n', eval_count, err.message);
    eqe = NaN;
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
end
end

%% ===== surrogateopt 목적+제약 결합 함수 =====
function out = surrogate_objconstr(x)
x = x(:).';
if ~isValidPoints(x)
    out.Ineq = 1;  out.Fval = 1;  return;   % infeasible: 시뮬 없이 즉시 반환
end
e = simulate_EQE(x);
if ~isfinite(e)
    out.Ineq = 1;  out.Fval = 1;            % 시뮬 실패도 회피
else
    out.Ineq = -1; out.Fval = -e;           % EQE 최대화 == -EQE 최소화
end
end

%% ===== patternsearch 정련용 목적함수 =====
function f = polish_objective(x)
x = x(:).';
if ~isValidPoints(x), f = 0; return; end
e = simulate_EQE(x);
if isnan(e), e = 0; end
f = -e;
end

%% ===== 무작위 valid 시드 생성 (rejection sampling) =====
function P = genValidPoints(K, lb, ub)
dim = numel(lb);  P = zeros(K, dim);
for i = 1:K
    ok = false;
    while ~ok
        p = lb + rand(1, dim) .* (ub - lb);
        if isValidPoints(p), ok = true; P(i, :) = p; end
    end
end
end


%% Objective Function (EQE_total, multi-fidelity: ray 수 + 파장 step 가변)
function output = objFcn_EQEtotal(point)
global ID_LT ID_swept ltml ltloc count ray_nums_current wave_n_current
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

d_sub=1.295;
r_OLED=1;
x_pattern=15;
y_pattern=15;
Lensheight=0.01;
wavelength_start=453;
wavelength_end=753;

% multi-fidelity: 파장 step 과 ray 수를 단계별로 가변
if isempty(wave_n_current), n = 10;    else, n = wave_n_current;    end
if isempty(ray_nums_current), ray_nums = 10000; else, ray_nums = ray_nums_current; end

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
dETL = point(11); dHTL = point(12); stretchZ=point(13);

% Create spline control points
xy = zeros(7,2);
xy(1,:) = [0, 1];  xy(7,:) = [1, 0];
xy(2,:) = [x2, y2];  xy(3,:) = [x3, y3];  xy(4,:) = [x4, y4];
xy(5,:) = [x5, y5];  xy(6,:) = [x6, y6];

lt = ltloc.GetLTAPI(ID_swept);
ltx= getltpointer(ID_swept);
lt2 = ltloc.GetLTAPI(ID_LT);

Curve="LENS_MANAGER[1].COMPONENTS[Components].SWEPT_SOLID[SweptEntity].SWEPT_PRIMITIVE[SweptPrimitive].SWEPT_PROFILE[SweptProfile].FITTED_CURVE[SweptSurface_1]";
ltx.SetSweptProfilePoints(Curve,xy,7);
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
ltx.SetSweptProfilePoints(Curve,xy,7);
ltx.DbSet(Curve,'StartSlopeMode',"Auto");
ltx.DbSet(Curve,'EndSlopeMode',"Auto");

xy_l = zeros(7,2);
for j=1:7
    xy_l(j,1) = ltml.LTDbGet(lt, Key, 'YAt', j);
    xy_l(j,2) = ltml.LTDbGet(lt, Key, 'ZAt', j);
end

% [버그수정] float 완전일치(isequal) 대신 tolerance 비교
tol = 1e-4;
if max(abs(xy(:) - xy_l(:))) > tol
    output = struct('EQE_0_20',0,'EQE_20_40',0,'EQE_40_60',0,'EQE_60_80',0,'EQE_total',0);
    return;
end

% File name and path configuration
rng('shuffle')
strLength = 10;  charSet = ['a':'z' 'A':'Z' '0':'9'];
randIndices = randi(length(charSet), 1, strLength);
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

EML_position=4;  z0=12.5;  u_data_num=499;  max_u=3;

CPS_result=CPS_for_Isub(no_bar,ne_bar,thickness,emission_spectrum,eta_rad,horizontal_dipole_ratio,bottom_air_refractive_index,EML_position,z0,u_data_num,max_u,wavelength);
EQE_air_CPS=CPS_result.EQE_air; %#ok<NASGU>
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
Transmittance=(T_p_bottom+T_s_bottom)/2; %#ok<NASGU>
Reflectance=(R_p_bottom+R_s_bottom)/2;

%% Coating (.mat to .coa)
lt = ltloc.GetLTAPI(ID_LT);
fileID = fopen(sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count), 'w');
fprintf(fileID,'%s\n%s%d\n%s\n%s\n%s\n%s\n ','DFAT Version 1.0', 'DATANAME: R_Bottom_',count, 'ABSORBING: YES', 'INDEX: 1.51', 'DATAITEMS: TAVG RAVG');
for i=wavelength_start:wavelength_end
    fprintf(fileID,'%s  %d\n','wv',i);
    for j=0:89
        fprintf(fileID,'%s  %d  %d  %.3f\n', 'AOI',j, 0, Reflectance(i-wavelength_start+1,j+1));
    end
end
fclose(fileID);   % [버그수정] LightTools 가 읽기 전에 버퍼 플러시 + 잠금 해제

ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings]" LoadFileName="' sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count) '"']);
List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');
Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');
Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName',sprintf('R_Bottom_%d', count));

%%
I_white=0.5*(CPS_result.I_sub_s+CPS_result.I_sub_p);
sin089=sind(0:89);
P_white=I_white.*repmat(sin089,wavelength_num,1);
weight_factor=sum(P_white,2);
I_white_ang=sum(P_white); %#ok<NASGU>
wavelength_num=length(wavelength);

I_air_1_2=zeros(90,(wavelength_num+n-1)/n);
for wv=1:n:wavelength_num
    fileID = fopen('C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','w');
    fprintf(fileID,'%s  %d  %d  %d  %d  %d  %d','SPHEREMESH:',1, 90, 0, 0, 360, 90);
    writematrix(flip(I_white(wv,:).'),'C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','Delimiter','tab','WriteMode','append');
    fclose(fileID);
    SRList=ltml.LTDbList(lt, 'Lens_manager[1]','DISK_SOURCE');
    SRKey=ltml.LTListAtPos(lt,SRList,1);
    ltml.LTDbSet(lt,SRKey,'Radiant_Power', weight_factor(wv));
    for k=1:1
        SRList=ltml.LTDbList(lt, 'Lens_manager[1]','Spectral_region');
        SRKey=ltml.LTListAtPos(lt,SRList,k+1);
        ltml.LTDbSet(lt,SRKey,'Spectral_Definition', 'Monochromatic');
        ltml.LTDbSet(lt,SRKey,'Single_Wavelength', wv+wavelength_start-1);
        List=ltml.LTDbList(lt,'lens_manager[1]','DIRECTION_GRID_APODIZER');
        Key=ltml.LTListAtPos(lt,List,k);
        pathname='C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\';
        ltml.LTDbSet(lt,Key,'LoadFileName',[pathname sprintf('AI_temp.txt')]);
    end
    ltml.LTBegin(lt);
    ltml.LTCmd(lt,'\V3D BeginAllSimulations');
    ltml.LTEnd(lt);
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,1);
    Power_output(wv)=ltml.LTDbGet(lt,Key,'TotalPower');
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,2);
    Power_output_30(wv)=ltml.LTDbGet(lt,Key,'TotalPower'); %#ok<NASGU>
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,3);
    for j=1:90
        I_air_1_JH(91-j,:)=ltml.LTDbGet(lt,Key,'CellValue_UI',1,91-j);
    end
    I_air_1_2(:,(wv+n-1)/n)=smooth(I_air_1_JH);
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
EQE_wv_matrix = Power_output_2 ./ weight_factor_2;
EQE_sub_matrix_2 = EQE_sub_matrix_2 / sum(EQE_sub_matrix_2) * EQE_sub_CPS;
EQE_total = sum(EQE_wv_matrix .* EQE_sub_matrix_2);

% (참고용) 각도별 EQE 도 계산해 함께 반환
EQE_0_20=0; EQE_20_40=0; EQE_40_60=0; EQE_60_80=0;
sin_col = sin089(:);
for k = 1:K
    contrib_k = EQE_wv_matrix(k) * EQE_sub_matrix_2(k);
    I_theta = I_air_1_2(:,k);
    W_theta = I_theta .* sin_col;  W_tot = sum(W_theta);
    EQE_0_20  = EQE_0_20  + contrib_k * sum(W_theta(1:20))  / W_tot;
    EQE_20_40 = EQE_20_40 + contrib_k * sum(W_theta(21:40)) / W_tot;
    EQE_40_60 = EQE_40_60 + contrib_k * sum(W_theta(41:60)) / W_tot;
    EQE_60_80 = EQE_60_80 + contrib_k * sum(W_theta(61:80)) / W_tot;
end

output = struct('EQE_0_20',EQE_0_20,'EQE_20_40',EQE_20_40, ...
    'EQE_40_60',EQE_40_60,'EQE_60_80',EQE_60_80,'EQE_total',EQE_total);

List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');
Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');
Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName','R_temp');
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings].COATING[' sprintf('R_Bottom_%d', count) ']" Delete= \Q']);
fclose('all');
end

%% Spline Constraints Function (원본 그대로)
function TF = isValidPoints(X)
numRows = size(X,1);  numPts = 7;  TF = true(numRows,1);
for k = 1:numRows
    x = [0, X(k,1:5), 1];
    y = [1, X(k,6:10), 0];
    violates = false;
    for i = 1:numPts - 1
        for j = i + 2:numPts - 1
            if i == 1 && j == numPts - 1, continue; end
            if checkIntersection([x(i),y(i)],[x(i+1),y(i+1)],[x(j),y(j)],[x(j+1),y(j+1)])
                violates = true; break;
            end
        end
        if violates, break; end
    end
    if ~violates
        for i = 1:numPts - 2
            if isCollinear([x(i),y(i)],[x(i+1),y(i+1)],[x(i+2),y(i+2)])
                violates = true; break;
            end
        end
    end
    if ~violates
        minD = 0.05; maxD = 1.0;
        d = hypot(diff(x), diff(y));
        if any(d < minD | d > maxD), violates = true; end
    end
    if ~violates
        maxAng = 2 * pi / 3;
        for i = 2:numPts - 1
            v1 = [x(i),y(i)] - [x(i-1),y(i-1)];
            v2 = [x(i+1),y(i+1)] - [x(i),y(i)];
            ang = atan2(norm(cross([v1,0],[v2,0])), dot(v1,v2));
            if ang > maxAng, violates = true; break; end
        end
    end
    TF(k) = ~violates;
end
    function isCol = isCollinear(p1, p2, p3)
        area = 0.5 * det([p1 1; p2 1; p3 1]);
        isCol = abs(area) < 1e-5;
    end
    function intersects = checkIntersection(p1, p2, p3, p4)
        function o = orientation(p, q, r)
            o = (q(2)-p(2))*(r(1)-q(1)) - (q(1)-p(1))*(r(2)-q(2));
        end
        o1 = orientation(p1,p2,p3);  o2 = orientation(p1,p2,p4);
        o3 = orientation(p3,p4,p1);  o4 = orientation(p3,p4,p2);
        intersects = (o1*o2 < 0) && (o3*o4 < 0);
    end
end


function RenewLightTools()
global ID_LT ID_swept ltml ltloc lt
lt_exe_path = 'C:\Program Files\Optical Research Associates\LightTools 2023.03\lt.exe';
model_file_path_swept = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\SweptEntity.2.lts';
model_file_path_LT = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\Lens_size_effect_for_PSO_bump_modified_v1.1.lts';

fprintf('--- Restarting LightTools ---\n');
target_user = 'jhkim';
kill_cmd = sprintf('taskkill /F /FI "USERNAME eq %s" /IM lt.exe', target_user);
[~, ~] = system(kill_cmd);
pause(2);

cmd = sprintf('"%s" "%s" &', lt_exe_path, model_file_path_swept);
status = system(cmd); %#ok<NASGU>
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
        ID_swept = str2double(tokens{1}{1});
        fprintf('PID found for user %s: %d\n', target_user, ID_swept);
    else
        error('프로세스는 찾았으나 PID 추출 실패.');
    end
else
    error('사용자 %s 로 실행된 LightTools 를 찾을 수 없습니다.', target_user);
end
cmd = sprintf('"%s" "%s" &', lt_exe_path, model_file_path_LT);
status = system(cmd); %#ok<NASGU>
find_cmd = sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);
[status, cmdout] = system(find_cmd);
if status == 0 && contains(cmdout, 'lt.exe')
    tokens = regexp(cmdout, '"(\d+)"', 'tokens');
    if ~isempty(tokens)
        ID_LT = str2double(tokens{3}{1});
        fprintf('PID found for user %s: %d\n', target_user, ID_LT);
    else
        error('프로세스는 찾았으나 PID 추출 실패.');
    end
else
    error('사용자 %s 로 실행된 LightTools 를 찾을 수 없습니다.', target_user);
end
pause(5);
end

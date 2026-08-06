% ============================================================
%  pareto_front_freeform.m
%
%  목적: "효율(EQE_total)과 방향성(EQE_40_60)이 교환될 뿐인가"를 판정하기 위한
%        Pareto front 추적 + achievable region 산점도 데이터 수집.
%
%  [왜 이 코드인가]
%   - 단일 목적 최적화 두 번(각 극단)만으로는 트레이드오프 '곡선의 모양'을 모른다.
%     가중합 w 를 스윕하며 최적화하면 Pareto front 가 직접 얻어지고,
%       곡선이 오목하게 휘면 -> 트레이드오프 실재 (원고 2.3절 유지)
%       한 점에 뭉치면     -> 두 목표가 정렬됨 (2.3절 '교환' 주장 철회)
%   - 무작위 valid 설계는 frontier 를 못 찾지만(13차원, 차원의 저주)
%     achievable region '내부'를 편향 없이 채우는 데 필수. 둘을 함께 수집한다.
%
%  [수집물]  EVAL_LOG = [x(1:13), EQE_total, EQE_40_60, phase, w]
%     phase 1 = random sampling, 2 = weighted-sum optimization, 3 = 고정밀 재평가
%     -> 산점도(내부+포락선), Pareto front, 형상 다양성 비교에 모두 사용.
%
%  [주의] 밀도를 확률로 해석하지 말 것: phase 2 점들은 최적점 근처에 몰린다.
%         내부 밀도 해석은 phase 1(무작위)만 사용.
%
%  기반: opt_freeform_EQEtotal_fast_v1.m (동일한 LightTools 연동/기하/스택)
% ============================================================
clear;
%% For LightTools Connection
global ID_swept ID_LT ltml ltloc count eval_count restart_interval ...
       ray_nums_current wave_n_current EVAL_LOG EVAL_PHASE EVAL_W
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

%% ===== Multi-fidelity =====
% [주의] 파장 step 은 objFcn_both 안의 파장창 폭에 맞춰 정해야 한다.
%   넓은 창(453-753, 301개): SEARCH=30 -> K=11,  FINAL=10 -> K=31
%   좁은 창(593-603,  11개): SEARCH=10 -> K=2,   FINAL=2  -> K=6
%   * 인덱스는 wv_list = 1:n:wavelength_num 로 생성되므로 나누어떨어지지 않아도
%     오류는 나지 않지만, n 이 파장 개수보다 크면 K=1 이 되어 스펙트럼 정보가 사라진다.
WAVE_N_SEARCH = 10;      % 탐색용 파장 step
WAVE_N_FINAL  = 2;       % 검증용 파장 step
RAY_SEARCH    = 10000;
RAY_FINAL     = 50000;
N_FINAL_REP   = 3;

%% ===== Pareto 스윕 설정 =====
W_LIST          = [0, 0.25, 0.5, 0.75, 1.0];   % w=1 -> EQE_total만, w=0 -> EQE_40_60만
EVALS_PER_W     = 120;   % w 하나당 surrogateopt 평가 예산
MIN_SURR_POINTS = 25;
N_SEED_VALID    = 30;    % w 하나당 valid 초기 시드
POLISH_EVALS    = 15;

%% ===== 무작위 표본 (achievable region 내부 채우기) =====
N_RANDOM = 150;          % 무작위 valid 설계 개수 (비용 감안해 조정)
                         % * 이 점들만 밀도 해석에 사용 가능

%% ===== 정규화 상수 =====
% 두 목적의 스케일이 크게 다르므로(EQE_total ~0.5, EQE_40_60 ~0.03) 정규화 필수.
% 아래 값은 기존 결과 기반 초기값이며, phase 1(무작위) 결과로 자동 갱신된다.
REF_TOTAL = 0.55;
REF_BAND  = 0.035;
AUTO_CALIBRATE = true;   % phase 1 최댓값의 1.1배로 REF 갱신

%% Optimization Variables (13-dim, 기존과 동일)
varNames = {'x2','x3','x4','x5','x6', 'y2','y3','y4','y5','y6', 'dETL','dHTL','stretchZ'};
lb = [0, 0, 0, 0, 0, 0,   0,   0,   0,   0,   10, 10, 0.1];
ub = [1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 1.5, 1.5, 150,150, 3];
nvar = numel(lb);

EVAL_LOG = [];           % [x(1:13), EQE_total, EQE_4060, phase, w]

psOpts = optimoptions('patternsearch', ...
    'MaxFunctionEvaluations', POLISH_EVALS, ...
    'InitialMeshSize', 0.1, 'MeshTolerance', 1e-3, ...
    'Cache', 'on', 'Display', 'off');

%% =====================================================================
%  PHASE 1 — 무작위 valid 설계 (achievable region 내부, 편향 없음)
%% =====================================================================
fprintf('\n########## PHASE 1: random valid designs (N=%d) ##########\n', N_RANDOM);
EVAL_PHASE = 1;  EVAL_W = NaN;
ray_nums_current = RAY_SEARCH;  wave_n_current = WAVE_N_SEARCH;
eval_count = 0;
Prand = genValidPoints(N_RANDOM, lb, ub);
for i = 1:N_RANDOM
    [et, eb] = simulate_both(Prand(i,:));
    if mod(i,10)==0
        fprintf('  random %3d/%d : EQE_total=%.4f  EQE_4060=%.4f\n', i, N_RANDOM, et, eb);
    end
end
save('pareto_log_partial.mat','EVAL_LOG','varNames','lb','ub');

% 정규화 상수 자동 보정
if AUTO_CALIBRATE && ~isempty(EVAL_LOG)
    mt = max(EVAL_LOG(:,nvar+1), [], 'omitnan');
    mb = max(EVAL_LOG(:,nvar+2), [], 'omitnan');
    if isfinite(mt) && mt>0, REF_TOTAL = 1.1*mt; end
    if isfinite(mb) && mb>0, REF_BAND  = 1.1*mb; end
end
fprintf('  정규화: REF_TOTAL=%.4f, REF_BAND=%.4f\n', REF_TOTAL, REF_BAND);

%% =====================================================================
%  PHASE 2 — 가중합 스윕으로 Pareto front 추적
%% =====================================================================
nW = numel(W_LIST);
pareto_x    = nan(nW, nvar);
pareto_tot  = nan(nW,1);
pareto_band = nan(nW,1);

for k = 1:nW
    w = W_LIST(k);
    fprintf('\n########## PHASE 2: w = %.2f  (%d/%d) ##########\n', w, k, nW);
    fprintf('  목적 = %.2f*EQE_total/%.3f + %.2f*EQE_4060/%.4f\n', ...
        w, REF_TOTAL, 1-w, REF_BAND);

    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);
    ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    eval_count = 0;
    EVAL_PHASE = 2;  EVAL_W = w;
    ray_nums_current = RAY_SEARCH;  wave_n_current = WAVE_N_SEARCH;

    seedMat = genValidPoints(N_SEED_VALID, lb, ub);
    sopts = optimoptions('surrogateopt', ...
        'MaxFunctionEvaluations', EVALS_PER_W, ...
        'MinSurrogatePoints',     MIN_SURR_POINTS, ...
        'InitialPoints',          struct('X', seedMat), ...
        'UseParallel', false, 'PlotFcn', [], 'Display', 'iter');
    [xS, ~] = surrogateopt(@(x) scalar_objconstr(x, w, REF_TOTAL, REF_BAND), lb, ub, sopts);

    % 국소 정련
    if ~isempty(xS) && isValidPoints(xS(:).')
        x0 = xS(:).';
        try
            xP = patternsearch(@(x) scalar_polish(x, w, REF_TOTAL, REF_BAND), ...
                               x0, [],[],[],[], lb, ub, [], psOpts);
            xP = xP(:).';
        catch
            xP = x0;
        end
    else
        fprintf('  [Warn] w=%.2f: feasible 해 미반환\n', w);
        continue;
    end

    % 후보 중 가중합이 더 좋은 쪽 채택 (저정밀 기준)
    cands = {x0};
    if ~isequal(xP, x0), cands{end+1} = xP; end
    bestScore = -inf; bestX = [];
    for c = 1:numel(cands)
        [et, eb] = simulate_both(cands{c});
        if ~isfinite(et), continue; end
        sc = w*et/REF_TOTAL + (1-w)*eb/REF_BAND;
        if sc > bestScore, bestScore = sc; bestX = cands{c}; end
    end
    if isempty(bestX), continue; end

    % 고정밀 재평가
    EVAL_PHASE = 3;
    ray_nums_current = RAY_FINAL;  wave_n_current = WAVE_N_FINAL;
    et_r = nan(1,N_FINAL_REP);  eb_r = nan(1,N_FINAL_REP);
    for r = 1:N_FINAL_REP
        [et_r(r), eb_r(r)] = simulate_both(bestX);
    end
    pareto_x(k,:)   = bestX;
    pareto_tot(k)   = mean(et_r,'omitnan');
    pareto_band(k)  = mean(eb_r,'omitnan');
    fprintf('  >>> w=%.2f : EQE_total=%.5f  EQE_4060=%.5f  (선택성 %.1f%%)\n', ...
        w, pareto_tot(k), pareto_band(k), 100*pareto_band(k)/pareto_tot(k));

    save('pareto_front_result.mat','EVAL_LOG','pareto_x','pareto_tot','pareto_band', ...
         'W_LIST','varNames','lb','ub','REF_TOTAL','REF_BAND');
end

%% =====================================================================
%  결과 정리 + 판정
%% =====================================================================
fprintf('\n################ PARETO FRONT ################\n');
fprintf('%6s | %12s | %12s | %10s\n','w','EQE_total','EQE_40_60','선택성');
for k = 1:nW
    if isfinite(pareto_tot(k))
        fprintf('%6.2f | %12.5f | %12.5f | %9.1f%%\n', W_LIST(k), ...
            pareto_tot(k), pareto_band(k), 100*pareto_band(k)/pareto_tot(k));
    end
end

% Test 1 — 트레이드오프 존재 판정
iT = find(W_LIST==1, 1);  iB = find(W_LIST==0, 1);
if ~isempty(iT) && ~isempty(iB) && isfinite(pareto_tot(iT)) && isfinite(pareto_tot(iB))
    dropTot  = 1 - pareto_tot(iB)/pareto_tot(iT);     % 방향성 최적화 시 효율 손실
    gainBand = pareto_band(iB)/pareto_band(iT) - 1;   % 그 대가로 얻은 방향성 이득
    fprintf('\n[Test 1] 방향성 최적화(w=0)로 전환 시:\n');
    fprintf('   EQE_total  %.5f -> %.5f  (%+.1f%%)\n', ...
        pareto_tot(iT), pareto_tot(iB), -100*dropTot);
    fprintf('   EQE_40_60  %.5f -> %.5f  (%+.1f%%)\n', ...
        pareto_band(iT), pareto_band(iB), 100*gainBand);
    if dropTot > 0.05
        fprintf('   => 트레이드오프 실재 (효율 %.1f%% 희생). 원고 2.3절 유지.\n', 100*dropTot);
    else
        fprintf('   => 트레이드오프 미미 (효율 손실 %.1f%%). 2.3절 "교환" 주장 재검토 필요.\n', 100*dropTot);
    end
end

save('pareto_front_result.mat','EVAL_LOG','pareto_x','pareto_tot','pareto_band', ...
     'W_LIST','varNames','lb','ub','REF_TOTAL','REF_BAND');
fprintf('\nsaved -> pareto_front_result.mat  (EVAL_LOG %d points)\n', size(EVAL_LOG,1));

%% ===== 그림 =====
% EVAL_LOG = [x(1:nvar) | EQE_total | b0_20 | b20_40 | b40_60 | b60_80 | phase | w]
nv = nvar;
Et = EVAL_LOG(:,nv+1);          % EQE_total
Eb = EVAL_LOG(:,nv+4);          % EQE_40_60 (네 구간 중 세 번째)
Ph = EVAL_LOG(:,nv+6);          % phase
ok = isfinite(Et) & isfinite(Eb) & Et>0;

figure('Name','Achievable region & Pareto front','Color','w','Position',[100 100 1000 420]);

subplot(1,2,1);
scatter(Et(ok&Ph==1), Eb(ok&Ph==1), 14, [.6 .6 .6], 'filled'); hold on;
scatter(Et(ok&Ph==2), Eb(ok&Ph==2), 14, [.2 .5 .8], 'filled');
plot(pareto_tot, pareto_band, 'r-o', 'LineWidth', 2, 'MarkerFaceColor','r');
xlabel('EQE_{total}'); ylabel('EQE_{40-60}'); grid on;
legend({'random (내부)','optimization','Pareto front'},'Location','best','FontSize',8);
title('(a) achievable region');

subplot(1,2,2);
Sel = Eb./Et;
scatter(Et(ok&Ph==1), Sel(ok&Ph==1), 14, [.6 .6 .6], 'filled'); hold on;
scatter(Et(ok&Ph==2), Sel(ok&Ph==2), 14, [.2 .5 .8], 'filled');
plot(pareto_tot, pareto_band./pareto_tot, 'r-o','LineWidth',2,'MarkerFaceColor','r');
xlabel('EQE_{total}'); ylabel('selectivity  EQE_{40-60}/EQE_{total}'); grid on;
title('(b) efficiency vs selectivity');
saveas(gcf,'pareto_front.png');
fprintf('saved -> pareto_front.png\n');


%% ===== 평가 래퍼 (네 각도구간을 모두 로그에 기록) =====
%  EVAL_LOG 열 구성:
%    [ x(1:13) | EQE_total | EQE_0_20 | EQE_20_40 | EQE_40_60 | EQE_60_80 | phase | w ]
%  -> 한 번의 시뮬로 네 구간 선택성을 동시에 얻어, Fig 3(FoM 지도)에서
%     S = sin^2(th2) - sin^2(th1) 예측을 네 구간 독립 검증할 수 있다.
function [eqe_total, eqe_band] = simulate_both(pt)
global ID_swept ltml ltloc eval_count restart_interval EVAL_LOG EVAL_PHASE EVAL_W
eval_count = eval_count + 1;
if mod(eval_count, restart_interval) == 0
    fprintf('\n[Refresh] 시뮬 %d회. LightTools 재시작...\n', eval_count);
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    pause(2);
end
bins = [NaN NaN NaN NaN];
try
    out = objFcn_both(pt);
    eqe_total = out.EQE_total;
    eqe_band  = out.EQE_40_60;
    bins = [out.EQE_0_20, out.EQE_20_40, out.EQE_40_60, out.EQE_60_80];
    if eqe_total == 0
        eqe_total = NaN; eqe_band = NaN; bins = [NaN NaN NaN NaN];
    end
catch err
    fprintf('\n[Error] eval %d LightTools 충돌: %s\n', eval_count, err.message);
    eqe_total = NaN;  eqe_band = NaN;
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
end
EVAL_LOG(end+1,:) = [pt(:).', eqe_total, bins, EVAL_PHASE, EVAL_W];
end

%% ===== 가중합 목적 (surrogateopt: 제약 결합형) =====
function out = scalar_objconstr(x, w, refT, refB)
x = x(:).';
if ~isValidPoints(x)
    out.Ineq = 1;  out.Fval = 1;  return;   % infeasible: 시뮬 없이 반환
end
[et, eb] = simulate_both(x);
if ~isfinite(et) || ~isfinite(eb)
    out.Ineq = 1;  out.Fval = 1;
else
    out.Ineq = -1;
    out.Fval = -( w*et/refT + (1-w)*eb/refB );   % 최대화 -> 부호 반전
end
end

%% ===== 가중합 목적 (patternsearch) =====
function f = scalar_polish(x, w, refT, refB)
x = x(:).';
if ~isValidPoints(x), f = 0; return; end
[et, eb] = simulate_both(x);
if ~isfinite(et) || ~isfinite(eb), f = 0; return; end
f = -( w*et/refT + (1-w)*eb/refB );
end

%% ===== 무작위 valid 시드 생성 =====
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


%% ===== Objective (EQE_total 과 EQE_40_60 동시 산출) =====
function output = objFcn_both(point)
global ID_LT ID_swept ltml ltloc count ray_nums_current wave_n_current
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

d_sub=1.295;  r_OLED=1;  x_pattern=15;  y_pattern=15;  Lensheight=0.01;
wavelength_start=453;  wavelength_end=753;

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

x2 = point(1);  x3 = point(2);  x4 = point(3);  x5 = point(4);  x6 = point(5);
y2 = point(6);  y3 = point(7);  y4 = point(8);  y5 = point(9);  y6 = point(10);
dETL = point(11); dHTL = point(12); stretchZ=point(13);

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
tol = 1e-4;   % float 완전일치(isequal) 대신 tolerance 비교
if max(abs(xy(:) - xy_l(:))) > tol
    output = struct('EQE_0_20',0,'EQE_20_40',0,'EQE_40_60',0,'EQE_60_80',0,'EQE_total',0);
    return;
end

rng('shuffle')
charSet = ['a':'z' 'A':'Z' '0':'9'];
index = charSet(randi(length(charSet), 1, 10));
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

%% CPS
load('nk_JH33.mat');  load('Photopic_400_800.mat');
load('CIE_1931.mat'); load('R_pd.mat');
wavelength=(wavelength_start:wavelength_end).';
wavelength_num=length(wavelength);
emission_spectrum=spectrum.l_I_Irdmppyph2tmd(wavelength_start-399:wavelength_end-399,:);
eta_rad=0.98;  horizontal_dipole_ratio=0.865;
bottom_air_refractive_index=ones(wavelength_num,1);

no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
layer_num=size(no_bar,2);
sin089=sind(0:89);  cos089=cosd(0:89);
no_bar=no_bar(wavelength_start-399:wavelength_end-399,:);
ne_bar=ne_bar(wavelength_start-399:wavelength_end-399,:);
thickness=[100 dETL 25 10 dHTL 150];
EML_position=4;  z0=12.5;  u_data_num=499;  max_u=3;

CPS_result=CPS_for_Isub(no_bar,ne_bar,thickness,emission_spectrum,eta_rad,horizontal_dipole_ratio,bottom_air_refractive_index,EML_position,z0,u_data_num,max_u,wavelength);
EQE_sub_CPS=CPS_result.EQE_sub;

%% bottom reflectance
TMF_OLED_bottom_p=TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_OLED_bottom_s=TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
R_p_bottom=abs(TMF_OLED_bottom_p.r_p).^2;
R_s_bottom=abs(TMF_OLED_bottom_s.r_s).^2;
Reflectance=(R_p_bottom+R_s_bottom)/2;

%% Coating
lt = ltloc.GetLTAPI(ID_LT);
fileID = fopen(sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count), 'w');
fprintf(fileID,'%s\n%s%d\n%s\n%s\n%s\n%s\n ','DFAT Version 1.0', 'DATANAME: R_Bottom_',count, 'ABSORBING: YES', 'INDEX: 1.51', 'DATAITEMS: TAVG RAVG');
for i=wavelength_start:wavelength_end
    fprintf(fileID,'%s  %d\n','wv',i);
    for j=0:89
        fprintf(fileID,'%s  %d  %d  %.3f\n', 'AOI',j, 0, Reflectance(i-wavelength_start+1,j+1));
    end
end
fclose(fileID);   % LightTools 가 읽기 전에 플러시

ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings]" LoadFileName="' sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count) '"']);
List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');
Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');
Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName',sprintf('R_Bottom_%d', count));

%% 파장 루프
I_white=0.5*(CPS_result.I_sub_s+CPS_result.I_sub_p);
sin089=sind(0:89);
P_white=I_white.*repmat(sin089,wavelength_num,1);
weight_factor=sum(P_white,2);
% [수정] 파장 샘플 인덱스를 명시적으로 생성.
%   기존 (wavelength_num+n-1)/n 식은 (wavelength_num-1)이 n으로 나누어떨어질 때만
%   정수가 되어, 좁은 파장창(예: 593-603, span 10)에 n=30 을 쓰면 zeros() 에서
%   "크기 입력값은 정수여야 합니다" 오류가 났다. 아래처럼 인덱스 벡터를 쓰면
%   파장창과 n 의 조합에 무관하게 안전하다.
wv_list = 1:n:wavelength_num;
K = numel(wv_list);
I_air_1_2 = zeros(90, K);
Power_output = zeros(1, wavelength_num);
for kk = 1:K
    wv = wv_list(kk);
    fileID = fopen('C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','w');
    fprintf(fileID,'%s  %d  %d  %d  %d  %d  %d','SPHEREMESH:',1, 90, 0, 0, 360, 90);
    writematrix(flip(I_white(wv,:).'),'C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','Delimiter','tab','WriteMode','append');
    fclose(fileID);
    SRList=ltml.LTDbList(lt, 'Lens_manager[1]','DISK_SOURCE');
    SRKey=ltml.LTListAtPos(lt,SRList,1);
    ltml.LTDbSet(lt,SRKey,'Radiant_Power', weight_factor(wv));
    SRList=ltml.LTDbList(lt, 'Lens_manager[1]','Spectral_region');
    SRKey=ltml.LTListAtPos(lt,SRList,2);
    ltml.LTDbSet(lt,SRKey,'Spectral_Definition', 'Monochromatic');
    ltml.LTDbSet(lt,SRKey,'Single_Wavelength', wv+wavelength_start-1);
    List=ltml.LTDbList(lt,'lens_manager[1]','DIRECTION_GRID_APODIZER');
    Key=ltml.LTListAtPos(lt,List,1);
    ltml.LTDbSet(lt,Key,'LoadFileName',['C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\' sprintf('AI_temp.txt')]);

    ltml.LTBegin(lt);
    ltml.LTCmd(lt,'\V3D BeginAllSimulations');
    ltml.LTEnd(lt);
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,1);
    Power_output(wv)=ltml.LTDbGet(lt,Key,'TotalPower');
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,3);
    for j=1:90
        I_air_1_JH(91-j,:)=ltml.LTDbGet(lt,Key,'CellValue_UI',1,91-j);
    end
    I_air_1_2(:,kk)=smooth(I_air_1_JH);
end

weight_factor_2  = zeros(K,1);
Power_output_2   = zeros(K,1);
EQE_sub_matrix_2 = zeros(K,1);
for k = 1:K
    idx = wv_list(k);
    weight_factor_2(k)  = weight_factor(idx);
    Power_output_2(k)   = Power_output(idx);
    EQE_sub_matrix_2(k) = CPS_result.EQE_sub_matrix(idx);
end
EQE_wv_matrix = Power_output_2 ./ weight_factor_2;
EQE_sub_matrix_2 = EQE_sub_matrix_2 / sum(EQE_sub_matrix_2) * EQE_sub_CPS;
EQE_total = sum(EQE_wv_matrix .* EQE_sub_matrix_2);

EQE_0_20=0; EQE_20_40=0; EQE_40_60=0; EQE_60_80=0;
sin_col = sin089(:);
for k = 1:K
    contrib_k = EQE_wv_matrix(k) * EQE_sub_matrix_2(k);
    W_theta = I_air_1_2(:,k) .* sin_col;  W_tot = sum(W_theta);
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

%% ===== Spline 제약 (기존과 동일) =====
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

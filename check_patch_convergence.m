% ============================================================
%  check_patch_convergence.m
%
%  목적: "25 x 25 mm 패치가 충분히 큰가" 에 대한 수렴 증거 확보.
%
%  [왜 필요한가] patch 15 캠페인과 25 캠페인의 EQE_total 이 다르다는 것은
%    15 mm 패치 가장자리 밖으로 기판광이 새어 나갔다는 뜻이다. 그렇다면
%    심사자는 당연히 "25 는 충분한가, 35 로 늘리면 또 오르는 것 아닌가" 를
%    묻는다. 임계각 광선의 왕복 측면 변위가 2.32 mm 이고 광원 반경이 1 mm
%    이므로 25 mm 면 왕복 ~5회 이상을 담지만, 이 산술보다 직접 측정이 강하다.
%
%  [무엇을 하는가] 최고 총 EQE freeform 설계 하나를 고정하고, 패치 크기만
%    15 / 25 / 35 mm 로 바꿔 고정밀 N_REP 회씩 잰다. 탐색 없음.
%    25 -> 35 변화가 Monte-Carlo 산포(~0.5%) 이내면 수렴 진술을 원고에 넣는다.
%
%  출력: patch_convergence_result.mat + 콘솔 표
%  실행 시간: 고정밀 3 패치 x N_REP 회. N_REP=3 이면 약 2 h.
%
%  기반: opt_4band_freeform.m 헬퍼 (동일 기하/스택/정밀도)
% ============================================================
clear;

DIARY_FILE = sprintf('patchconv_log_%s.txt', datestr(now,'yyyymmdd_HHMMSS'));
diary(DIARY_FILE); diary on;
fprintf('[Log] -> %s\n', DIARY_FILE);

%% For LightTools Connection
global ID_swept ID_LT ltml ltloc count eval_count restart_interval ...
       ray_nums_current wave_n_current EVAL_LOG EVAL_PHASE EVAL_W ...
       GEOM_TOL GEOM_MISMATCH_LOG REQUIRE_MONOTONIC_X PATCH_XY

GEOM_TOL = 1e-4;  GEOM_MISMATCH_LOG = [];
REQUIRE_MONOTONIC_X = true;

RenewLightTools();
try
    ltml.LTCmd(ltml.GetLTAPI(ID_LT), 'Message "Check Connection"');
catch
    ltml = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
end
count = 1;  restart_interval = 20;
lt = ltloc.GetLTAPI(ID_swept);
ltx= getltpointer(ID_swept);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

%% ===== 설정 =====
PATCHES = [15 25 35];      % mm. 25 가 기준, 15 는 과소(트렌드 앵커), 35 는 수렴 검사
N_REP   = 3;
RAY_FINAL = 50000;  WAVE_N_FINAL = 2;
ray_nums_current = RAY_FINAL;  wave_n_current = WAVE_N_FINAL;
EVAL_LOG = [];  EVAL_PHASE = 11;  EVAL_W = -2;
eval_count = 0;

%% ===== 시험 설계: 가중합 스윕의 고정밀 최고 총 EQE 설계 =====
Dp = load('pareto_front_result.mat');
[~, ib] = max(Dp.pareto_tot);
xTest = Dp.pareto_x(ib, :);
fprintf('시험 설계: pareto w-스윕 %d번째 (고정밀 EQE_total = %.4f @ patch 25)\n', ...
        ib, Dp.pareto_tot(ib));

%% ===== 패치 스윕 =====
Et  = nan(numel(PATCHES), N_REP);
Bn  = nan(numel(PATCHES), 4, N_REP);
for ip = 1:numel(PATCHES)
    PATCH_XY = PATCHES(ip);           % objFcn_both 가 전역으로 읽는다
    fprintf('\n--- patch %d x %d mm ---\n', PATCH_XY, PATCH_XY);
    for r = 1:N_REP
        [Et(ip,r), Bn(ip,:,r)] = simulate_bands(xTest);
        fprintf('  rep %d/%d : EQE_total = %.5f\n', r, N_REP, Et(ip,r));
        save('patch_convergence_result.mat','PATCHES','Et','Bn','xTest','N_REP','EVAL_LOG');
    end
end
PATCH_XY = 25;   % 복원 (이후 다른 스크립트 대비)

%% ===== 요약 =====
fprintf('\n%-8s | %10s | %8s | %s\n','patch','EQE_total','+-sd','상대변화(25 기준)');
m25 = mean(Et(PATCHES==25,:),'omitnan');
for ip = 1:numel(PATCHES)
    mi = mean(Et(ip,:),'omitnan');  si = std(Et(ip,:),'omitnan');
    fprintf('%-8d | %10.5f | %8.5f | %+7.2f%%\n', PATCHES(ip), mi, si, 100*(mi-m25)/m25);
end
i35 = find(PATCHES==35,1);
if ~isempty(i35)
    d35 = 100*(mean(Et(i35,:),'omitnan') - m25)/m25;
    if abs(d35) < 0.5
        fprintf(['\n=> 25 -> 35 변화 %+.2f%% 로 Monte-Carlo 산포(~0.5%%) 이내.\n' ...
                 '   "패치 크기 수렴" 진술을 원고 4.2절에 추가할 수 있다.\n'], d35);
    else
        fprintf(['\n=> 25 -> 35 변화 %+.2f%% 로 산포를 넘는다.\n' ...
                 '   25 패치의 EQE 는 하한으로 서술해야 하며, 결론(순위/비율)에는\n' ...
                 '   영향이 없는지 band 구성 변화를 함께 확인할 것.\n'], d35);
    end
end
S = squeeze(mean(Bn,3)) ./ mean(Et,2);   % 패치별 선택성
fprintf('\n선택성 (0-20/20-40/40-60/60-80):\n');
for ip = 1:numel(PATCHES)
    fprintf('  patch %-3d : %.3f %.3f %.3f %.3f\n', PATCHES(ip), S(ip,:));
end
save('patch_convergence_result.mat','PATCHES','Et','Bn','S','xTest','N_REP','EVAL_LOG');
fprintf('saved -> patch_convergence_result.mat\n');
diary off;

%% =====================================================================
%  이하 헬퍼 — opt_4band_freeform.m 과 동일 (MATLAB local function 규칙상 복제)
%% =====================================================================
function [eqe_total, bins] = simulate_bands(pt)
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
    bins = [out.EQE_0_20, out.EQE_20_40, out.EQE_40_60, out.EQE_60_80];
    if eqe_total == 0
        eqe_total = NaN; bins = [NaN NaN NaN NaN];
    end
catch err
    fprintf('\n[Error] eval %d LightTools 충돌: %s\n', eval_count, err.message);
    eqe_total = NaN;
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
end
EVAL_LOG(end+1,:) = [pt(:).', eqe_total, bins, EVAL_PHASE, EVAL_W];
end

%% ===== 구간 단일목적 (surrogateopt: 제약 결합형) =====
%  pareto 의 scalar_objconstr 에서 가중합 대신 단일 구간 열(bcol)만 최대화.
function out = band_objconstr(x, bcol, refB)
global REQUIRE_MONOTONIC_X
x = x(:).';
% 비단조 x 는 LightTools 표현 불가로 어차피 거부되므로, 제약 단계에서 배제해
% surrogateopt 가 예산을 NaN 에 쓰지 않게 한다.
if ~isempty(REQUIRE_MONOTONIC_X) && REQUIRE_MONOTONIC_X && any(diff(x(1:5)) < 0)
    out.Ineq = 1;  out.Fval = 1;  return;
end
if ~isValidPoints(x)
    out.Ineq = 1;  out.Fval = 1;  return;   % infeasible: 시뮬 없이 반환
end
[et, bins] = simulate_bands(x);
% bcol = 0 -> 총 EQE 대조군 (control arm). 그 외에는 해당 구간 EQE.
if bcol == 0, val = et; else, val = bins(bcol); end
if ~isfinite(et) || ~isfinite(val)
    out.Ineq = 1;  out.Fval = 1;
else
    out.Ineq = -1;
    out.Fval = -val/refB;   % 최대화 -> 부호 반전
end
end

%% ===== 구간 단일목적 (patternsearch) =====
function f = band_polish(x, bcol, refB)
global REQUIRE_MONOTONIC_X
x = x(:).';
if ~isempty(REQUIRE_MONOTONIC_X) && REQUIRE_MONOTONIC_X && any(diff(x(1:5)) < 0)
    f = 0; return;
end
if ~isValidPoints(x), f = 0; return; end
[et, bins] = simulate_bands(x);
if bcol == 0, val = et; else, val = bins(bcol); end
if ~isfinite(et) || ~isfinite(val), f = 0; return; end
f = -val/refB;
end

%% ===== 기하 거부 진단 (pareto_front_freeform.m 과 동일) =====
%  NaN 이 자주 뜨는 원인이 '수치 오차'인지 '형상 왜곡'인지 분포로 판별한다.
function report_geom_rejection(mismLog, tol)
if isempty(mismLog), return; end
mism = mismLog(:,1);  maxlen = mismLog(:,2);  resc = mismLog(:,3) > 0;
rej = mism > tol;
fprintf('\n--- 기하 거부 진단 (NaN 원인) ---\n');
fprintf('  평가 %d회 중 거부 %d회 (%.1f%%), tol=%.1e\n', ...
    numel(mism), sum(rej), 100*mean(rej), tol);

% (1) 재스케일 발동 여부와 거부의 상관 -> 원인 특정
if any(resc) || any(~resc)
    r1 = mean(rej(resc));   n1 = sum(resc);
    r0 = mean(rej(~resc));  n0 = sum(~resc);
    fprintf('  재스케일 발동(곡선 x>1): %d회 (%.1f%%), 이때 거부율 %.1f%%\n', ...
        n1, 100*mean(resc), 100*r1);
    fprintf('  재스케일 미발동        : %d회, 이때 거부율 %.1f%%\n', n0, 100*r0);
    if n1 > 0 && n0 > 0 && r1 > 0.5 && r0 < 0.1
        fprintf(['  => 거부는 거의 전부 재스케일에서 발생. 원인 확정:\n' ...
                 '     비단조 x 프로파일이 스플라인 overshoot 을 일으켜 재설정 불일치를 만든다.\n' ...
                 '     REQUIRE_MONOTONIC_X = true 로 두면 수율이 크게 오른다.\n']);
    elseif n1 > 0 && r0 > 0.2
        fprintf(['  => 재스케일과 무관하게도 거부가 발생. 단조성만으로는 부족하며\n' ...
                 '     GEOM_TOL(현재 %.1e) 상향을 함께 검토할 것.\n'], tol);
    end
    if any(maxlen > 1)
        fprintf('  overshoot 크기: median max_length=%.3f, max=%.3f\n', ...
            median(maxlen(maxlen>1)), max(maxlen));
    end
end

% (2) 불일치 크기 분포 -> tolerance 조정 판단
if any(rej)
    r = mism(rej);
    fprintf('  거부 시 불일치: median=%.2e, p90=%.2e, max=%.2e\n', ...
        median(r), prctile(r,90), max(r));
    if median(r) < 1e-2
        fprintf('     (불일치가 작음 -> GEOM_TOL 1e-3 상향도 유효)\n');
    else
        fprintf('     (불일치가 큼 -> 형상 왜곡. GEOM_TOL 올리지 말 것)\n');
    end
end
if any(~rej)
    fprintf('  통과 시 불일치: median=%.2e, max=%.2e\n', ...
        median(mism(~rej)), max(mism(~rej)));
end
end

%% ===== 무작위 valid 시드 생성 (pareto_front_freeform.m 과 동일) =====
%  REQUIRE_MONOTONIC_X 가 true 면 x2..x6 를 오름차순으로 정렬해 생성한다.
%  (기하 제약을 만족하면서도 LightTools 가 표현하지 못하던 비단조 프로파일을 배제)
function P = genValidPoints(K, lb, ub)
global REQUIRE_MONOTONIC_X
mono = ~isempty(REQUIRE_MONOTONIC_X) && REQUIRE_MONOTONIC_X;
dim = numel(lb);  P = zeros(K, dim);
for i = 1:K
    ok = false;
    while ~ok
        p = lb + rand(1, dim) .* (ub - lb);
        if mono
            p(1:5) = sort(p(1:5));          % x2..x6 오름차순
        end
        if isValidPoints(p), ok = true; P(i, :) = p; end
    end
end
end


%% ===== Objective (EQE_total 과 네 구간 EQE 동시 산출) =====
%  pareto_front_freeform.m 의 objFcn_both 를 그대로 복사 (수정 금지 구역).
function output = objFcn_both(point)
global ID_LT ID_swept ltml ltloc count ray_nums_current wave_n_current PATCH_XY
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

d_sub=1.295;  r_OLED=1;  Lensheight=0.01;
x_pattern = PATCH_XY;  y_pattern = PATCH_XY;   % 스크립트 상단 상수에서 받는다
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
rescaled = (max_length > 1);        % 재스케일 발동 여부 (거부와의 상관 진단용)
if rescaled
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
% [기하 검증] LightTools 에 설정한 제어점과 읽어온 값의 불일치 검사.
%   불일치가 크면 의도한 형상이 아니므로 거부한다(EQE_total=0 -> 상위에서 NaN).
global GEOM_TOL GEOM_MISMATCH_LOG
if isempty(GEOM_TOL), GEOM_TOL = 1e-4; end
mism = max(abs(xy(:) - xy_l(:)));

% COM 왕복 글리치 가능성 -> 1회 재설정 후 재확인
if mism > GEOM_TOL
    ltx.SetSweptProfilePoints(Curve,xy,7);
    ltx.DbSet(Curve,'StartSlopeMode',"Auto");
    ltx.DbSet(Curve,'EndSlopeMode',"Auto");
    for j=1:7
        xy_l(j,1) = ltml.LTDbGet(lt, Key, 'YAt', j);
        xy_l(j,2) = ltml.LTDbGet(lt, Key, 'ZAt', j);
    end
    mism = max(abs(xy(:) - xy_l(:)));
end

GEOM_MISMATCH_LOG(end+1,:) = [mism, max_length, double(rescaled)];
if mism > GEOM_TOL
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
% [주의] 파장 샘플 인덱스는 반드시 아래처럼 명시적으로 생성할 것.
%   나눗셈 기반 배열 크기 계산은 파장창/step 조합에 따라 비정수가 되어 오류를 냈다.
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

% ============================================================
%  convergence_check.m
%
%  목적: pareto_front_freeform.m 의 결론(구간별 선택성의 회전 추세, R = +0.60 … -0.57)이
%        실제 물리인지, 아니면 저정밀 시뮬의 인공물인지 판정한다.
%
%  [의심 요인 3가지]
%   (1) Ray 수 10,000 -> Monte-Carlo 노이즈. 선택성의 이항 노이즈는
%         sigma_S = sqrt(S(1-S)/N_esc)  (40-60 대역에서 상대 ~2%)
%       노이즈 자체는 상관을 만들지 않지만, 산포의 상당 부분을 차지할 수 있다.
%   (2) smooth(I_air_1_JH) -> MATLAB 기본 5점 이동평균. 90개 각도 빈에 적용되어
%       경계(0도 근처)에서 비대칭 왜곡을 만든다. 저효율 설계일수록 배광이 노이즈해
%       왜곡이 커진다면 EQE_total 과 '상관된 계통오차'가 생길 수 있다.
%   (3) 파장 590-600nm 만 사용 -> 캐비티의 I_sub(theta) 는 파장에 따라 이동하므로
%       광대역에서는 입력 분포가 달라질 수 있다.
%
%  [설계] 전체 재실행 없이, EQE_total 범위를 고루 덮는 설계 N_SEL 개만 골라
%         고정밀로 재평가하고 저정밀 결과와 비교한다.
%    STAGE A : 좁은 파장 유지 + ray 대폭 증가 + 반복  -> (1) 노이즈 판정
%    STAGE B : 광대역 파장                              -> (3) 스펙트럼 일반성 판정
%    두 스테이지 모두 smooth on/off 를 함께 계산        -> (2) smoothing 판정
%
%  [판정] 고정밀에서도 R 값이 유지되면 회전 추세는 실재. R 이 0 으로 붕괴하면 인공물.
%
%  [입력] pareto_front_result.mat   [출력] convergence_check_result.mat, .png, .txt
% ============================================================
clear;
global ID_swept ID_LT ltml ltloc count eval_count restart_interval ...
       ray_nums_current wave_n_current wave_range_current SMOOTH_ON

%% ===== 설정 =====
N_SEL      = 20;      % 재평가할 설계 개수 (EQE_total 구간별 층화 추출)
N_REP      = 3;       % 설계당 반복 (노이즈를 직접 측정)
RUN_STAGE_A = true;   % ray 수렴 (좁은 파장). 저비용·노이즈 판정에 결정적
RUN_STAGE_B = true;   % 광대역 파장. 비용 큼 -> 필요 없으면 false

RAY_HI     = 50000;  % STAGE A 고정밀 ray 수
RAY_B      = 20000;   % STAGE B ray 수 (광대역이라 파장 수가 많아 낮춤)
WAVE_NARROW = [590 600 10];   % [start end step]  (K=2)
WAVE_BROAD  = [450 750 10];   % [start end step]  (K=31)

restart_interval = 20;
count = 1;

%% ===== 데이터 로드 + 층화 추출 =====
D = load('pareto_front_result.mat');
LOG = D.EVAL_LOG;  nvar = numel(D.lb);
X    = LOG(:, 1:nvar);
EtLo = LOG(:, nvar+1);
BLo  = LOG(:, nvar+2:nvar+5);
ok   = isfinite(EtLo) & EtLo > 0.05 & all(isfinite(BLo),2);
X = X(ok,:);  EtLo = EtLo(ok);  BLo = BLo(ok,:);
fprintf('로그에서 유효 설계 %d개 확보\n', numel(EtLo));

% EQE_total 을 N_SEL 개 구간으로 나눠 각 구간에서 1개씩 -> 범위를 고루 덮음
edgesQ = linspace(min(EtLo), max(EtLo), N_SEL+1);
sel = [];
for i = 1:N_SEL
    idx = find(EtLo >= edgesQ(i) & EtLo < edgesQ(i+1));
    if isempty(idx), continue; end
    [~,j] = min(abs(EtLo(idx) - mean(edgesQ(i:i+1))));   % 구간 중앙에 가까운 것
    sel(end+1) = idx(j); %#ok<SAGROW>
end
sel = unique(sel);
Xs = X(sel,:);  EtLo_s = EtLo(sel);  BLo_s = BLo(sel,:);
nS = numel(sel);
fprintf('재평가 대상 %d개 (EQE_total %.3f ~ %.3f)\n\n', nS, min(EtLo_s), max(EtLo_s));

%% ===== LightTools 연결 =====
RenewLightTools();
lt = ltloc.GetLTAPI(ID_swept);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
eval_count = 0;

edgesA = [0 20; 20 40; 40 60; 60 80];
names  = {'0-20','20-40','40-60','60-80'};
S_pred = sind(edgesA(:,2)).^2 - sind(edgesA(:,1)).^2;

%% ===== 재평가 루틴 =====
%  smooth on/off 를 각각 계산해 (2) 를 판정
run_stage = @(rays, wrange, tag) evaluate_set(Xs, rays, wrange, N_REP, tag);

R = struct();
if RUN_STAGE_A
    fprintf('########## STAGE A: ray 수렴 (%d rays, %d-%d nm) ##########\n', ...
        RAY_HI, WAVE_NARROW(1), WAVE_NARROW(2));
    [R.A_sm, R.A_raw, R.A_std] = run_stage(RAY_HI, WAVE_NARROW, 'A');
end
if RUN_STAGE_B
    fprintf('\n########## STAGE B: 광대역 (%d rays, %d-%d nm) ##########\n', ...
        RAY_B, WAVE_BROAD(1), WAVE_BROAD(2));
    [R.B_sm, R.B_raw, R.B_std] = run_stage(RAY_B, WAVE_BROAD, 'B');
end

%% ===== 분석 + 판정 =====
fid = fopen('convergence_check.txt','w');
pr = @(varargin) (fprintf(varargin{:}) + fprintf(fid, varargin{:}));

pr('\n================ 수렴 검사 결과 ================\n');
pr('설계 %d개, 설계당 반복 %d회\n\n', nS, N_REP);

% (기준) 저정밀 로그에서의 상관
R_lo = corr_by_bin(EtLo_s, BLo_s);
pr('%-8s | %8s | %10s', 'band', '예측 S', '저정밀 R');
if RUN_STAGE_A, pr(' | %10s %10s', 'A:R(smooth)', 'A:R(raw)'); end
if RUN_STAGE_B, pr(' | %10s', 'B:R(광대역)'); end
pr('\n%s\n', repmat('-',1,78));

R_A_sm = nan(4,1); R_A_raw = nan(4,1); R_B_sm = nan(4,1);
if RUN_STAGE_A
    R_A_sm  = corr_by_bin(R.A_sm(:,1),  R.A_sm(:,2:5));
    R_A_raw = corr_by_bin(R.A_raw(:,1), R.A_raw(:,2:5));
end
if RUN_STAGE_B
    R_B_sm  = corr_by_bin(R.B_sm(:,1),  R.B_sm(:,2:5));
end
for b = 1:4
    pr('%-8s | %8.3f | %+10.2f', names{b}, S_pred(b), R_lo(b));
    if RUN_STAGE_A, pr(' | %+10.2f %+10.2f', R_A_sm(b), R_A_raw(b)); end
    if RUN_STAGE_B, pr(' | %+10.2f', R_B_sm(b)); end
    pr('\n');
end

% 노이즈 직접 측정 (반복 간 표준편차)
if RUN_STAGE_A
    pr('\n[노이즈 직접 측정] STAGE A 반복 간 선택성 표준편차 (상대):\n');
    for b = 1:4
        s_mean = mean(R.A_sm(:,1+b)./R.A_sm(:,1), 'omitnan');
        pr('  %-8s : %.2f%%  (설계 간 산포와 비교해 판단)\n', ...
            names{b}, 100*mean(R.A_std(:,b),'omitnan')/max(s_mean,eps));
    end
end

% --- 자동 판정 ---
pr('\n[판정]\n');
if RUN_STAGE_A
    keep = abs(R_A_sm([1 4])) > 0.3;      % 양 끝 구간에서 추세가 살아있는가
    if all(keep)
        pr('  (1) ray 증가 후에도 R 유지 => 회전 추세는 실재. 원고에 사용 가능.\n');
    elseif ~any(keep)
        pr('  (1) ray 증가 시 R 붕괴 => 저정밀 노이즈 인공물. "불변량" 주장이 오히려 성립.\n');
    else
        pr('  (1) R 이 부분적으로만 유지 => 추세 약함. 원고에서 강한 주장 금지.\n');
    end
    dsm = max(abs(R_A_sm - R_A_raw));
    if dsm > 0.25
        pr('  (2) smooth on/off 로 R 이 %.2f 변동 => smoothing 인공물 유의. raw 기준으로 보고할 것.\n', dsm);
    else
        pr('  (2) smooth 영향 작음 (max dR = %.2f).\n', dsm);
    end
end
if RUN_STAGE_B
    db = max(abs(R_B_sm - R_A_sm));
    if isfinite(db) && db > 0.3
        pr('  (3) 광대역에서 R 이 %.2f 변동 => 좁은 파장 결과를 일반화 금지. 광대역 기준으로 다시 보고할 것.\n', db);
    else
        pr('  (3) 광대역에서도 유사 => 스펙트럼 일반성 확보.\n');
    end
end
fclose(fid);
fprintf('\nsaved -> convergence_check.txt\n');

save('convergence_check_result.mat','R','Xs','EtLo_s','BLo_s','R_lo', ...
     'R_A_sm','R_A_raw','R_B_sm','S_pred','edgesA','names');

%% ===== 그림 =====
figure('Name','Convergence check','Color','w','Position',[80 80 1150 430]);
subplot(1,2,1);
bh = bar([R_lo(:), R_A_sm(:), R_B_sm(:)]);
set(gca,'XTickLabel',names); ylabel('상관 R (선택성 vs EQE_{total})');
legend({'저정밀(10k, 좁은 \lambda)','A: 고 ray','B: 광대역'},'Location','best','FontSize',8);
yline(0.3,'--'); yline(-0.3,'--'); grid on;
title('(a) R 값의 정밀도 의존성');

subplot(1,2,2);
if RUN_STAGE_A
    cols = lines(4);
    for b = 1:4
        plot(EtLo_s, BLo_s(:,b)./EtLo_s, 'o', 'Color', [cols(b,:) 0.4], ...
            'MarkerSize',5); hold on;
        plot(R.A_sm(:,1), R.A_sm(:,1+b)./R.A_sm(:,1), 's', 'Color', cols(b,:), ...
            'MarkerFaceColor', cols(b,:), 'MarkerSize',6);
        yline(S_pred(b),'--','Color',cols(b,:));
    end
    xlabel('EQE_{total}'); ylabel('선택성 S');
    title('(b) 저정밀(o) vs 고정밀(\Box)'); grid on;
end
saveas(gcf,'convergence_check.png');
fprintf('saved -> convergence_check.png\n');


%% ===== 설계 집합 재평가 =====
function [M_sm, M_raw, S_std] = evaluate_set(Xs, rays, wrange, nrep, tag)
global ray_nums_current wave_n_current wave_range_current SMOOTH_ON
nS = size(Xs,1);
M_sm  = nan(nS,5);   % [EQE_total, b0_20, b20_40, b40_60, b60_80]  (smooth 적용)
M_raw = nan(nS,5);   % (smooth 미적용)
S_std = nan(nS,4);   % 반복 간 선택성 표준편차
ray_nums_current = rays;
wave_range_current = wrange;
wave_n_current = wrange(3);
for i = 1:nS
    acc_sm = nan(nrep,5);  acc_raw = nan(nrep,5);
    for r = 1:nrep
        SMOOTH_ON = true;   o1 = safe_eval(Xs(i,:));
        SMOOTH_ON = false;  o2 = safe_eval(Xs(i,:));
        acc_sm(r,:)  = o1;  acc_raw(r,:) = o2;
    end
    M_sm(i,:)  = mean(acc_sm, 1,'omitnan');
    M_raw(i,:) = mean(acc_raw,1,'omitnan');
    sel_rep = acc_sm(:,2:5) ./ acc_sm(:,1);
    S_std(i,:) = std(sel_rep, 0, 1, 'omitnan');
    fprintf('  [%s] %2d/%d : EQE_total=%.4f  S(40-60)=%.4f\n', ...
        tag, i, nS, M_sm(i,1), M_sm(i,4)/M_sm(i,1));
end
end

function v = safe_eval(pt)
global ID_swept ltml ltloc eval_count restart_interval
eval_count = eval_count + 1;
if mod(eval_count, restart_interval) == 0
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    pause(2);
end
v = nan(1,5);
try
    o = objFcn_conv(pt);
    if o.EQE_total > 0
        v = [o.EQE_total, o.EQE_0_20, o.EQE_20_40, o.EQE_40_60, o.EQE_60_80];
    end
catch err
    fprintf('    [Error] %s\n', err.message);
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
end
end

function R = corr_by_bin(Et, B)
R = nan(4,1);
for b = 1:4
    s = B(:,b)./Et;
    g = isfinite(s) & isfinite(Et);
    if sum(g) > 3
        c = corrcoef(Et(g), s(g));  R(b) = c(1,2);
    end
end
end


%% ===== Objective (파장창·ray·smooth 를 외부에서 지정) =====
function output = objFcn_conv(point)
global ID_LT ID_swept ltml ltloc count ray_nums_current wave_n_current ...
       wave_range_current SMOOTH_ON
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

d_sub=1.295;  r_OLED=1;  x_pattern=25;  y_pattern=25;  Lensheight=0.01;
wavelength_start = wave_range_current(1);
wavelength_end   = wave_range_current(2);
n                = wave_range_current(3);
ray_nums         = ray_nums_current;

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

x2=point(1); x3=point(2); x4=point(3); x5=point(4); x6=point(5);
y2=point(6); y3=point(7); y4=point(8); y5=point(9); y6=point(10);
dETL=point(11); dHTL=point(12); stretchZ=point(13);

xy = zeros(7,2);
xy(1,:) = [0, 1];  xy(7,:) = [1, 0];
xy(2,:) = [x2,y2]; xy(3,:) = [x3,y3]; xy(4,:) = [x4,y4];
xy(5,:) = [x5,y5]; xy(6,:) = [x6,y6];

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
for a=1:101, x_values(a)=ltml.LTDbGet(lt,Key,'YFacetsAt',a); end
max_length = max(x_values);
if max_length > 1, xy = xy / max_length; end
ltx.SetSweptProfilePoints(Curve,xy,7);
ltx.DbSet(Curve,'StartSlopeMode',"Auto");
ltx.DbSet(Curve,'EndSlopeMode',"Auto");
xy_l = zeros(7,2);
for j=1:7
    xy_l(j,1) = ltml.LTDbGet(lt, Key, 'YAt', j);
    xy_l(j,2) = ltml.LTDbGet(lt, Key, 'ZAt', j);
end
if max(abs(xy(:) - xy_l(:))) > 1e-4
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
ltml.LTCmd(lt, 'Undo');  ltml.LTCmd(lt, 'Undo');
totalpathmod = [pathname index '.1.ent"'];
List = ltml.LTDbList(lt2, 'LENS_MANAGER[1]', 'LIBRARY_ELEMENT_UNIT_CELL');
Key = ltml.LTListByName(lt2, List, 'LibraryElement');
ltml.LTDbSet(lt2, Key, 'Filename', totalpathmod);
List = ltml.LTDbList(lt2, 'LENS_MANAGER[1]', 'TEXTURE_PARAMETER');
Key = ltml.LTListByName(lt2, List, 'StretchZ');
ltml.LTDbSet(lt2, Key, 'Value', stretchZ);

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
CPS_result=CPS_for_Isub(no_bar,ne_bar,thickness,emission_spectrum,eta_rad, ...
    horizontal_dipole_ratio,bottom_air_refractive_index,4,12.5,499,3,wavelength);
EQE_sub_CPS=CPS_result.EQE_sub;

TMF_p=TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_s=TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
Reflectance=(abs(TMF_p.r_p).^2 + abs(TMF_s.r_s).^2)/2;

lt = ltloc.GetLTAPI(ID_LT);
fileID = fopen(sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count), 'w');
fprintf(fileID,'%s\n%s%d\n%s\n%s\n%s\n%s\n ','DFAT Version 1.0','DATANAME: R_Bottom_',count,'ABSORBING: YES','INDEX: 1.51','DATAITEMS: TAVG RAVG');
for i=wavelength_start:wavelength_end
    fprintf(fileID,'%s  %d\n','wv',i);
    for j=0:89
        fprintf(fileID,'%s  %d  %d  %.3f\n','AOI',j,0,Reflectance(i-wavelength_start+1,j+1));
    end
end
fclose(fileID);
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings]" LoadFileName="' sprintf('C:\\Users\\jhkim\\Desktop\\Green_CE_Calculation\\TRA_temp\\R_Al_%d.coa', count) '"']);
List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');
Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');
Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName',sprintf('R_Bottom_%d', count));

I_white=0.5*(CPS_result.I_sub_s+CPS_result.I_sub_p);
P_white=I_white.*repmat(sin089,wavelength_num,1);
weight_factor=sum(P_white,2);
wv_list = 1:n:wavelength_num;   K = numel(wv_list);
I_air_1_2 = zeros(90,K);  Power_output = zeros(1,wavelength_num);
for kk = 1:K
    wv = wv_list(kk);
    fileID = fopen('C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','w');
    fprintf(fileID,'%s  %d  %d  %d  %d  %d  %d','SPHEREMESH:',1,90,0,0,360,90);
    writematrix(flip(I_white(wv,:).'),'C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt','Delimiter','tab','WriteMode','append');
    fclose(fileID);
    SRList=ltml.LTDbList(lt,'Lens_manager[1]','DISK_SOURCE');
    SRKey=ltml.LTListAtPos(lt,SRList,1);
    ltml.LTDbSet(lt,SRKey,'Radiant_Power', weight_factor(wv));
    SRList=ltml.LTDbList(lt,'Lens_manager[1]','Spectral_region');
    SRKey=ltml.LTListAtPos(lt,SRList,2);
    ltml.LTDbSet(lt,SRKey,'Spectral_Definition','Monochromatic');
    ltml.LTDbSet(lt,SRKey,'Single_Wavelength', wv+wavelength_start-1);
    List=ltml.LTDbList(lt,'lens_manager[1]','DIRECTION_GRID_APODIZER');
    Key=ltml.LTListAtPos(lt,List,1);
    ltml.LTDbSet(lt,Key,'LoadFileName','C:\Users\jhkim\Desktop\Green_CE_Calculation\Angular_temp\AI_temp.txt');
    ltml.LTBegin(lt);
    ltml.LTCmd(lt,'\V3D BeginAllSimulations');
    ltml.LTEnd(lt);
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,1);
    Power_output(wv)=ltml.LTDbGet(lt,Key,'TotalPower');
    List=ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key=ltml.LTListAtPos(lt,List,3);
    I_raw = zeros(90,1);
    for j=1:90
        I_raw(91-j) = ltml.LTDbGet(lt,Key,'CellValue_UI',1,91-j);
    end
    % [검증 대상] smooth 적용 여부를 외부에서 토글
    if isempty(SMOOTH_ON) || SMOOTH_ON
        I_air_1_2(:,kk) = smooth(I_raw);
    else
        I_air_1_2(:,kk) = I_raw;
    end
end

weight_factor_2=zeros(K,1); Power_output_2=zeros(K,1); EQE_sub_matrix_2=zeros(K,1);
for k = 1:K
    idx = wv_list(k);
    weight_factor_2(k)=weight_factor(idx);
    Power_output_2(k)=Power_output(idx);
    EQE_sub_matrix_2(k)=CPS_result.EQE_sub_matrix(idx);
end
EQE_wv_matrix = Power_output_2 ./ weight_factor_2;
EQE_sub_matrix_2 = EQE_sub_matrix_2 / sum(EQE_sub_matrix_2) * EQE_sub_CPS;
EQE_total = sum(EQE_wv_matrix .* EQE_sub_matrix_2);

E=zeros(1,4); sin_col=sin089(:);
for k = 1:K
    contrib = EQE_wv_matrix(k)*EQE_sub_matrix_2(k);
    W = I_air_1_2(:,k).*sin_col;  Wt = sum(W);
    E(1)=E(1)+contrib*sum(W(1:20))/Wt;
    E(2)=E(2)+contrib*sum(W(21:40))/Wt;
    E(3)=E(3)+contrib*sum(W(41:60))/Wt;
    E(4)=E(4)+contrib*sum(W(61:80))/Wt;
end
output = struct('EQE_0_20',E(1),'EQE_20_40',E(2),'EQE_40_60',E(3), ...
                'EQE_60_80',E(4),'EQE_total',EQE_total);

List=ltml.LTDbList(lt,'lens_manager[1]','PROPERTY');
Key=ltml.LTListByName(lt,List,'R_Al');
List=ltml.LTDbList(lt,Key,'USER_COATING_AMPLITUDE_ZONE');
Key=ltml.LTListNext(lt,List);
ltml.LTDbSet(lt,Key,'SelectedCoatingName','R_temp');
ltml.LTCmd(lt,['\O"LENS_MANAGER[1].USER_COATINGS[User Coatings].COATING[' sprintf('R_Bottom_%d', count) ']" Delete= \Q']);
fclose('all');
end


function RenewLightTools()
global ID_LT ID_swept ltml ltloc lt
lt_exe_path = 'C:\Program Files\Optical Research Associates\LightTools 2023.03\lt.exe';
model_file_path_swept = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\SweptEntity.2.lts';
model_file_path_LT = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\Lens_size_effect_for_PSO_bump_modified_v1.1.lts';
fprintf('--- Restarting LightTools ---\n');
target_user = 'jhkim';
[~,~] = system(sprintf('taskkill /F /FI "USERNAME eq %s" /IM lt.exe', target_user));
pause(2);
system(sprintf('"%s" "%s" &', lt_exe_path, model_file_path_swept));
try
    ltml = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
catch
    error('LightTools 재시작 실패.');
end
find_cmd = sprintf('tasklist /fi "imagename eq lt.exe" /fi "username eq %s" /fo csv /nh', target_user);
[st,out] = system(find_cmd);
tok = regexp(out,'"(\d+)"','tokens');
if st==0 && ~isempty(tok), ID_swept = str2double(tok{1}{1});
else, error('PID 추출 실패'); end
system(sprintf('"%s" "%s" &', lt_exe_path, model_file_path_LT));
[st,out] = system(find_cmd);
tok = regexp(out,'"(\d+)"','tokens');
if st==0 && numel(tok)>=3, ID_LT = str2double(tok{3}{1});
else, error('PID 추출 실패'); end
pause(5);
end

% ============================================================
%  stress_random_mla.m
%
%  목적: 주기적 freeform MLA 에서 확인된 "practical saturation"
%        (효율과 모든 각도구간이 함께 상승; 선택성이 Lambertian 분할
%         [0.117 0.296 0.337 0.220] 근처; drift 상관 +0.6/+0.7/+0.05/-0.7)
%        이 비주기·무작위 조립 MLA (예: 무작위 microsphere assembly) 에서도
%        성립하는지 스트레스 테스트.
%
%  [방법] 15x15 mm 를 무작위 렌즐릿 수백만 개로 직접 채울 수 없으므로
%   pseudo-random supercell: 기존 텍스처 unit-cell 크기의 슈퍼셀 하나에
%   렌즐릿 위치/크기/프로파일을 무작위로 채우고(generate_random_supercell_ent.m),
%   기존 타일링(x_pattern x y_pattern = 15 x 15)이 그 슈퍼셀을 반복한다.
%   상관길이 << 슈퍼셀 이므로 통계적으로 무작위 어레이와 동등.
%
%  [설계공간] 개별 렌즈 형상이 아니라 무작위 조립의 '통계' = 6차원 하이퍼파라미터:
%     [fill, rJitter, posJitter, aspect, aspectJitter, profileMix]
%   따라서 프로토콜이 다른 가족과 약간 다르다:
%     Phase A: N_RANDOM=100 realization (매번 새 seed + 무작위 하이퍼벡터, 탐색 정밀도)
%     Phase B: surrogateopt 60 + patternsearch 15, 하이퍼벡터만 최적화
%              (평가마다 seed 고정 -> 목적함수 결정론적)
%     Phase C: 승자를 서로 다른 seed 3개로 고정밀 재평가 (seed-강건성; mean±std)
%
%  [수집물] EVAL_LOG (stress_random_result.mat):
%     [ hyp(1:6) | NaN(1:7) | EQE_total | b0_20 | b20_40 | b40_60 | b60_80 | phase | w ]
%     * 하이퍼파라미터는 6차원이지만 다른 가족의 13-var 로그와 열을 맞추기 위해
%       NaN 7개로 13열까지 패딩한다 (열 7..13 = NaN). phase=7, w=-1 (이 가족 표식).
%     * 각 행의 (seed, 서브페이즈 1=random/2=opt/3=final) 는 EVAL_META 에 병행 기록.
%  [그림] stress_random_check.png — (a) b40_60 vs EQE_total, (b) 승자 선택성 vs
%     Lambertian, (c) 구간별 drift 상관 R vs freeform 기준선. + 텍스트 판정.
%
%  기반: pareto_front_freeform.m (LightTools 연동/CPS/파장루프/코팅 전부 동일.
%        렌즈 주입부만 swept-curve 대신 슈퍼셀 .ent 로 교체.)
% ============================================================
clear;
%% For LightTools Connection
global ID_swept ID_LT ltml ltloc count eval_count restart_interval ...
       ray_nums_current wave_n_current EVAL_LOG EVAL_META EVAL_PHASE EVAL_W

RenewLightTools();
try
    ltml.LTCmd(ltml.GetLTAPI(ID_LT), 'Message "Check Connection"');
catch
    ltml = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
end
count = 1;
restart_interval = 20;
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

%% ===== Multi-fidelity (pareto_front_freeform.m 과 동일) =====
WAVE_N_SEARCH = 10;
WAVE_N_FINAL  = 2;
RAY_SEARCH    = 10000;
RAY_FINAL     = 50000;
N_FINAL_REP   = 3;      % Phase C: 서로 다른 seed 3개 (seed-강건성 검정)

%% ===== 프로토콜 예산 =====
N_RANDOM     = 100;     % Phase A: 무작위 realization 수
OPT_EVALS    = 60;      % Phase B: surrogateopt 예산
POLISH_EVALS = 15;      % Phase B: patternsearch 예산
MIN_SURR     = 20;

SEED_RANDOM_BASE = 1000;            % Phase A seed = 1000+i
SEED_OPT         = 777;             % Phase B: 모든 평가 동일 seed (결정론적 목적함수)
SEEDS_FINAL      = [2001 2002 2003];% Phase C: seed-강건성

%% ===== 하이퍼파라미터 벡터 (6차원) =====
%  generate_random_supercell_ent.m 의 params 필드와 1:1 대응.
hypNames = {'fill','rJitter','posJitter','aspect','aspectJitter','profileMix'};
lb = [0.35, 0.00, 0.0, 0.3, 0.00, 0.0];
ub = [0.90, 0.30, 1.0, 1.5, 0.40, 1.0];
nhyp = numel(lb);
NPAD = 13;                          % 다른 가족 로그(13-var)와 열 폭 맞춤 (NaN 패딩)

EVAL_LOG  = [];   % [hyp(1:6) NaN(1:7) | EQE_total b0_20 b20_40 b40_60 b60_80 | phase w]
EVAL_META = [];   % [seed, subphase]  (1=random, 2=opt, 3=final) — EVAL_LOG 와 행 대응
EVAL_PHASE = 7;   % 이 가족(무작위 조립 MLA)의 phase 표식
EVAL_W     = -1;

% freeform(주기) 가족의 drift 상관 기준선 (원고 값): R(EQE_total, S_j)
R_FREEFORM = [0.6, 0.7, 0.05, -0.7];
edges  = [0 20; 20 40; 40 60; 60 80];
names  = {'0-20 deg','20-40 deg','40-60 deg','60-80 deg'};
S_pred = sind(edges(:,2)).^2 - sind(edges(:,1)).^2;   % [0.117 0.296 0.337 0.220]

psOpts = optimoptions('patternsearch', ...
    'MaxFunctionEvaluations', POLISH_EVALS, ...
    'InitialMeshSize', 0.1, 'MeshTolerance', 1e-3, ...
    'Cache', 'on', 'Display', 'off');

%% =====================================================================
%  PHASE A — 무작위 realization (새 seed + 무작위 하이퍼벡터)
%% =====================================================================
fprintf('\n########## PHASE A: random realizations (N=%d) ##########\n', N_RANDOM);
ray_nums_current = RAY_SEARCH;  wave_n_current = WAVE_N_SEARCH;
eval_count = 0;
for i = 1:N_RANDOM
    hyp  = lb + rand(1, nhyp) .* (ub - lb);
    seed = SEED_RANDOM_BASE + i;
    [et, bins] = simulate_supercell(hyp, seed, 1);
    if mod(i,10)==0
        fprintf('  random %3d/%d : EQE_total=%.4f  b40_60=%.4f\n', i, N_RANDOM, et, bins(3));
        save('stress_random_result.mat','EVAL_LOG','EVAL_META','hypNames','lb','ub');  % 중간저장
    end
end
save('stress_random_result.mat','EVAL_LOG','EVAL_META','hypNames','lb','ub');

%% =====================================================================
%  PHASE B — 하이퍼파라미터 최적화 (seed 고정 -> 결정론적)
%  목적: EQE_total 최대화. saturation 이면 최적화가 모든 구간을 함께 끌어올린다.
%% =====================================================================
fprintf('\n########## PHASE B: surrogateopt over hyperparameters ##########\n');
RenewLightTools();
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
eval_count = 0;
ray_nums_current = RAY_SEARCH;  wave_n_current = WAVE_N_SEARCH;

sopts = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', OPT_EVALS, ...
    'MinSurrogatePoints',     MIN_SURR, ...
    'UseParallel', false, 'PlotFcn', [], 'Display', 'iter');
[xS, ~] = surrogateopt(@(x) opt_obj(x, SEED_OPT), lb, ub, sopts);

if ~isempty(xS)
    x0 = xS(:).';
    try
        xP = patternsearch(@(x) opt_obj(x, SEED_OPT), x0, [],[],[],[], lb, ub, [], psOpts);
        xP = xP(:).';
    catch
        xP = x0;
    end
else
    error('surrogateopt 이 해를 반환하지 않음.');
end

% 후보 중 저정밀 기준 우위 채택
cands = {x0};
if ~isequal(xP, x0), cands{end+1} = xP; end
bestScore = -inf;  bestHyp = [];
for c = 1:numel(cands)
    [et, ~] = simulate_supercell(cands{c}, SEED_OPT, 2);
    if isfinite(et) && et > bestScore, bestScore = et; bestHyp = cands{c}; end
end
if isempty(bestHyp), bestHyp = x0; end
fprintf('  승자 하이퍼벡터: [%s]\n', num2str(bestHyp, '%.3f '));
save('stress_random_result.mat','EVAL_LOG','EVAL_META','hypNames','lb','ub','bestHyp');

%% =====================================================================
%  PHASE C — 고정밀 재평가, 서로 다른 seed 3개 (seed-강건성)
%% =====================================================================
fprintf('\n########## PHASE C: high-fidelity, %d seeds ##########\n', numel(SEEDS_FINAL));
ray_nums_current = RAY_FINAL;  wave_n_current = WAVE_N_FINAL;
et_f   = nan(1, N_FINAL_REP);
bins_f = nan(N_FINAL_REP, 4);
for r = 1:N_FINAL_REP
    [et_f(r), bins_f(r,:)] = simulate_supercell(bestHyp, SEEDS_FINAL(r), 3);
    fprintf('  seed %d : EQE_total=%.5f  bands=[%.4f %.4f %.4f %.4f]\n', ...
        SEEDS_FINAL(r), et_f(r), bins_f(r,:));
end
et_mean = mean(et_f,'omitnan');   et_std = std(et_f,'omitnan');
fprintf('  >>> 승자 (seed 3개): EQE_total = %.5f ± %.5f  (CV %.2f%%)\n', ...
    et_mean, et_std, 100*et_std/max(et_mean,eps));
save('stress_random_result.mat','EVAL_LOG','EVAL_META','hypNames','lb','ub', ...
     'bestHyp','et_f','bins_f','SEEDS_FINAL','R_FREEFORM','S_pred');

%% =====================================================================
%  분석 + 그림 stress_random_check.png (3패널) + 판정
%% =====================================================================
nv = NPAD;
Et   = EVAL_LOG(:, nv+1);
Bins = EVAL_LOG(:, nv+2 : nv+5);
sub  = EVAL_META(:,2);
ok = isfinite(Et) & Et > 0.05 & all(isfinite(Bins),2);   % 저효율은 노이즈 지배

% 구간별 drift 상관 R(EQE_total, S_j) — random phase 만 (편향 없는 표본)
okR = ok & sub==1;
R_meas = nan(1,4);  S_meas = nan(1,4);  S_std_m = nan(1,4);
for b = 1:4
    s = Bins(okR,b) ./ Et(okR);
    S_meas(b) = mean(s,'omitnan');  S_std_m(b) = std(s,'omitnan');
    c = corrcoef(Et(okR), s);  R_meas(b) = c(1,2);
end
% 승자의 선택성 (고정밀, seed 평균)
S_best = mean(bins_f ./ et_f(:), 1, 'omitnan');

figure('Name','stress random MLA — saturation check','Color','w','Position',[80 80 1400 420]);

subplot(1,3,1);   % (a) b40_60 vs EQE_total + 선형적합
xa = Et(ok);  ya = Bins(ok,3);
scatter(xa, ya, 16, [.5 .5 .5], 'filled', 'MarkerFaceAlpha', 0.5); hold on;
pf = polyfit(xa, ya, 1);
xx = linspace(min(xa), max(xa), 50);
plot(xx, polyval(pf, xx), 'r-', 'LineWidth', 1.8);
plot(et_mean, mean(bins_f(:,3)), 'p', 'MarkerSize', 14, ...
    'MarkerFaceColor',[.9 .6 .1], 'MarkerEdgeColor','k');
xlabel('EQE_{total}'); ylabel('EQE_{40-60}'); grid on;
title(sprintf('(a) b_{40-60} vs EQE_{total}  (slope=%.3f)', pf(1)));
legend({'realizations','linear fit','winner (hi-fi)'},'Location','best','FontSize',8);

subplot(1,3,2);   % (b) 승자 선택성 vs Lambertian
bar([S_pred(:), S_best(:)]); grid on;
set(gca,'XTickLabel',names,'FontSize',8);
ylabel('selectivity  S_j = EQE_{band}/EQE_{total}');
legend({'Lambertian  sin^2\theta_2 - sin^2\theta_1','best random MLA'}, ...
    'Location','northwest','FontSize',8);
title('(b) winner selectivity vs Lambertian');

subplot(1,3,3);   % (c) drift 상관 vs freeform 기준선
bar([R_FREEFORM(:), R_meas(:)]); grid on;
set(gca,'XTickLabel',names,'FontSize',8);  ylim([-1 1]);  yline(0,'k-');
ylabel('R( EQE_{total},  S_j )');
legend({'periodic freeform (baseline)','random MLA (this test)'}, ...
    'Location','southwest','FontSize',8);
title('(c) per-band drift correlation');

%% --- 텍스트 판정 ---
% saturation 서명 3가지 + seed 강건성:
%  (1) 선택성이 Lambertian 분할의 15% 이내
%  (2) 모든 구간이 효율과 함께 상승: corr(EQE_total, EQE_band_j) > 0, all j
%  (3) drift 상관 부호 패턴이 freeform 기준선과 일치 (구간 3은 |R|<0.3 = ~0)
%  (4) seed 간 CV < 5% (슈퍼셀 통계가 realization 에 강건)
rise = nan(1,4);
for b = 1:4
    c = corrcoef(Et(okR), Bins(okR,b));  rise(b) = c(1,2);
end
dev_ok   = all(abs(S_meas(:) - S_pred(:)) ./ S_pred(:) < 0.15);
rise_ok  = all(rise > 0);
sign_ok  = (sign(R_meas(1))==sign(R_FREEFORM(1))) && ...
           (sign(R_meas(2))==sign(R_FREEFORM(2))) && ...
           (abs(R_meas(3)) < 0.3) && ...
           (sign(R_meas(4))==sign(R_FREEFORM(4)));
seed_ok  = (et_std / max(et_mean,eps)) < 0.05;

verdict = sprintf(['[판정 — 무작위 조립 MLA 의 saturation]\n' ...
    ' (1) 선택성 Lambertian 15%% 이내 : %s   S=[%.3f %.3f %.3f %.3f]\n' ...
    ' (2) 전 구간 동반 상승 (R>0)    : %s   R_rise=[%+.2f %+.2f %+.2f %+.2f]\n' ...
    ' (3) drift 상관 패턴 재현       : %s   R=[%+.2f %+.2f %+.2f %+.2f] (기준 [%+.2f %+.2f %+.2f %+.2f])\n' ...
    ' (4) seed 강건성 (CV<5%%)       : %s   EQE=%.4f±%.4f\n' ...
    ' => %s\n'], ...
    tern(dev_ok), S_meas, tern(rise_ok), rise, tern(sign_ok), R_meas, R_FREEFORM, ...
    tern(seed_ok), et_mean, et_std, ...
    tern2(dev_ok && rise_ok && sign_ok && seed_ok, ...
      'practical saturation 은 비주기·무작위 조립 MLA 에서도 재현된다.', ...
      'saturation 서명 일부 미재현 — 위 실패 항목을 원고에서 별도 논의할 것.'));
fprintf('\n%s\n', verdict);
annotation('textbox',[0.01 0.0 0.98 0.08],'String',verdict,'FontSize',7, ...
    'Interpreter','none','EdgeColor','none');

saveas(gcf,'stress_random_check.png');
save('stress_random_result.mat','EVAL_LOG','EVAL_META','hypNames','lb','ub', ...
     'bestHyp','et_f','bins_f','SEEDS_FINAL','R_FREEFORM','S_pred', ...
     'S_meas','S_std_m','R_meas','rise','S_best','verdict');
fprintf('saved -> stress_random_result.mat / stress_random_check.png\n');
fprintf('\n########## 완료 ##########\n');


%% ===== 최적화 목적 (seed 고정 -> 결정론적) =====
function f = opt_obj(x, seed)
[et, ~] = simulate_supercell(x(:).', seed, 2);
if ~isfinite(et), f = 0; else, f = -et; end   % EQE_total 최대화
end

%% ===== 평가 래퍼 (pareto_front_freeform.m 의 simulate_both 미러) =====
%  EVAL_LOG 행: [hyp(1:6) NaN(1:7) | EQE_total b0_20 b20_40 b40_60 b60_80 | 7 | -1]
function [eqe_total, bins] = simulate_supercell(hyp, seed, subphase)
global ID_LT ltml ltloc eval_count restart_interval EVAL_LOG EVAL_META EVAL_PHASE EVAL_W
eval_count = eval_count + 1;
if mod(eval_count, restart_interval) == 0
    fprintf('\n[Refresh] 시뮬 %d회. LightTools 재시작...\n', eval_count);
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_LT);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    pause(2);
end
bins = [NaN NaN NaN NaN];
try
    out = objFcn_supercell(hyp, seed);
    eqe_total = out.EQE_total;
    bins = [out.EQE_0_20, out.EQE_20_40, out.EQE_40_60, out.EQE_60_80];
    if eqe_total == 0
        eqe_total = NaN; bins = [NaN NaN NaN NaN];
    end
catch err
    fprintf('\n[Error] eval %d LightTools 충돌: %s\n', eval_count, err.message);
    eqe_total = NaN;
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_LT);  ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
end
pad = nan(1, 13 - numel(hyp));                   % 13열 패딩 (문서화: 열 7..13 = NaN)
EVAL_LOG(end+1,:)  = [hyp(:).', pad, eqe_total, bins, EVAL_PHASE, EVAL_W];
EVAL_META(end+1,:) = [seed, subphase];
end

%% ===== 판정 문구 헬퍼 =====
function s = tern(c)
if c, s = '통과'; else, s = '실패'; end
end
function s = tern2(c, a, b)
if c, s = a; else, s = b; end
end

%% ===== Objective — objFcn_both 미러, 렌즈 주입부만 슈퍼셀 .ent 로 교체 =====
%  [변경점 요약]
%   - swept-curve 제어점 설정/왕복검증/SaveLibrary 블록 삭제 (곡선이 없으므로
%     GEOM_TOL 검사도 해당 없음). 대신 MATLAB 이 직접 슈퍼셀 .ent 를 쓰고
%     LIBRARY_ELEMENT_UNIT_CELL 'LibraryElement' 의 Filename 만 바꾼다
%     (objFcn_both 이 swept_XXX.1.ent 를 물리던 것과 동일한, 검증된 COM 패턴).
%   - 높이가 .ent 에 이미 구워져 있으므로 텍스처 StretchZ = 1 고정.
%   - 그 외 (기판/존/소스/CPS/코팅/파장루프/각도적분) 전부 objFcn_both 와 동일.
function output = objFcn_supercell(hyp, seed)
global ID_LT ltml ltloc count ray_nums_current wave_n_current
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

%% --- 슈퍼셀 .ent 생성 + unit-cell 교체 (objFcn_both 의 swept 블록 대체) ---
BASE = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\';
scDir = [BASE 'supercell_ents\'];
if ~exist(scDir,'dir'), mkdir(scDir); end

params = struct();
params.fill         = hyp(1);
params.rJitter      = hyp(2);
params.posJitter    = hyp(3);
params.aspect       = hyp(4);
params.aspectJitter = hyp(5);
params.profileMix   = hyp(6);
params.templatePath = [BASE 'freeform_template_v2.ent'];   % 원점정렬 검증 템플릿
% 파일명: objFcn_both 의 swept_<tag> 관례를 따르되 결정론성은 (seed,hyp) 가 보장
rng('shuffle');
charSet = ['a':'z' 'A':'Z' '0':'9'];
index = charSet(randi(length(charSet), 1, 10));
params.outPath = [scDir 'supercell_' index '.1.ent'];

generate_random_supercell_ent(seed, params);   % 내부에서 rng(seed) -> 결정론적

lt2 = ltloc.GetLTAPI(ID_LT);
totalpathmod = ['"' params.outPath '"'];       % objFcn_both 와 동일하게 따옴표 포함
%% !!! VERIFY IN LIGHTTOOLS: supercell entity loading — check on first run !!!
%  아래 3 call 은 objFcn_both 이 swept_XXX.1.ent 를 물리던 검증된 패턴 그대로다.
%  다만 (i) SaveLibrary 산출물이 아닌 MATLAB 직접 작성 freeform .ent (SurfacePairLens/
%  FreeformEntity) 가 같은 unit-cell 에 물리는지, (ii) RepairEntities 없이 유효
%  solid 로 로드되는지, (iii) StretchZ=1 이 이 엔티티에 항등으로 작용하는지를
%  첫 실행에서 GUI 로 육안 확인할 것. (generate_freeform_ent.m 헤더가 같은 방식
%  연결을 전제하며, ent_tests/ 의 무작위 .ent 들이 수동 로드 검증된 전례.)
List = ltml.LTDbList(lt2, 'LENS_MANAGER[1]', 'LIBRARY_ELEMENT_UNIT_CELL');
Key = ltml.LTListByName(lt2, List, 'LibraryElement');
ltml.LTDbSet(lt2, Key, 'Filename', totalpathmod);
List = ltml.LTDbList(lt2, 'LENS_MANAGER[1]', 'TEXTURE_PARAMETER');
Key = ltml.LTListByName(lt2, List, 'StretchZ');
ltml.LTDbSet(lt2, Key, 'Value', 1);            % 높이는 .ent 에 이미 구워짐

%% CPS  (objFcn_both 와 동일; dETL/dHTL 은 이 가족의 변수가 아니므로 기준값 고정)
dETL = 60;  dHTL = 60;    % 스택 두께 기준값 (하이퍼벡터에 포함되지 않음)
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

%% bottom reflectance  (objFcn_both 와 동일)
TMF_OLED_bottom_p=TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
TMF_OLED_bottom_s=TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
R_p_bottom=abs(TMF_OLED_bottom_p.r_p).^2;
R_s_bottom=abs(TMF_OLED_bottom_s.r_s).^2;
Reflectance=(R_p_bottom+R_s_bottom)/2;

%% Coating  (objFcn_both 와 동일)
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

%% 파장 루프  (objFcn_both 와 동일)
I_white=0.5*(CPS_result.I_sub_s+CPS_result.I_sub_p);
sin089=sind(0:89);
P_white=I_white.*repmat(sin089,wavelength_num,1);
weight_factor=sum(P_white,2);
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

% 임시 슈퍼셀 .ent 정리 (디스크 누적 방지; 시뮬 종료 후이므로 안전)
try, delete(params.outPath); catch, end
end


%% ===== LightTools 재시작 (pareto_front_freeform.m 과 동일) =====
%  swept 모델 인스턴스도 그대로 띄운다: ID_LT 추출이 tasklist 의 3번째 토큰
%  (= 두 번째 lt.exe) 을 가정하므로, 프로세스 개수를 바꾸면 인덱싱이 깨진다.
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

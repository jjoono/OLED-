% ============================================================
%  stress_inverted_mla.m
%
%  목적: "practical saturation" (효율과 모든 각도구간이 함께 오르고,
%        선택성이 Lambertian 분할 [0.117 0.296 0.337 0.220] 근처에 고정)이
%        MLA 패밀리를 바꿔도 재현되는지 검사하는 per-family 경량 stress test.
%
%  이 파일의 패밀리: INVERTED (concave) MLA
%   - 기준(reference) freeform 은 볼록(convex) 렌즈렛이 공기 쪽으로 돌출.
%   - 여기서는 동일한 13차원 파라미터화의 렌즈렛을 오목(concave) 대응형으로
%     바꾼다. 기본 경로는 3D texture 의 양각/음각(bump vs hole) 토글이며,
%     프로파일 자체를 뒤집는 대안도 둔다 (아래 [반전 구현 — 두 경로] 참조).
%
%  프로토콜 (per-family stress protocol, 고정):
%   1. N_RANDOM = 100 개의 무작위 valid 설계 (search fidelity).
%   2. 단일 목적(EQE_total) 최적화 1회:
%        surrogateopt (MaxFunctionEvaluations=60, MinSurrogatePoints=25)
%        + patternsearch polish (15 evals)
%        + 승자 고정밀 재평가 (N_FINAL_REP=3).
%   3. 출력:
%        stress_inverted_result.mat
%          EVAL_LOG = [x(1:13) | EQE_total | b0_20 | b20_40 | b40_60 | b60_80 | phase | w]
%          (10회 평가마다 incremental save -> crash-safe)
%        stress_inverted_check.png  (3-panel):
%          (a) b40_60 vs EQE_total 산점도 + 선형 fit (near-linear collapse 검사)
%          (b) 최적 설계의 선택성 S_j vs Lambertian 예측 (grouped bars)
%          (c) 구간별 상관 R(EQE_total, S_j) 막대 + freeform baseline
%              (+0.6 / +0.7 / +0.05 / -0.7) overlay
%        + 세 가지 saturation signature 재현 여부 텍스트 판정.
%   4. Phase 마커: EVAL_PHASE = 6 (inverted family), EVAL_W = -1 (전 평가 공통).
%      random / optimization / high-precision 구분은 EVAL_LOG 행 인덱스 경계
%      (idx_random_end, idx_opt_end) 로 저장한다.
%
%  [반전 구현 — 두 경로]
%   (A) INVERT_METHOD = 'texture'  [기본값, 권장]
%       .lts 모델의 3D texture 는 unit cell 형상을 표면에 양각(bump)으로 얹을지
%       음각(hole)으로 파낼지를 자체 속성으로 지정한다. 이 경로에서는 unit cell
%       프로파일을 기준(볼록) 그대로 만들고 그 속성만 hole 로 돌린다.
%       -> 13변수 설계공간과 기하 정의가 freeform 계열과 완전히 동일하므로
%          계열 간 비교가 가장 공정하고, 모델이 원래 지원하는 경로라 안전하다.
%       속성 이름/값은 스크립트 상단 TEX_RELIEF_PARAM / TEX_RELIEF_VALUE_HOLE 에
%       두었다. 레포의 검증된 스크립트에 'StretchZ' 외 텍스처 속성 선례가 없으므로
%       실제 이름은 probe_texture_keys.m 으로 한 번 확인한 뒤 적어 넣을 것.
%       설정에 실패하면 즉시 error 를 던진다 (볼록 형상으로 조용히 계산되는 것 방지).
%
%   (B) INVERT_METHOD = 'profile'  [대안]
%       텍스처 토글 이름을 확인하지 못했을 때 쓰는 경로. 제어점 준비 단계
%       (정규화 공간) 에서 높이를 뒤집는다:
%           z_max   = max(높이 제어점)            (apex 포함, y_i 는 최대 1.5)
%           z_inv_i = z_max - z_i
%       stretchZ 가 uniform scale 이므로 정규화 공간 반전은 stretch 후 반전과
%       동치다. 결과 형상: 이전 apex 높이(z_max)가 슬래브 윗면이 되고 그 슬래브에
%       같은 프로파일의 접시(dish)가 파인다 — 축 높이 z_max-1, 림 높이 z_max.
%       인접 unit cell 과는 림에서 만나므로 타일링도 연속적이다. 이 경로를 쓰면
%       첫 실행에서 3D 뷰로 오목 여부와 RepairEntities 의 림 닫힘을 확인할 것.
%       stretchZ 부호 반전 / sweep 축 미러는 쓰지 않았다: 음수 stretchZ 는 렌즈
%       재료를 기판 '아래'로 뒤집어 층 순서를 깨고, RepairEntities 와 제어점 왕복
%       검증이 검증된 적 없는 경로가 되기 때문이다.
%
%   [주의 — 재료 방향은 물리적으로 유지되어야 함]
%   기판(substrate) 아래 / 텍스처된 렌즈 재료 위 / 그 위 공기. 반전은 렌즈
%   재료 층 '내부'의 윗면 형상만 dome -> dish 로 바꾸며, 층 순서·재질·소스·
%   리시버·각도구간 판독은 기준 코드와 완전히 동일하다.
%
%   [GEOM_TOL 검증과의 관계]
%   왕복 검증은 "SetSweptProfilePoints 로 넣은 배열" 과 "LightTools 에서
%   다시 읽은 제어점" 을 비교한다. 반전은 넣기 '전에' 적용되므로 검증은
%   자동으로 변환된(TRANSFORMED) 제어점 기준으로 수행된다 — 별도 수정 불필요.
%   반경방향 overshoot 재스케일(xy/max_length)도 두 열을 함께 나누는 uniform
%   scale 이라 반전된 프로파일을 그대로 보존한다.
%
%  기반: pareto_front_freeform.m (검증된 LightTools 연동/기하/스택을 그대로 복사,
%        위 반전 변환 한 가지만 추가)
% ============================================================
clear;
%% For LightTools Connection
global ID_swept ID_LT ltml ltloc count eval_count restart_interval ...
       ray_nums_current wave_n_current EVAL_LOG EVAL_PHASE EVAL_W ...
       GEOM_TOL GEOM_MISMATCH_LOG REQUIRE_MONOTONIC_X ...
       INVERT_METHOD TEX_RELIEF_PARAM TEX_RELIEF_VALUE_HOLE

%% ===== 음각(hole) 구현 방식 =====
%  'texture'  : LightTools 3D texture 자체의 양각/음각(bump vs hole) 토글을 쓴다.
%               unit cell 형상은 기준(볼록) 프로파일 그대로 두고 텍스처가
%               표면에서 그 형상을 파낸다 — .lts 모델에서 원래 지원하는 경로.
%               [권장] 기하를 손대지 않으므로 freeform 계열과 형상 정의가 동일하다.
%  'profile'  : 제어점 공간에서 z -> z_max - z 로 프로파일 자체를 뒤집는다.
%               텍스처 토글의 DB 속성 이름을 못 찾은 경우의 대안.
INVERT_METHOD = 'texture';

%  [!] TEXTURE_PARAMETER 안에서 양각/음각을 지정하는 항목 이름과, 음각에 해당하는
%      값. 레포의 검증된 스크립트에는 'StretchZ' 외 텍스처 속성 선례가 없으므로,
%      실제 이름은 probe_texture_keys.m 을 한 번 돌려 확인한 뒤 여기에 적는다.
%      (GUI Database Browser 에서 직접 확인해도 된다.)
TEX_RELIEF_PARAM      = 'BumpOrHole';   % <- probe 결과로 교체할 것
TEX_RELIEF_VALUE_HOLE = 1;              % <- 음각에 해당하는 값 (0=bump, 1=hole 로 가정)

% [기하 검증 tolerance] LightTools 제어점 왕복 불일치 허용치.
GEOM_TOL = 1e-4;
GEOM_MISMATCH_LOG = [];   % 열: [mismatch, max_length, rescale_triggered]

% [수율 개선] x 좌표 단조성 요구 (기준 스크립트와 동일).
%   비단조 x 프로파일은 스플라인 overshoot -> 재스케일 -> 재설정 불일치 ->
%   거부(NaN) 를 만들므로 생성 단계에서 배제한다.
REQUIRE_MONOTONIC_X = true;
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

%% ===== Multi-fidelity (기준 스크립트와 동일) =====
WAVE_N_SEARCH = 10;      % 탐색용 파장 step
WAVE_N_FINAL  = 2;       % 검증용 파장 step
RAY_SEARCH    = 10000;
RAY_FINAL     = 50000;
N_FINAL_REP   = 3;

%% ===== Stress-test 예산 (프로토콜 고정값) =====
N_RANDOM        = 100;   % 무작위 valid 설계 개수
OPT_EVALS       = 60;    % surrogateopt 평가 예산
MIN_SURR_POINTS = 25;
N_SEED_VALID    = 20;    % surrogateopt valid 초기 시드
POLISH_EVALS    = 15;    % patternsearch polish 예산

RESULT_MAT = 'stress_inverted_result.mat';

%% ===== Freeform baseline (비교 대상 saturation signatures) =====
% 기준 freeform 패밀리에서 관측된 구간별 상관 R(EQE_total, S_j):
edges      = [0 20; 20 40; 40 60; 60 80];
band_names = {'0-20 deg','20-40 deg','40-60 deg','60-80 deg'};
S_lamb     = sind(edges(:,2)).^2 - sind(edges(:,1)).^2;   % [0.117 0.296 0.337 0.220]
R_baseline = [0.6, 0.7, 0.05, -0.7];                      % freeform 관측치

%% Optimization Variables (13-dim, 기준과 동일)
varNames = {'x2','x3','x4','x5','x6', 'y2','y3','y4','y5','y6', 'dETL','dHTL','stretchZ'};
lb = [0, 0, 0, 0, 0, 0,   0,   0,   0,   0,   10, 10, 0.1];
ub = [1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 1.5, 1.5, 150,150, 3];
nvar = numel(lb);

EVAL_LOG = [];           % [x(1:13) | EQE_total | b0_20 | b20_40 | b40_60 | b60_80 | phase | w]
EVAL_PHASE = 6;          % inverted-family 마커 (전 평가 공통)
EVAL_W     = -1;

psOpts = optimoptions('patternsearch', ...
    'MaxFunctionEvaluations', POLISH_EVALS, ...
    'InitialMeshSize', 0.1, 'MeshTolerance', 1e-3, ...
    'Cache', 'on', 'Display', 'off');

%% =====================================================================
%  STEP 1 — 무작위 valid 설계 100개 (search fidelity)
%% =====================================================================
fprintf('\n########## INVERTED MLA STRESS — STEP 1: random valid designs (N=%d) ##########\n', N_RANDOM);
ray_nums_current = RAY_SEARCH;  wave_n_current = WAVE_N_SEARCH;
eval_count = 0;
Prand = genValidPoints(N_RANDOM, lb, ub);
for i = 1:N_RANDOM
    [et, eb] = simulate_both(Prand(i,:));
    if mod(i,10)==0
        fprintf('  random %3d/%d : EQE_total=%.4f  b40_60=%.4f\n', i, N_RANDOM, et, eb);
    end
end
idx_random_end = size(EVAL_LOG,1);
save(RESULT_MAT,'EVAL_LOG','varNames','lb','ub','idx_random_end');

% --- 기하 거부 진단 (NaN 의 주원인 판별) ---
report_geom_rejection(GEOM_MISMATCH_LOG, GEOM_TOL);

%% =====================================================================
%  STEP 2 — 단일 목적 EQE_total 최적화 (surrogateopt + polish)
%% =====================================================================
fprintf('\n########## STEP 2: single-objective EQE_total optimization ##########\n');
RenewLightTools();
lt = ltloc.GetLTAPI(ID_swept);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
ray_nums_current = RAY_SEARCH;  wave_n_current = WAVE_N_SEARCH;

seedMat = genValidPoints(N_SEED_VALID, lb, ub);
sopts = optimoptions('surrogateopt', ...
    'MaxFunctionEvaluations', OPT_EVALS, ...
    'MinSurrogatePoints',     MIN_SURR_POINTS, ...
    'InitialPoints',          struct('X', seedMat), ...
    'UseParallel', false, 'PlotFcn', [], 'Display', 'iter');
[xS, ~] = surrogateopt(@(x) scalar_objconstr_total(x), lb, ub, sopts);

% 국소 정련 (patternsearch, 15 evals)
if ~isempty(xS) && isValidPoints(xS(:).')
    x0 = xS(:).';
    try
        xP = patternsearch(@(x) scalar_polish_total(x), ...
                           x0, [],[],[],[], lb, ub, [], psOpts);
        xP = xP(:).';
    catch
        xP = x0;
    end
else
    error('surrogateopt 가 feasible 해를 반환하지 못했습니다. 시드/제약을 확인하세요.');
end

% 후보 중 EQE_total 이 더 좋은 쪽 채택 (저정밀 기준)
cands = {x0};
if ~isequal(xP, x0), cands{end+1} = xP; end
bestScore = -inf; bestX = [];
for c = 1:numel(cands)
    [et, ~] = simulate_both(cands{c});
    if ~isfinite(et), continue; end
    if et > bestScore, bestScore = et; bestX = cands{c}; end
end
if isempty(bestX), error('후보 재평가가 모두 NaN. LightTools 상태를 확인하세요.'); end
idx_opt_end = size(EVAL_LOG,1);
save(RESULT_MAT,'EVAL_LOG','varNames','lb','ub','idx_random_end','idx_opt_end','bestX');

%% =====================================================================
%  STEP 3 — 승자 고정밀 재평가 (N_FINAL_REP=3)
%% =====================================================================
fprintf('\n########## STEP 3: high-precision re-evaluation of winner ##########\n');
ray_nums_current = RAY_FINAL;  wave_n_current = WAVE_N_FINAL;
best_tot_r  = nan(1,N_FINAL_REP);
best_bins_r = nan(N_FINAL_REP,4);
for r = 1:N_FINAL_REP
    [best_tot_r(r), ~] = simulate_both(bestX);
    best_bins_r(r,:) = EVAL_LOG(end, nvar+2 : nvar+5);
end
best_tot  = mean(best_tot_r, 'omitnan');
best_bins = mean(best_bins_r, 1, 'omitnan');
best_S    = best_bins / best_tot;      % 최적 설계의 선택성 S_j
fprintf('  >>> best (inverted): EQE_total=%.5f  bins=[%.4f %.4f %.4f %.4f]\n', ...
    best_tot, best_bins);

save(RESULT_MAT,'EVAL_LOG','varNames','lb','ub','idx_random_end','idx_opt_end', ...
     'bestX','best_tot','best_bins','best_S','S_lamb','R_baseline','edges');

%% =====================================================================
%  분석 — 세 가지 saturation signature
%% =====================================================================
nv = nvar;
Et   = EVAL_LOG(:,nv+1);              % EQE_total
Bins = EVAL_LOG(:, nv+2 : nv+5);      % [0-20, 20-40, 40-60, 60-80]
ok   = isfinite(Et) & Et > 0.05 & all(isfinite(Bins),2);   % 저효율 설계는 노이즈 지배
fprintf('\n유효 설계 %d개 (전체 %d개 중)\n', sum(ok), size(EVAL_LOG,1));

% --- Signature 1: b40_60 vs EQE_total 의 near-linear collapse ---
p_fit = polyfit(Et(ok), Bins(ok,3), 1);
c1 = corrcoef(Et(ok), Bins(ok,3));  R_collapse = c1(1,2);

% --- Signature 3: 구간별 상관 R(EQE_total, S_j) ---
R_corr = nan(1,4);  S_meas = nan(1,4);  S_std = nan(1,4);
for b = 1:4
    s = Bins(ok,b) ./ Et(ok);
    S_meas(b) = mean(s,'omitnan');
    S_std(b)  = std(s,'omitnan');
    c = corrcoef(Et(ok), s);  R_corr(b) = c(1,2);
end

% --- 판정 ---
sig1 = R_collapse > 0.95;                                  % near-linear collapse
sig2 = all(abs(best_S(:) - S_lamb(:)) ./ S_lamb(:) < 0.15); % 최적 설계 S_j가 Lambertian 15% 이내
sig3 = all(abs(R_corr(:) - R_baseline(:)) < 0.25);          % 구간별 상관이 freeform baseline 재현

fprintf('\n################ SATURATION VERDICT — INVERTED (concave) MLA ################\n');
fprintf('[Sig 1] b40_60 vs EQE_total 선형 붕괴: R = %+.3f  (기준 >0.95) : %s\n', ...
    R_collapse, ternary(sig1,'재현','불일치'));
fprintf('[Sig 2] 최적 설계 선택성 vs Lambertian [%.3f %.3f %.3f %.3f]:\n', S_lamb);
for b = 1:4
    fprintf('        %-10s S=%.3f (예측 %.3f, 편차 %+.1f%%)\n', band_names{b}, ...
        best_S(b), S_lamb(b), 100*(best_S(b)-S_lamb(b))/S_lamb(b));
end
fprintf('        15%% 이내: %s\n', ternary(sig2,'재현','불일치'));
fprintf('[Sig 3] 구간별 R(EQE_total,S_j) vs freeform baseline (+0.6/+0.7/+0.05/-0.7):\n');
for b = 1:4
    fprintf('        %-10s R=%+.2f (baseline %+.2f)\n', band_names{b}, R_corr(b), R_baseline(b));
end
fprintf('        |dR|<0.25 전 구간: %s\n', ternary(sig3,'재현','불일치'));
if sig1 && sig2 && sig3
    fprintf('\n=> 세 signature 모두 재현. practical saturation 은 inverted(concave)\n');
    fprintf('   패밀리로 일반화된다 (각도 조향 없음, 효율-구간 동반 상승).\n');
else
    fprintf('\n=> 일부 signature 불일치 (%d/3 재현). inverted 패밀리에서 saturation\n', ...
        sum([sig1 sig2 sig3]));
    fprintf('   가설이 깨졌거나 표본/노이즈 문제. 위 세부 수치를 검토할 것.\n');
end

%% =====================================================================
%  3-panel 그림 — stress_inverted_check.png
%% =====================================================================
figure('Name','Inverted MLA stress check','Color','w','Position',[80 80 1350 420]);

subplot(1,3,1);   % (a) near-linear collapse
scatter(Et(ok), Bins(ok,3), 16, [.35 .45 .75], 'filled', 'MarkerFaceAlpha', 0.5); hold on;
xf = linspace(min(Et(ok)), max(Et(ok)), 50);
plot(xf, polyval(p_fit, xf), 'r-', 'LineWidth', 1.8);
xlabel('EQE_{total}'); ylabel('EQE_{40-60}'); grid on;
title(sprintf('(a) collapse check  R=%+.3f', R_collapse));
legend({'designs', sprintf('fit: %.3f x %+.4f', p_fit(1), p_fit(2))}, ...
    'Location','northwest','FontSize',8);

subplot(1,3,2);   % (b) best-design selectivity vs Lambertian
bar([best_S(:), S_lamb(:)], 'grouped'); grid on;
set(gca,'XTickLabel', band_names);
ylabel('selectivity  S_j = EQE_{band}/EQE_{total}');
legend({'best inverted design','Lambertian'},'Location','northwest','FontSize',8);
title('(b) best design vs Lambertian partition');

subplot(1,3,3);   % (c) per-band correlation vs freeform baseline
bar([R_corr(:), R_baseline(:)], 'grouped'); grid on;
set(gca,'XTickLabel', band_names);  ylim([-1 1]);  yline(0,'k-');
ylabel('R( EQE_{total}, S_j )');
legend({'inverted family','freeform baseline'},'Location','southwest','FontSize',8);
title('(c) per-band correlation');

saveas(gcf,'stress_inverted_check.png');
fprintf('\nsaved -> stress_inverted_check.png\n');

save(RESULT_MAT,'EVAL_LOG','varNames','lb','ub','idx_random_end','idx_opt_end', ...
     'bestX','best_tot','best_bins','best_S','S_lamb','R_baseline','edges', ...
     'R_collapse','p_fit','R_corr','S_meas','S_std','sig1','sig2','sig3');
report_geom_rejection(GEOM_MISMATCH_LOG, GEOM_TOL);
fprintf('\n########## 완료 ##########\n');
fprintf('  %s / stress_inverted_check.png\n', RESULT_MAT);


%% ===== 평가 래퍼 (네 각도구간을 모두 로그에 기록, crash-safe 저장 포함) =====
%  EVAL_LOG 열 구성:
%    [ x(1:13) | EQE_total | b0_20 | b20_40 | b40_60 | b60_80 | phase | w ]
%  phase = EVAL_PHASE = 6 (inverted family), w = EVAL_W = -1.
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
% [crash-safe] 10회 평가마다 incremental save.
%   파일이 이미 있으면 -append 로 EVAL_LOG 만 갱신 -> 앞선 milestone save 가
%   기록한 메타데이터(varNames/lb/ub/idx_*)를 지우지 않는다.
if mod(size(EVAL_LOG,1), 10) == 0
    try
        if isfile('stress_inverted_result.mat')
            save('stress_inverted_result.mat','EVAL_LOG','-append');
        else
            save('stress_inverted_result.mat','EVAL_LOG');
        end
    catch
    end
end
end

%% ===== 단일 목적 (surrogateopt: 제약 결합형, EQE_total 최대화) =====
function out = scalar_objconstr_total(x)
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
[et, ~] = simulate_both(x);
if ~isfinite(et)
    out.Ineq = 1;  out.Fval = 1;
else
    out.Ineq = -1;
    out.Fval = -et;                          % 최대화 -> 부호 반전
end
end

%% ===== 단일 목적 (patternsearch polish) =====
function f = scalar_polish_total(x)
global REQUIRE_MONOTONIC_X
x = x(:).';
if ~isempty(REQUIRE_MONOTONIC_X) && REQUIRE_MONOTONIC_X && any(diff(x(1:5)) < 0)
    f = 0; return;
end
if ~isValidPoints(x), f = 0; return; end
[et, ~] = simulate_both(x);
if ~isfinite(et), f = 0; return; end
f = -et;
end

%% ===== 판정 문구 헬퍼 =====
function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

%% ===== 기하 거부 진단 (기준 스크립트와 동일) =====
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

%% ===== 무작위 valid 시드 생성 (기준 스크립트와 동일) =====
%  REQUIRE_MONOTONIC_X 가 true 면 x2..x6 를 오름차순으로 정렬해 생성한다.
%  [주의] 제약은 반전 '전'의 원(原) 프로파일 다각형에 걸린다 — 13차원
%  파라미터화(설계 공간)는 freeform 기준과 완전히 동일해야 패밀리 간 비교가
%  성립하기 때문이다. 반전은 objFcn_both 내부의 기하 구성 단계에서만 적용된다.
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


%% ===== Objective (EQE_total + 네 각도구간; INVERTED profile) =====
%  기준 objFcn_both 를 그대로 복사하고, 제어점 준비 단계에 프로파일 반전
%  (z -> z_max - z) 한 가지만 추가했다. 소스/기판/리시버/코팅/CPS/각도구간
%  판독은 전부 동일.
function output = objFcn_both(point)
global ID_LT ID_swept ltml ltloc count ray_nums_current wave_n_current ...
       INVERT_METHOD TEX_RELIEF_PARAM TEX_RELIEF_VALUE_HOLE
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

% [반전 — 이 패밀리의 유일한 물리 변경]
%   INVERT_METHOD = 'texture' 일 때는 여기서 프로파일을 건드리지 않는다.
%   unit cell 은 기준(볼록) 형상 그대로 만들어지고, 아래쪽 텍스처 설정에서
%   LightTools 의 양각/음각(bump vs hole) 토글로 표면에서 파내진다.
%   -> 설계공간(13변수)과 기하 정의가 freeform 계열과 완전히 동일해지므로
%      계열 간 비교가 가장 공정하다.
%
%   INVERT_METHOD = 'profile' 은 텍스처 토글의 DB 속성 이름을 확인하지 못한
%   경우의 대안이다: 정규화 제어점 공간에서 z -> z_max - z 로 프로파일 자체를
%   뒤집어, 이전 apex 높이가 윗면이 되는 슬래브에 같은 프로파일의 접시가
%   파인 형상을 만든다. stretchZ 는 uniform scale 이라 정규화 공간 반전과
%   순서 교환이 가능하고, 아래의 재스케일/GEOM_TOL 왕복 검증은 변환된 xy 에
%   대해 그대로 수행된다. 이 경로를 쓸 때는 첫 실행에서 3D 뷰로 오목 여부와
%   RepairEntities 의 림 닫힘을 반드시 확인할 것.
if strcmpi(INVERT_METHOD, 'profile')
    z_max = max(xy(:,2));
    xy(:,2) = z_max - xy(:,2);
end

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
    % uniform scale 이므로 반전된 프로파일 형상은 그대로 보존된다.
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
%   여기의 xy 는 이미 '반전된(TRANSFORMED)' 제어점이므로, 왕복 검증은
%   변환된 형상 기준으로 수행된다 (기준 스크립트와 동일한 tol/재시도 로직).
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

% [음각 설정] 텍스처의 양각/음각 토글을 hole 로 돌린다.
%   StretchZ 와 동일한 TEXTURE_PARAMETER 목록/설정 패턴을 그대로 쓴다.
%   속성 이름을 못 찾으면 조용히 넘어가지 않고 즉시 멈춘다 — 그대로 두면
%   볼록(bump) 형상으로 계산이 끝나 freeform 계열과 구분되지 않기 때문이다.
if strcmpi(INVERT_METHOD, 'texture')
    try
        KeyRelief = ltml.LTListByName(lt2, List, TEX_RELIEF_PARAM);
        ltml.LTDbSet(lt2, KeyRelief, 'Value', TEX_RELIEF_VALUE_HOLE);
    catch ME
        error(['텍스처 음각 파라미터 ''%s'' 설정 실패: %s\n' ...
               'probe_texture_keys.m 을 실행해 실제 속성 이름과 음각 값을 확인한 뒤 ' ...
               '스크립트 상단의 TEX_RELIEF_PARAM / TEX_RELIEF_VALUE_HOLE 을 고칠 것. ' ...
               '(대안: INVERT_METHOD = ''profile'')'], TEX_RELIEF_PARAM, ME.message);
    end
end

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
% [주의] 파장 샘플 인덱스는 명시적으로 생성 (division-based 크기 계산 금지).
%   wv_list = 1:n:wavelength_num 은 파장창과 n 의 조합에 무관하게 안전하다.
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

%% ===== Spline 제약 (기준과 동일 — 반전 '전' 원 프로파일에 적용) =====
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

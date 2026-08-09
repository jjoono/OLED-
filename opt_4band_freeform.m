% ============================================================
%  opt_4band_freeform.m
%
%  목적: 원고 abstract 의 "each of the four angular bands was independently
%        optimized" 주장을 데이터로 채우기 위한 네 구간 단일목적 최적화.
%        (40-60도는 pareto_front_freeform.m 의 w=0 에서도 수행됐지만, 동일 조건
%         비교를 위해 여기서 한 번에 전부 다시 돌린다)
%
%  [실행 내용] 네 구간 + 총 EQE 대조군을 각각 독립적으로 최대화:
%        0-20도, 20-40도, 40-60도, 60-80도   (목적 = -EQE_band, 단일목적)
%        + EQE_total (control arm)           (목적 = -EQE_total)
%    대조군을 같은 스크립트에서 같은 예산/정밀도/기하로 돌리는 이유:
%    "구간 전용 최적화가 무조향 대비 이득인가" 의 기준선을 외부 숫자에 의존하지
%    않게 하기 위함. 다른 스크립트의 max EQE_total 을 기준으로 쓰면 예산,
%    patch 크기(x_pattern/y_pattern), fidelity 차이가 결론에 섞인다.
%    각 구간마다:
%      (i)  surrogateopt 탐색 (저정밀: RAY_SEARCH / WAVE_N_SEARCH, 예산 60회)
%      (ii) patternsearch 국소 정련 (예산 15회)
%      (iii) 승자 고정밀 재평가 (RAY_FINAL / WAVE_N_FINAL, N_FINAL_REP회 반복 평균)
%
%  [warm start] pareto_front_result.mat 이 pwd 에 있으면 그 EVAL_LOG 에서
%        해당 구간 EQE 상위 10개 설계를 surrogateopt InitialPoints 시드로 재사용.
%        파일이 없어도 무작위 valid 시드만으로 그대로 동작한다 (standalone).
%
%  [수집물] EVAL_LOG (pareto_front_freeform.m 과 동일한 열 구성):
%        [ x(1:13) | EQE_total | EQE_0_20 | EQE_20_40 | EQE_40_60 | EQE_60_80 | phase | w ]
%     phase 4 = 이 스크립트 표식,  w 열 = 구간 하한 각도 (0/20/40/60), 대조군 = -1.
%     구간이 하나 끝날 때마다 opt_4band_result.mat 에 저장 (crash-safe).
%     최종 요약은 opt_4band_result_merged.mat 에 저장한다 (건너뛴 arm 을 이전
%     결과에서 병합하므로, 원본 opt_4band_result.mat 을 덮어쓰지 않는다).
%
%  [부분 실행] ARMS_TO_RUN 으로 일부 arm 만 돌릴 수 있다. 예: 네 band 결과가
%     이미 있고 대조군만 없으면  ARMS_TO_RUN = CTRL_IDX;  (약 1/5 시간)
%
%  기반: pareto_front_freeform.m (동일한 LightTools 연동/기하/스택/제약)
% ============================================================
clear;
%% For LightTools Connection
global ID_swept ID_LT ltml ltloc count eval_count restart_interval ...
       ray_nums_current wave_n_current EVAL_LOG EVAL_PHASE EVAL_W ...
       GEOM_TOL GEOM_MISMATCH_LOG REQUIRE_MONOTONIC_X PATCH_XY

%% ===== 텍스처 패치 크기 (mm) =====
%  다른 스크립트와 반드시 같아야 하는 값. 서로 다르면 EQE_total 이 통째로
%  달라져 스크립트 간 비교가 무의미해진다 (repo 기준: pareto_front_freeform.m
%  = 15, convergence_check.m = 25).
%  Fig.2/Fig.3 데이터를 만든 pareto 실행과 반드시 일치시킬 것.
PATCH_XY = 15;

% [기하 검증 tolerance] LightTools 제어점 왕복 불일치 허용치.
GEOM_TOL = 1e-4;
GEOM_MISMATCH_LOG = [];   % 열: [mismatch, max_length, rescale_triggered]

% [수율 개선] x 좌표 단조성 요구 (pareto_front_freeform.m 과 동일).
%   비단조 x2..x6 프로파일은 스플라인 overshoot -> 재스케일 -> 제어점 불일치로
%   어차피 거부(NaN)되므로 생성 단계에서 배제한다 (NaN 40% -> 수율 대폭 개선).
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

%% ===== Multi-fidelity (pareto_front_freeform.m 과 동일) =====
WAVE_N_SEARCH = 10;      % 탐색용 파장 step
WAVE_N_FINAL  = 2;       % 검증용 파장 step
RAY_SEARCH    = 10000;
RAY_FINAL     = 50000;
N_FINAL_REP   = 3;

%% ===== 구간별 단일목적 최적화 설정 =====
%  네 구간 + 총 EQE 대조군(control arm)을 한 스크립트에서 같은 예산/정밀도/기하로
%  돌린다. 대조군이 있어야 "구간 전용 최적화가 무조향 대비 이득인가" 를 외부
%  숫자(다른 스크립트의 max EQE_total) 없이 내부적으로 판정할 수 있다 —
%  예산, patch 크기(x_pattern/y_pattern), fidelity 차이가 결론에 섞이지 않는다.
BAND_LIST  = {[0 20], [20 40], [40 60], [60 80], [0 90]};
BAND_NAMES = {'0-20 deg', '20-40 deg', '40-60 deg', '60-80 deg', 'EQE_total (control)'};
BAND_COL   = [1, 2, 3, 4, 0];   % bins 열 인덱스. 0 = 총 EQE 대조군
CTRL_IDX   = find(BAND_COL == 0, 1);
nBand = numel(BAND_LIST);

%% ===== 어느 arm 을 돌릴지 =====
%  한 번 다 돌리는 데 오래 걸리므로, 이미 결과가 있는 arm 은 건너뛸 수 있게 한다.
%    []            -> 전부 (기본)
%    CTRL_IDX      -> 대조군만 (band 결과를 이미 갖고 있을 때. 약 1/5 시간)
%    [1 3]         -> 해당 arm 만
%  건너뛴 arm 은 SKIP_MERGE_FILE 의 이전 결과로 채워 요약표/그림을 그대로 낸다.
ARMS_TO_RUN     = [];
SKIP_MERGE_FILE = 'opt_4band_result.mat';   % 건너뛴 arm 을 가져올 이전 결과 파일

EVALS_PER_BAND  = 60;    % 구간당 surrogateopt 평가 예산
MIN_SURR_POINTS = 25;
N_SEED_VALID    = 15;    % 무작위 valid 시드 (warm start 10개와 합쳐 예산 내 유지)
POLISH_EVALS    = 15;
N_WARM          = 10;    % warm start 로 가져올 상위 설계 수

WARMSTART_FILE = 'pareto_front_result.mat';   % pareto_front_freeform.m 이 저장하는 이름

%% ===== 정규화 상수 =====
% 단일목적이라 정규화가 최적점을 바꾸지는 않지만, surrogateopt 의 Fval 스케일을
% O(1) 로 맞춰준다. 초기값 = Lambertian 선택성 x EQE_total 최고치(0.56) 추정.
% warm start 로그가 있으면 해당 구간 최댓값의 1.1배로 자동 갱신.
S_LAMB_ALL = [0.117, 0.296, 0.337, 0.220];    % [0-20, 20-40, 40-60, 60-80] Lambertian 선택성
MAX_EQE_TOTAL_REF = 0.56;                     % 정규화 스케일용 초기 추정치일 뿐.
                                              % no-steering 기준은 대조군 결과로 대체된다.
REF_BAND_LIST = nan(1, nBand);
for b = 1:nBand
    if BAND_COL(b) == 0
        REF_BAND_LIST(b) = MAX_EQE_TOTAL_REF;                       % 대조군: 총 EQE
    else
        REF_BAND_LIST(b) = S_LAMB_ALL(BAND_COL(b)) * MAX_EQE_TOTAL_REF;
    end
end

%% Optimization Variables (13-dim, 기존과 동일)
varNames = {'x2','x3','x4','x5','x6', 'y2','y3','y4','y5','y6', 'dETL','dHTL','stretchZ'};
lb = [0, 0, 0, 0, 0, 0,   0,   0,   0,   0,   10, 10, 0.1];
ub = [1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 1.5, 1.5, 150,150, 3];
nvar = numel(lb);

EVAL_LOG = [];           % [x(1:13), EQE_total, b0_20, b20_40, b40_60, b60_80, phase, w]

psOpts = optimoptions('patternsearch', ...
    'MaxFunctionEvaluations', POLISH_EVALS, ...
    'InitialMeshSize', 0.1, 'MeshTolerance', 1e-3, ...
    'Cache', 'on', 'Display', 'off');

%% ===== Warm start 로그 로드 (있으면) =====
%  pareto_front_freeform.m 의 EVAL_LOG 는 이 스크립트와 동일한 열 구성이므로
%  구간별 EQE 열을 그대로 랭킹에 쓸 수 있다. 파일이 없으면 그냥 건너뛴다.
WARM_LOG = [];
if exist(WARMSTART_FILE, 'file')
    Dw = load(WARMSTART_FILE);
    if isfield(Dw, 'EVAL_LOG') && size(Dw.EVAL_LOG,2) >= nvar+5
        WARM_LOG = Dw.EVAL_LOG;
        fprintf('[Warm start] %s 로드: EVAL_LOG %d행\n', WARMSTART_FILE, size(WARM_LOG,1));
        % 정규화 상수 자동 보정 (구간별)
        for b = 1:nBand
            if BAND_COL(b) == 0
                mb = max(WARM_LOG(:, nvar+1), [], 'omitnan');        % 총 EQE 열
            else
                mb = max(WARM_LOG(:, nvar+1+BAND_COL(b)), [], 'omitnan');
            end
            if isfinite(mb) && mb > 0, REF_BAND_LIST(b) = 1.1*mb; end
        end
    end
else
    fprintf('[Warm start] %s 없음 -> 무작위 시드만 사용 (standalone 실행)\n', WARMSTART_FILE);
end

%% ===== 결과 컨테이너 =====
band_x           = nan(nBand, nvar);   % 구간별 최적 설계
band_eqe_coarse  = nan(nBand, 1);      % 저정밀(탐색 fidelity) 구간 EQE
band_eqe_hi      = nan(nBand, 1);      % 고정밀 구간 EQE (N_FINAL_REP 평균)
band_tot_hi      = nan(nBand, 1);      % 고정밀 EQE_total
band_bins_hi     = nan(nBand, 4);      % 고정밀 [b0_20 b20_40 b40_60 b60_80]

%% ===== 건너뛸 arm 을 이전 결과에서 채운다 =====
if isempty(ARMS_TO_RUN)
    ARMS_TO_RUN = 1:nBand;
end
skipArms = setdiff(1:nBand, ARMS_TO_RUN);
if ~isempty(skipArms)
    if ~exist(SKIP_MERGE_FILE, 'file')
        error(['arm %s 을(를) 건너뛰려면 이전 결과 파일 %s 이 필요하다. ' ...
               'ARMS_TO_RUN 을 [] 로 두고 전부 돌리거나 파일을 준비할 것.'], ...
               mat2str(skipArms), SKIP_MERGE_FILE);
    end
    Dp = load(SKIP_MERGE_FILE);
    for k = skipArms
        if isfield(Dp,'band_eqe_hi') && numel(Dp.band_eqe_hi) >= k && isfinite(Dp.band_eqe_hi(k))
            band_x(k,:)       = Dp.band_x(k,:);
            band_eqe_coarse(k)= Dp.band_eqe_coarse(k);
            band_eqe_hi(k)    = Dp.band_eqe_hi(k);
            band_tot_hi(k)    = Dp.band_tot_hi(k);
            band_bins_hi(k,:) = Dp.band_bins_hi(k,:);
            fprintf('[Skip] %-20s : 이전 결과 사용 (EQE_band=%.5f, EQE_total=%.5f)\n', ...
                BAND_NAMES{k}, band_eqe_hi(k), band_tot_hi(k));
        else
            fprintf('[Skip] %-20s : 이전 결과에 없음 -> NaN 유지\n', BAND_NAMES{k});
        end
    end
    % 이전 EVAL_LOG 도 이어 붙여 산점도/진단이 전체 표본을 보게 한다.
    if isfield(Dp,'EVAL_LOG') && size(Dp.EVAL_LOG,2) == nvar+7
        EVAL_LOG = Dp.EVAL_LOG;
        fprintf('[Skip] 이전 EVAL_LOG %d행 이어받음\n', size(EVAL_LOG,1));
    end
end

%% =====================================================================
%  구간 루프 — 각 구간을 독립 단일목적으로 최적화
%% =====================================================================
for k = ARMS_TO_RUN(:).'
    band  = BAND_LIST{k};
    bcol  = BAND_COL(k);
    if bcol == 0
        wtag = -1;                     % 대조군 표식 (구간 하한과 겹치지 않게)
    else
        wtag = band(1);                % w 열에 기록할 구간 하한 (0 / 20 / 40 / 60)
    end
    refB  = REF_BAND_LIST(k);
    fprintf('\n########## ARM %s  (%d/%d) ##########\n', BAND_NAMES{k}, k, nBand);
    if bcol == 0
        fprintf('  목적 = maximize EQE_total  (무조향 대조군, 정규화 ref %.4f)\n', refB);
    else
        fprintf('  목적 = maximize EQE_%d_%d  (정규화 ref %.4f)\n', band(1), band(2), refB);
    end

    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);
    ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    eval_count = 0;
    EVAL_PHASE = 4;  EVAL_W = wtag;
    ray_nums_current = RAY_SEARCH;  wave_n_current = WAVE_N_SEARCH;

    % --- 시드: warm start 상위 N_WARM + 무작위 valid ---
    seedMat = genValidPoints(N_SEED_VALID, lb, ub);
    if ~isempty(WARM_LOG)
        Xw = WARM_LOG(:, 1:nvar);
        if bcol == 0
            Bw = WARM_LOG(:, nvar+1);            % 대조군: 총 EQE 열
        else
            Bw = WARM_LOG(:, nvar+1+bcol);       % 현재 구간의 EQE 열
        end
        okw = isfinite(Bw) & isValidPoints(Xw);
        if REQUIRE_MONOTONIC_X                    % 비단조 x 는 objconstr 에서 거부되므로 제외
            okw = okw & all(diff(Xw(:,1:5), 1, 2) >= 0, 2);
        end
        Xw = Xw(okw,:);  Bw = Bw(okw);
        [Xw, iu] = unique(Xw, 'rows', 'stable');  Bw = Bw(iu);
        [~, is] = sort(Bw, 'descend');
        nTake = min(N_WARM, numel(is));
        if nTake > 0
            seedMat = [Xw(is(1:nTake), :); seedMat];
            fprintf('  warm start: 상위 %d개 시드 (최고 EQE_%d_%d = %.4f)\n', ...
                nTake, band(1), band(2), Bw(is(1)));
        end
    end

    sopts = optimoptions('surrogateopt', ...
        'MaxFunctionEvaluations', EVALS_PER_BAND, ...
        'MinSurrogatePoints',     MIN_SURR_POINTS, ...
        'InitialPoints',          struct('X', seedMat), ...
        'UseParallel', false, 'PlotFcn', [], 'Display', 'iter');
    [xS, ~] = surrogateopt(@(x) band_objconstr(x, bcol, refB), lb, ub, sopts);

    % 국소 정련
    if ~isempty(xS) && isValidPoints(xS(:).')
        x0 = xS(:).';
        try
            xP = patternsearch(@(x) band_polish(x, bcol, refB), ...
                               x0, [],[],[],[], lb, ub, [], psOpts);
            xP = xP(:).';
        catch
            xP = x0;
        end
    else
        fprintf('  [Warn] band %s: feasible 해 미반환\n', BAND_NAMES{k});
        continue;
    end

    % 후보 중 구간 EQE 가 더 좋은 쪽 채택 (저정밀 기준)
    cands = {x0};
    if ~isequal(xP, x0), cands{end+1} = xP; end
    bestScore = -inf; bestX = [];
    for c = 1:numel(cands)
        [et, bins] = simulate_bands(cands{c});
        if ~isfinite(et), continue; end
        if bcol == 0, sc = et; else, sc = bins(bcol); end
        if sc > bestScore, bestScore = sc; bestX = cands{c}; end
    end
    if isempty(bestX), continue; end
    band_eqe_coarse(k) = bestScore;

    % 고정밀 재평가 (N_FINAL_REP 반복 평균)
    ray_nums_current = RAY_FINAL;  wave_n_current = WAVE_N_FINAL;
    et_r   = nan(N_FINAL_REP, 1);
    bins_r = nan(N_FINAL_REP, 4);
    for r = 1:N_FINAL_REP
        [et_r(r), bins_r(r,:)] = simulate_bands(bestX);
    end
    band_x(k,:)      = bestX;
    band_tot_hi(k)   = mean(et_r, 'omitnan');
    band_bins_hi(k,:)= mean(bins_r, 1, 'omitnan');
    if bcol == 0
        band_eqe_hi(k) = band_tot_hi(k);        % 대조군의 '목적값' 은 총 EQE 자체
    else
        band_eqe_hi(k) = band_bins_hi(k, bcol);
    end
    fprintf('  >>> %s : EQE_band=%.5f (coarse %.5f)  EQE_total=%.5f  (선택성 %.1f%%)\n', ...
        BAND_NAMES{k}, band_eqe_hi(k), band_eqe_coarse(k), band_tot_hi(k), ...
        100*band_eqe_hi(k)/band_tot_hi(k));

    % crash-safe: 구간 하나 끝날 때마다 저장
    save('opt_4band_result.mat', 'EVAL_LOG', 'band_x', 'band_eqe_coarse', ...
         'band_eqe_hi', 'band_tot_hi', 'band_bins_hi', ...
         'BAND_LIST', 'BAND_NAMES', 'BAND_COL', 'REF_BAND_LIST', ...
         'varNames', 'lb', 'ub', 'GEOM_MISMATCH_LOG');
end

%% =====================================================================
%  결과 정리 — 요약표 + 그림
%% =====================================================================
% 구간별 Lambertian 선택성 (대조군 자리는 NaN)
S_LAMB_SEL = nan(nBand,1);
for k = 1:nBand
    if BAND_COL(k) ~= 0, S_LAMB_SEL(k) = S_LAMB_ALL(BAND_COL(k)); end
end

% no-steering 기준: 대조군(총 EQE 최적화)의 고정밀 EQE_total.
%   같은 스크립트/예산/정밀도/기하에서 나온 값이라 외부 숫자와 달리 교란이 없다.
if isfinite(band_tot_hi(CTRL_IDX))
    E_CTRL     = band_tot_hi(CTRL_IDX);
    CTRL_LABEL = sprintf('control EQE_{total} = %.4f', E_CTRL);
else
    E_CTRL     = MAX_EQE_TOTAL_REF;   % 대조군 실패 시에만 폴백
    CTRL_LABEL = sprintf('fallback %.2f (대조군 실패)', E_CTRL);
    fprintf('\n[Warn] 대조군 실패 -> no-steering 기준을 폴백 %.2f 로 사용\n', E_CTRL);
end

fprintf('\n################ 4-BAND + CONTROL SUMMARY ################\n');
fprintf('%-20s | %10s | %10s | %10s | %10s\n', ...
    'arm', 'EQE_band', 'EQE_total', '선택성', 'Lamb. S');
for k = 1:nBand
    if ~isfinite(band_eqe_hi(k))
        fprintf('%-20s | %10s | %10s | %10s | %10s\n', BAND_NAMES{k}, 'FAIL', '-', '-', '-');
    elseif BAND_COL(k) == 0
        fprintf('%-20s | %10s | %10.5f | %10s | %10s\n', ...
            BAND_NAMES{k}, '-', band_tot_hi(k), '-', '-');
    else
        fprintf('%-20s | %10.5f | %10.5f | %9.1f%% | %10.3f\n', ...
            BAND_NAMES{k}, band_eqe_hi(k), band_tot_hi(k), ...
            100*band_eqe_hi(k)/band_tot_hi(k), S_LAMB_SEL(k));
    end
end

% --- 이득 분해: 선택성 이득 x 총량 비 = 순 band power 비 ---
%   구간 전용 최적화가 무조향(대조군 + Lambertian 분배) 대비 실제로 이득인가?
fprintf('\n---- 이득 분해 (기준: %s) ----\n', CTRL_LABEL);
fprintf('%-20s | %8s | %8s | %8s | %8s\n', ...
    'arm', 'S/S_Lamb', 'E/E_ctrl', '순이득', '판정');
gain_sel = nan(nBand,1); gain_tot = nan(nBand,1); gain_net = nan(nBand,1);
for k = 1:nBand
    if BAND_COL(k) == 0 || ~isfinite(band_eqe_hi(k)), continue; end
    Sk = band_eqe_hi(k)/band_tot_hi(k);
    gain_sel(k) = Sk / S_LAMB_SEL(k);
    gain_tot(k) = band_tot_hi(k) / E_CTRL;
    gain_net(k) = band_eqe_hi(k) / (S_LAMB_SEL(k) * E_CTRL);
    if gain_net(k) > 1.05,      verdict = '이득';
    elseif gain_net(k) < 0.95,  verdict = '손해';
    else,                       verdict = '무이득';
    end
    fprintf('%-20s | %8.3f | %8.3f | %8.3f | %8s\n', ...
        BAND_NAMES{k}, gain_sel(k), gain_tot(k), gain_net(k), verdict);
end
fprintf(['\n순이득 = (전용 최적화 band EQE) / (Lambertian 선택성 x 대조군 EQE_total)\n' ...
         '  > 1 이면 조향이 순이득, ~1 이면 선택성 이득이 총량 손실과 상쇄,\n' ...
         '  < 1 이면 구간 전용 최적화가 오히려 손해.\n']);

% --- 그림: 달성치 vs no-steering 예측 (대조군 기준) ---
isBand  = BAND_COL(:) ~= 0;
noSteer = S_LAMB_SEL(isBand) * E_CTRL;
figure('Name','4-band independent optimization vs control','Color','w', ...
       'Position',[120 120 760 420]);
bar([band_eqe_hi(isBand), noSteer(:)]);
set(gca, 'XTickLabel', BAND_NAMES(isBand));
ylabel('EQE_{band}');
legend({'achieved (independent opt.)', ...
        sprintf('no-steering  S_{Lamb} x %.4f (control)', E_CTRL)}, ...
       'Location','best', 'FontSize', 9);
title('Band-EQE: independent optimization vs no-steering control');
grid on;
saveas(gcf, 'opt_4band_summary.png');
fprintf('saved -> opt_4band_summary.png\n');

save('opt_4band_result_merged.mat', 'EVAL_LOG', 'band_x', 'band_eqe_coarse', ...
     'band_eqe_hi', 'band_tot_hi', 'band_bins_hi', ...
     'BAND_LIST', 'BAND_NAMES', 'BAND_COL', 'REF_BAND_LIST', ...
     'S_LAMB_ALL', 'MAX_EQE_TOTAL_REF', 'CTRL_IDX', 'E_CTRL', 'PATCH_XY', ...
     'gain_sel', 'gain_tot', 'gain_net', ...
     'varNames', 'lb', 'ub', 'GEOM_MISMATCH_LOG');
fprintf('saved -> opt_4band_result_merged.mat  (EVAL_LOG %d points)\n', size(EVAL_LOG,1));
report_geom_rejection(GEOM_MISMATCH_LOG, GEOM_TOL);   % 전체 실행 기준 최종 진단
fprintf('\n########## 완료 ##########\n');


%% ===== 평가 래퍼 (네 각도구간을 모두 로그에 기록) =====
%  EVAL_LOG 열 구성 (pareto_front_freeform.m 과 동일):
%    [ x(1:13) | EQE_total | EQE_0_20 | EQE_20_40 | EQE_40_60 | EQE_60_80 | phase | w ]
%  phase = 4 (이 스크립트), w = 구간 하한 각도.
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

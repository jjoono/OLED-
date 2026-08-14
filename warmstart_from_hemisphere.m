% ============================================================
%  warmstart_from_hemisphere.m
%
%  목적: "13변수 탐색이 3변수 탐색에 진 것은 형상 자유도가 쓸모없어서가 아니라
%        탐색이 수렴하지 않아서다" 라는 심사 반론을 차단한다.
%
%  [왜 필요한가] opt_4band_freeform.m 결과에서 G_j = 0.946 / 0.966 / 1.070 /
%    1.020 로, 0-20도와 20-40도에서 freeform 이 반구보다 낮게 나왔다. 반구는
%    13변수 feasible set 의 한 점이므로 예산이 무한하면 G_j >= 1 이 자명하다.
%    즉 현재 결과는 "형상 자유도가 없다" 가 아니라 "탐색이 못 찾았다" 로도
%    읽힌다. 이 스크립트는 그 해석을 배제한다.
%
%  [무엇을 하는가] 각 arm 의 탐색을 **그 arm 의 반구 최적해에서 출발**시킨다.
%    반구 해가 이미 시드에 들어 있으므로, 탐색이 반환하는 값은 정의상
%    반구 이상이다. 따라서 결과는 둘 중 하나로만 읽힌다:
%      (a) 유의미하게 개선됨  -> 형상 자유도가 실재하고 원래 예산이 부족했다.
%                               이 경우 원고의 saturation 주장을 수정해야 한다.
%      (b) 노이즈 안에서 정체 -> 반구가 13차원 공간의 국소 최적점이다.
%                               "탐색 실패" 반론이 성립하지 않는다.
%    (b) 가 나오면 원고에서 가장 약한 문장(2.2절 "차원의 비용")을 가장 강한
%    문장("반구에서 출발시켜도 넘지 못한다")으로 교체할 수 있다.
%
%  [세 갈래로 민다] 한 가지 탐색기의 실패로 (b) 가 나오는 것을 막기 위해:
%      (i)   surrogateopt : 반구 해 + 그 근방 섭동점들을 InitialPoints 로
%      (ii)  patternsearch: surrogateopt 승자에서 국소 정련
%      (iii) patternsearch: **반구 해에서 직접** 국소 정련
%    (iii) 이 핵심이다. 반구 점에서 시작한 국소 탐색이 개선하지 못하면
%    그 점이 국소 최적이라는 가장 직접적인 증거가 된다.
%
%  [노이즈 처리] 개선 판정을 Monte-Carlo 노이즈와 구분해야 한다. 반구 해와
%    최종 승자를 **둘 다** 같은 고정밀 설정으로 N_FINAL_REP 회 반복 평가하고,
%    차이가 두 표준오차의 합보다 클 때만 '개선' 으로 판정한다. 반구 기준값을
%    opt_hemisphere_result.mat 에서 그냥 읽어오지 않고 여기서 다시 재는 이유는,
%    두 스크립트 사이의 사소한 설정 차이가 차이값에 섞이지 않게 하기 위함이다.
%    (읽어온 값과는 교차 검증만 한다.)
%
%  [수집물] warmstart_hemisphere_result.mat
%      base_val/base_tot/base_bins/base_sd : 반구 해의 고정밀 재측정
%      ws_val/ws_tot/ws_bins/ws_sd         : warm start 탐색 승자
%      ws_x                                 : 승자 설계 (13변수)
%      ws_src                               : 승자가 어느 갈래에서 나왔는지
%      verdict                              : arm 별 판정 문자열
%    EVAL_LOG 열 구성은 다른 캠페인과 동일. phase 9 = 이 스크립트 표식.
%
%  [실행 시간] arm 당 고정밀 6회 + 저정밀 약 70회. 5 arm 전부면
%    opt_4band_freeform.m 한 번과 비슷하거나 조금 짧다. 시간이 부족하면
%    ARMS_TO_RUN = [1 2] 로 G_j < 1 인 두 arm 만 돌려도 반론 차단에는 충분하다.
%
%  기반: opt_4band_freeform.m (동일한 LightTools 연동/기하/스택/제약/헬퍼)
% ============================================================
clear;
%% For LightTools Connection
global ID_swept ID_LT ltml ltloc count eval_count restart_interval ...
       ray_nums_current wave_n_current EVAL_LOG EVAL_PHASE EVAL_W ...
       GEOM_TOL GEOM_MISMATCH_LOG REQUIRE_MONOTONIC_X PATCH_XY

%% ===== 텍스처 패치 크기 (mm) =====
%  [통일 규칙] 모든 캠페인 스크립트에서 25 x 25 mm. 이 값이 다르면 반구 기준값과
%  비교 자체가 무의미해진다.
PATCH_XY = 25;

GEOM_TOL = 1e-4;
GEOM_MISMATCH_LOG = [];
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

%% ===== Multi-fidelity (다른 캠페인과 동일하게 유지할 것) =====
WAVE_N_SEARCH = 10;      % 탐색용 파장 step
WAVE_N_FINAL  = 2;       % 검증용 파장 step
RAY_SEARCH    = 10000;
RAY_FINAL     = 50000;
N_FINAL_REP   = 3;       % >= 2 여야 한다. 판정이 반복 간 표준편차에 의존한다.
if N_FINAL_REP < 2
    error('N_FINAL_REP 는 2 이상이어야 한다 (노이즈 대비 판정에 std 가 필요).');
end

%% ===== arm 정의 (opt_4band_freeform.m / opt_hemisphere_arms.m 과 동일 순서) =====
BAND_LIST  = {[0 20], [20 40], [40 60], [60 80], [0 90]};
BAND_NAMES = {'0-20 deg', '20-40 deg', '40-60 deg', '60-80 deg', 'EQE_total'};
BAND_COL   = [1, 2, 3, 4, 0];     % bins 열 인덱스. 0 = 총 EQE
nBand = numel(BAND_LIST);

%  [부분 실행] G_j < 1 인 arm 만 돌리려면 ARMS_TO_RUN = [1 2];
ARMS_TO_RUN = [];                 % [] = 전부

%% ===== 탐색 예산 =====
%  원래 캠페인(60+15)보다 적게 잡아도 된다. 시드가 이미 최적점 근방이므로
%  탐색의 역할이 '찾기' 가 아니라 '근방에 더 나은 것이 있는지 확인' 이다.
EVALS_PER_BAND  = 40;
MIN_SURR_POINTS = 20;
POLISH_EVALS    = 15;

%  반구 해 주변 섭동 시드
N_PERTURB    = 8;        % 섭동 시드 개수
PERTURB_FRAC = 0.08;     % 각 변수 범위의 +-8% 이내로 흔든다
N_SEED_RANDOM = 0;       % 무작위 시드 개수. 0 을 기본으로 두는 이유는 이 시험의
                         % 논지가 "반구에서 출발했다" 이기 때문. (surrogateopt 는
                         % 어차피 전역 박스를 스스로 표집하므로, 이 설정은
                         % '반구 근방을 조밀하게 덮은' 탐색이 된다.)

HEMI_FILE = 'opt_hemisphere_result.mat';

%% ===== 개선 판정 문턱 =====
%  [주의] 이 시험은 구조적으로 '개선' 쪽에 기울어 있다. 반구 점이 후보에
%  들어가 있어 결과가 반구보다 의미 있게 낮게 나올 수 없고, 승자는 노이즈가
%  섞인 저정밀 평가 수십 회 중 최댓값으로 뽑힌다. 그래서 판정에서 조심해야 할
%  것은 '개선을 놓치는 것' 이 아니라 '노이즈를 개선으로 착각하는 것' 이다.
%
%  반복이 N_FINAL_REP=3 회뿐이라 표준편차 추정 자체가 부정확하므로, 정규분포의
%  2.0 이 아니라 pooled t 분포를 쓴다. df = 2*(N_FINAL_REP-1) = 4,
%  단측 95% -> t = 2.132. (df 가 커지면 2.0 으로 수렴한다.)
T_CRIT_TABLE = [6.314 2.920 2.353 2.132 2.015 1.943 1.895 1.860];  % df=1..8, 단측 95%
df_pool  = 2*(N_FINAL_REP-1);
T_CRIT   = T_CRIT_TABLE(min(max(df_pool,1), numel(T_CRIT_TABLE)));
fprintf('[판정] 문턱 t = %.3f (pooled df = %d, 단측 95%%)\n', T_CRIT, df_pool);

%% Optimization Variables (13-dim, 다른 캠페인과 동일)
varNames = {'x2','x3','x4','x5','x6', 'y2','y3','y4','y5','y6', 'dETL','dHTL','stretchZ'};
lb = [0, 0, 0, 0, 0, 0,   0,   0,   0,   0,   10, 10, 0.1];
ub = [1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 1.5, 1.5, 150,150, 3];
nvar = numel(lb);

EVAL_LOG = [];

psOpts = optimoptions('patternsearch', ...
    'MaxFunctionEvaluations', POLISH_EVALS, ...
    'InitialMeshSize', 0.1, 'MeshTolerance', 1e-3, ...
    'Cache', 'on', 'Display', 'off');

%% ===== 반구 결과 로드 =====
if ~exist(HEMI_FILE, 'file')
    error(['%s 이 pwd 에 없다. 이 스크립트는 반구 최적해를 시드로 쓰므로 ' ...
           '먼저 opt_hemisphere_arms.m 을 돌리거나 결과 파일을 준비할 것.'], HEMI_FILE);
end
Dh = load(HEMI_FILE);
for f = {'HEMI_X','HEMI_Y','hemi_x','hemi_val','hemi_tot'}
    if ~isfield(Dh, f{1})
        error('%s 에 %s 가 없다. 파일이 예상과 다르다.', HEMI_FILE, f{1});
    end
end
HEMI_X = Dh.HEMI_X(:).';   HEMI_Y = Dh.HEMI_Y(:).';
expand = @(v) [HEMI_X, HEMI_Y, v(1), v(2), v(3)];   % 3 -> 13 (opt_hemisphere_arms.m 과 동일)

fprintf('[Hemisphere] %s 로드. arm %d개\n', HEMI_FILE, size(Dh.hemi_x,1));

%% ===== 결과 컨테이너 =====
base_val  = nan(nBand,1);   base_tot = nan(nBand,1);
base_bins = nan(nBand,4);   base_sd  = nan(nBand,1);
ws_val    = nan(nBand,1);   ws_tot   = nan(nBand,1);
ws_bins   = nan(nBand,4);   ws_sd    = nan(nBand,1);
ws_x      = nan(nBand,nvar);
ws_src    = repmat({''}, nBand, 1);
verdict   = repmat({''}, nBand, 1);

if isempty(ARMS_TO_RUN), ARMS_TO_RUN = 1:nBand; end

%% =====================================================================
%  arm 루프
%% =====================================================================
for k = ARMS_TO_RUN(:).'
    band = BAND_LIST{k};
    bcol = BAND_COL(k);
    if bcol == 0, wtag = -1; else, wtag = band(1); end

    fprintf('\n########## ARM %s  (%d/%d) ##########\n', BAND_NAMES{k}, k, nBand);

    if size(Dh.hemi_x,1) < k || any(~isfinite(Dh.hemi_x(k,:)))
        fprintf('  [Skip] 반구 결과에 이 arm 이 없다.\n');
        continue;
    end

    x_hemi = expand(Dh.hemi_x(k,:));
    if ~isValidPoints(x_hemi)
        fprintf('  [Skip] 반구 해가 현재 제약에서 valid 하지 않다.\n');
        continue;
    end

    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);
    ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    eval_count = 0;
    EVAL_PHASE = 9;  EVAL_W = wtag;

    % ---------- (1) 반구 기준값을 여기서 다시 잰다 ----------
    %  아카이브 값을 그대로 쓰지 않는 이유: 두 스크립트 간 사소한 설정 차이가
    %  '개선폭' 에 섞이면 판정이 무의미해진다. 같은 세션, 같은 정밀도로 잰다.
    ray_nums_current = RAY_FINAL;  wave_n_current = WAVE_N_FINAL;
    et_r = nan(N_FINAL_REP,1);  bins_r = nan(N_FINAL_REP,4);
    for r = 1:N_FINAL_REP
        [et_r(r), bins_r(r,:)] = simulate_bands(x_hemi);
    end
    base_tot(k)    = mean(et_r,'omitnan');
    base_bins(k,:) = mean(bins_r,1,'omitnan');
    if bcol == 0
        base_val(k) = base_tot(k);
        base_sd(k)  = std(et_r,'omitnan');
    else
        base_val(k) = base_bins(k,bcol);
        base_sd(k)  = std(bins_r(:,bcol),'omitnan');
    end

    % 아카이브와 교차 검증 (판정에는 쓰지 않고, 설정 어긋남만 잡아낸다)
    ref = Dh.hemi_val(k);
    if isfinite(ref) && ref > 0
        dev = 100*abs(base_val(k)-ref)/ref;
        if dev > 2
            fprintf(['  [Warn] 반구 재측정 %.5f 이 아카이브 %.5f 와 %.1f%% 어긋난다.\n' ...
                     '         PATCH_XY / 스택 / 파장 설정이 opt_hemisphere_arms.m 과\n' ...
                     '         같은지 확인할 것. 이 상태의 비교는 신뢰할 수 없다.\n'], ...
                     base_val(k), ref, dev);
        else
            fprintf('  반구 재측정 %.5f (아카이브 %.5f, 차이 %.1f%%) OK\n', ...
                    base_val(k), ref, dev);
        end
    end
    fprintf('  기준값 = %.5f +- %.5f (n=%d),  EQE_total=%.5f\n', ...
            base_val(k), base_sd(k), N_FINAL_REP, base_tot(k));

    % 정규화 상수: 기준값 기반이라 arm 마다 자동으로 맞는다.
    refB = max(base_val(k), 1e-6);

    % ---------- (2) 시드 구성: 반구 해 + 근방 섭동 ----------
    ray_nums_current = RAY_SEARCH;  wave_n_current = WAVE_N_SEARCH;
    seedMat = x_hemi;
    nAdded = 0;  nTry = 0;
    while nAdded < N_PERTURB && nTry < 40*N_PERTURB
        nTry = nTry + 1;
        cand = x_hemi + (2*rand(1,nvar)-1) .* (PERTURB_FRAC*(ub-lb));
        cand = min(max(cand, lb), ub);
        cand(1:5) = sort(cand(1:5));            % x 단조성 복구
        if REQUIRE_MONOTONIC_X && any(diff(cand(1:5)) < 0), continue; end
        if ~isValidPoints(cand), continue; end
        seedMat(end+1,:) = cand;  %#ok<SAGROW>
        nAdded = nAdded + 1;
    end
    if N_SEED_RANDOM > 0
        seedMat = [seedMat; genValidPoints(N_SEED_RANDOM, lb, ub)];
    end
    fprintf('  시드 %d개 (반구 1 + 섭동 %d + 무작위 %d)\n', ...
            size(seedMat,1), nAdded, N_SEED_RANDOM);

    % ---------- (3) 세 갈래 탐색 ----------
    cands = {x_hemi};  srcs = {'hemisphere seed'};

    % (i) surrogateopt
    sopts = optimoptions('surrogateopt', ...
        'MaxFunctionEvaluations', EVALS_PER_BAND, ...
        'MinSurrogatePoints',     MIN_SURR_POINTS, ...
        'InitialPoints',          struct('X', seedMat), ...
        'UseParallel', false, 'PlotFcn', [], 'Display', 'iter');
    try
        xS = surrogateopt(@(x) band_objconstr(x, bcol, refB), lb, ub, sopts);
        xS = xS(:).';
        if ~isempty(xS) && isValidPoints(xS)
            cands{end+1} = xS; srcs{end+1} = 'surrogateopt';
        end
    catch ME
        fprintf('  [Warn] surrogateopt 실패: %s\n', ME.message);
        xS = [];
    end

    % (ii) surrogateopt 승자에서 국소 정련
    if ~isempty(xS) && isValidPoints(xS)
        try
            xP = patternsearch(@(x) band_polish(x, bcol, refB), xS, ...
                               [],[],[],[], lb, ub, [], psOpts);
            xP = xP(:).';
            if isValidPoints(xP), cands{end+1} = xP; srcs{end+1} = 'polish(surrogate)'; end
        catch ME
            fprintf('  [Warn] polish(surrogate) 실패: %s\n', ME.message);
        end
    end

    % (iii) **반구 해에서 직접** 국소 정련 — 이 시험의 핵심
    try
        xH = patternsearch(@(x) band_polish(x, bcol, refB), x_hemi, ...
                           [],[],[],[], lb, ub, [], psOpts);
        xH = xH(:).';
        if isValidPoints(xH), cands{end+1} = xH; srcs{end+1} = 'polish(hemisphere)'; end
    catch ME
        fprintf('  [Warn] polish(hemisphere) 실패: %s\n', ME.message);
    end

    % ---------- (4) 저정밀로 후보 선별 ----------
    bestScore = -inf; bestX = []; bestSrc = '';
    for c = 1:numel(cands)
        [et, bins] = simulate_bands(cands{c});
        if ~isfinite(et), continue; end
        if bcol == 0, sc = et; else, sc = bins(bcol); end
        fprintf('    후보 %-20s : %.5f\n', srcs{c}, sc);
        if sc > bestScore, bestScore = sc; bestX = cands{c}; bestSrc = srcs{c}; end
    end
    if isempty(bestX)
        fprintf('  [Warn] 유효 후보 없음\n');
        continue;
    end

    % ---------- (5) 승자 고정밀 재평가 ----------
    ray_nums_current = RAY_FINAL;  wave_n_current = WAVE_N_FINAL;
    et_r = nan(N_FINAL_REP,1);  bins_r = nan(N_FINAL_REP,4);
    for r = 1:N_FINAL_REP
        [et_r(r), bins_r(r,:)] = simulate_bands(bestX);
    end
    ws_x(k,:)    = bestX;
    ws_src{k}    = bestSrc;
    ws_tot(k)    = mean(et_r,'omitnan');
    ws_bins(k,:) = mean(bins_r,1,'omitnan');
    if bcol == 0
        ws_val(k) = ws_tot(k);   ws_sd(k) = std(et_r,'omitnan');
    else
        ws_val(k) = ws_bins(k,bcol);  ws_sd(k) = std(bins_r(:,bcol),'omitnan');
    end

    % ---------- (6) 노이즈 대비 판정 ----------
    %  두 평균의 차이를 두 표준오차의 결합과 비교한다. 표준오차 = sd/sqrt(n).
    se   = sqrt(base_sd(k)^2 + ws_sd(k)^2) / sqrt(N_FINAL_REP);
    dlt  = ws_val(k) - base_val(k);
    rel  = 100*dlt/base_val(k);
    tval = dlt / max(se, eps);
    if tval > T_CRIT
        verdict{k} = sprintf('개선 (%+.2f%%, t=%.1f > %.2f) — 형상 자유도 실재', rel, tval, T_CRIT);
    else
        verdict{k} = sprintf('정체 (%+.2f%%, t=%.1f <= %.2f) — 반구가 국소 최적', rel, tval, T_CRIT);
    end
    fprintf('  >>> 승자 %s : %.5f (기준 %.5f, %+.2f%%, t=%.1f)\n', ...
            bestSrc, ws_val(k), base_val(k), rel, tval);
    fprintf('      %s\n', verdict{k});

    % crash-safe 저장
    save('warmstart_hemisphere_result.mat', 'EVAL_LOG', 'ws_x', 'ws_val', 'ws_tot', ...
         'ws_bins', 'ws_sd', 'ws_src', 'base_val', 'base_tot', 'base_bins', 'base_sd', ...
         'verdict', 'BAND_LIST', 'BAND_NAMES', 'BAND_COL', 'HEMI_X', 'HEMI_Y', ...
         'varNames', 'lb', 'ub', 'GEOM_MISMATCH_LOG');
end

%% =====================================================================
%  요약
%% =====================================================================
fprintf('\n############ WARM START FROM HEMISPHERE ############\n');
fprintf('%-12s | %10s | %10s | %8s | %7s | %s\n', ...
        'arm', 'hemisphere', 'warm start', 'delta%', 't', 'winner');
anyImproved = false;  nRun = 0;  nFlag = 0;
for k = 1:nBand
    if ~isfinite(ws_val(k))
        fprintf('%-12s | %10s | %10s | %8s | %7s | %s\n', BAND_NAMES{k}, '-','-','-','-','(미실행)');
        continue;
    end
    se  = sqrt(base_sd(k)^2 + ws_sd(k)^2) / sqrt(N_FINAL_REP);
    dlt = ws_val(k) - base_val(k);
    nRun = nRun + 1;
    fprintf('%-12s | %10.5f | %10.5f | %+7.2f%% | %6.1f | %s\n', ...
            BAND_NAMES{k}, base_val(k), ws_val(k), 100*dlt/base_val(k), ...
            dlt/max(se,eps), ws_src{k});
    if dlt > T_CRIT*se, anyImproved = true; nFlag = nFlag + 1; end
end

%  [다중비교] arm 을 여러 개 돌리면 그중 하나가 우연히 문턱을 넘을 확률이
%  arm 수만큼 늘어난다 (arm 당 5%, 5개면 약 23%). 하나만 넘었다면 그 arm 을
%  반복 수를 늘려 재확인하기 전에는 '개선' 으로 보고하지 말 것.
if nRun > 1
    fprintf('\n[다중비교] arm %d개 실행 -> 우연히 하나 이상 넘을 확률 약 %.0f%%.\n', ...
            nRun, 100*(1-0.95^nRun));
end

fprintf('\n---- 판정 ----\n');
for k = 1:nBand
    if ~isempty(verdict{k}), fprintf('  %-12s : %s\n', BAND_NAMES{k}, verdict{k}); end
end

if anyImproved
    fprintf(['\n=> arm %d개에서 반구를 문턱 이상으로 넘었다.\n' ...
             '   원래 캠페인의 G_j < 1 은 탐색 예산 부족이었다는 뜻이므로,\n' ...
             '   원고 2.2절의 saturation 서술을 개선폭에 맞춰 수정해야 한다.\n' ...
             '   단, 넘은 arm 이 1개뿐이고 여러 arm 을 돌렸다면 다중비교를\n' ...
             '   의심하고 N_FINAL_REP 를 늘려 그 arm 만 재확인할 것.\n'], nFlag);
else
    fprintf(['\n=> 어느 arm 도 반구를 노이즈 이상으로 넘지 못했다.\n' ...
             '   이 시험은 구조상 개선 쪽에 기울어 있다 — 반구 점이 후보에\n' ...
             '   포함되고, 승자는 노이즈 섞인 저정밀 평가 수십 회 중 최댓값으로\n' ...
             '   뽑히며, 국소 정련을 반구에서 직접 한 번 더 돌린다. 그렇게 하고도\n' ...
             '   개선이 없다는 것이므로, 반구는 13차원 설계공간의 국소 최적점이다.\n' ...
             '   원고 2.2절의 "차원의 비용" 서술을 이 결과로 대체할 수 있다.\n']);
end

save('warmstart_hemisphere_result.mat', 'EVAL_LOG', 'ws_x', 'ws_val', 'ws_tot', ...
     'ws_bins', 'ws_sd', 'ws_src', 'base_val', 'base_tot', 'base_bins', 'base_sd', ...
     'verdict', 'BAND_LIST', 'BAND_NAMES', 'BAND_COL', 'HEMI_X', 'HEMI_Y', ...
     'varNames', 'lb', 'ub', 'GEOM_MISMATCH_LOG');
fprintf('\nsaved -> warmstart_hemisphere_result.mat\n');
report_geom_rejection(GEOM_MISMATCH_LOG, GEOM_TOL);

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

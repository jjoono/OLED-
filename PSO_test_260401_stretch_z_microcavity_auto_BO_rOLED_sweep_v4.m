% ============================================================
%  Bayesian Optimization (bayesopt) + r_pat sweep + warm-start  [v4]
%
%  [v4 변경점: local optimum 탈출력 강화 (계산량 증가 거의 없음)]
%   (A) 최종 고정밀 검증 후보에 "직전 r_pat 최적점(prev_best)" 추가
%       -> 직전 최적 형상을 현재 r_pat에서 직접 재평가. BO가 이번 r_pat에서
%          길을 잃어도 이 후보가 하한선 역할 -> 인위적 dip 원천 억제.
%          (조작이 아님: 실제 시뮬레이션 결과가 더 좋은 design을 채택하는 것)
%   (B) ExplorationRatio 0.5(기본) -> 0.7: warm-start가 exploitation을 이미
%       보장하므로 acquisition은 탐험 쪽으로 균형 이동
%   (C) 시드 구성 재조정: warm-start 상위점 N_WARM 5->3, fresh valid 시드
%       비중 확대 (이전 지역 3 : 새 지역 7) -> GP가 옛 branch에 덜 묶임
%   (D) 단조성 위반 재시도 순서 반전: resume(같은 GP 계속 파기)보다
%       fresh 시드 재초기화(새 지역 탐색)를 먼저 시도 -> local trap 탈출력 ↑
%
%  - PSO 대체: 비싼 블랙박스(레이트레이싱)에 표본 효율적
%  - isValidPoints는 XConstraintFcn으로 연결 -> 불가능 형상은 시뮬 전에 제외
%
%  [v2 변경점: 시간 대비 정확도/재현성 개선]
%   (1) 탐색/검증 ray 수 분리: 탐색은 저정밀(빠름), 최종 검증은 고정밀
%       -> 같은 시간에 더 많은 탐색 평가 + 보고값의 노이즈 대폭 감소
%   (2) 크래시/기하오류 시 0 대신 NaN 반환
%       -> bayesopt가 오류점으로 따로 모델링, GP 스케일 오염 방지
%   (3) warm-start를 "직전 최적 1점" -> "직전 상위 N_WARM점"으로 확장
%   (4) BO 종료 후 patternsearch 국소 정련 (BO=전역 탐색, PS=국소 수렴)
%   (5) 최종 best 후보를 고정밀 ray로 N_FINAL_REP회 반복 평가
%       -> mean±std 기록 (논문 error bar / 재현성 근거)
%   (6) [버그 수정] 코팅(.coa) 파일을 LightTools가 읽기 전에 fclose
%       (기존: 버퍼 미플러시 상태로 읽혀 평가값이 실행마다 흔들릴 수 있음)
%   (7) [v2.1] r_pat마다 수렴 판정 기반 적응 예산:
%       CONV_BLOCK회 추가 평가 후 추정 최적값(estimated)의 상대 개선이
%       CONV_TOL 미만인 상태가 CONV_PATIENCE 블록 연속이면 수렴으로 판정하고
%       다음 r_pat로 진행. 미수렴이면 resume()으로 평가를 계속 추가
%       (안전 상한 MAX_EVAL_PER_RPAT).
%       * rng('shuffle')은 기존 그대로 유지
%
%  - 결과: rpat_trend (best EQE_40_60), rpat_std (반복 평가 표준편차)
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

%% ===== 정확도/시간 트레이드오프 설정 (v2) =====
RAY_SEARCH   = 20000;    % BO/정련 탐색용 ray 수 (노이즈 sigma ~ 1/sqrt(N), 속도 우선)
RAY_FINAL    = 100000;   % 최종 검증용 ray 수 (보고값 정밀도 우선)
N_FINAL_REP  = 3;        % 최종 best 반복 평가 횟수 -> mean±std
N_WARM       = 3;        % (v4-C) 직전 r_pat warm-start 상위점 5->3 (옛 branch 종속성 완화)
N_FRESH      = 7;        % (v4-C) 매 r_pat fresh valid 시드 개수 (이전 3 : 새 7)
POLISH_EVALS = 15;       % BO 후 patternsearch 국소 정련 평가 예산
EXPLORATION_RATIO = 0.7; % (v4-B) EI-plus 탐험 비중 (기본 0.5 -> 0.7)

% --- 수렴 판정 파라미터 (v2.1) ---
INIT_EVAL_FIRST   = 60;    % 첫 r_pat 초기 BO 예산 (시드 20개 포함)
INIT_EVAL_NEXT    = 30;    % 이후 r_pat 초기 BO 예산 (warm-start 포함)
CONV_BLOCK        = 10;    % 수렴 판정 단위: 한 번에 추가하는 평가 횟수
CONV_TOL          = 0.002; % 블록당 추정 최적 EQE의 상대 개선 < 0.2% 면 "개선 없음"
CONV_PATIENCE     = 2;     % 연속 CONV_PATIENCE 블록 개선 없음 -> 수렴 판정
MAX_EVAL_PER_RPAT = 150;   % r_pat당 안전 상한 (수렴 안 해도 여기서 중단)

% --- (v2.2) 단조성 자기진단 + 재시뮬레이션 파라미터 ---
% [주의] 이 로직은 "r_pat_list(s) >= MONOTONIC_FROM_RPAT 구간에서는 물리적으로
%  단조증가해야 한다"는 사전 가정을 전제로 한다. 이 가정 자체가 틀렸을 가능성
%  (예: 실제로는 특정 구간에서 간섭/공진 효과로 비단조적일 수 있음)을 배제할 수
%  없으므로, 재시도 상한(MAX_RETRY_PER_RPAT)을 반드시 두고, 상한 도달 시
%  "강제로 단조증가하게 덮어쓰지 않고" rpat_flagged로 표시만 하여 나중에 사람이
%  직접 판단하게 한다. 이 로직이 결과를 마사지(p-hacking)하지 않도록 주의할 것.
MONOTONIC_FROM_RPAT = 3;     % 이 r_pat 값부터 단조증가 기대 (사용자 요청: 약 2~3부터)
MONOTONIC_TOL_K     = 2;     % 허용 오차 = K * sqrt(std_i^2 + std_{i-1}^2) (노이즈 감안)
MAX_RETRY_PER_RPAT  = 3;     % r_pat 1개당 재시뮬레이션 최대 횟수 (무한루프 방지)
RETRY_EXTRA_EVAL    = 20;    % 재시도 1회당 추가 BO 평가 수 (resume)

%% Optimization Variables (13-dim, 기존 lb/ub와 동일)
varNames = {'x2','x3','x4','x5','x6', 'y2','y3','y4','y5','y6', 'dETL','dHTL','stretchZ'};
lb = [0, 0, 0, 0, 0, 0,   0,   0,   0,   0,   10, 10, 0.1];
ub = [1, 1, 1, 1, 1, 1.5, 1.5, 1.5, 1.5, 1.5, 150,150, 3];

optVars = optimizableVariable.empty(0, numel(lb));
for i = 1:numel(lb)
    optVars(i) = optimizableVariable(varNames{i}, [lb(i), ub(i)]);
end

% [참고] isValidPoints의 feasible 영역이 매우 좁아(무작위 ~0.3%) bayesopt가
%   "XConstraintFcn을 충족한 점이 적다"는 경고를 띄울 수 있다. InitialX로 valid
%   시드를 직접 공급하므로 정상 동작하며, 이 경고는 무시해도 된다.

%% ===== r_pat Sweep Driver =====
r_pat_list   = 25;                   % 1 ~ 25, step 1 -> 1x25
n_pat        = numel(r_pat_list);
rpat_trend   = zeros(1, n_pat);          % 각 r_pat의 best EQE_40_60 (고정밀 재평가 mean)
rpat_std     = zeros(1, n_pat);          % 고정밀 반복 평가 표준편차 (논문 error bar)
rpat_evals   = zeros(1, n_pat);          % r_pat별 실제 사용한 BO 평가 횟수
rpat_conv    = false(1, n_pat);          % r_pat별 수렴 판정 통과 여부
rpat_retry   = zeros(1, n_pat);          % (v2.2) 단조성 위반으로 재시뮬레이션한 횟수
rpat_flagged = false(1, n_pat);          % (v2.2) 재시도 상한까지도 위반 지속 -> 수동 확인 필요
rpat_results = cell(1, n_pat);           % bayesopt 결과 객체 보관
prev_best    = [];                       % warm-start용 직전 최적해 (table row)
prev_results = [];                       % warm-start용 직전 bayesopt 결과

for s = 1:n_pat
    r_pat = r_pat_list(s);               % global -> objFcn이 이 값 사용
    eval_count = 0;                      % r_pat마다 재시작 카운터 리셋
    % (v4-A) 이 시점의 prev_best = "직전 r_pat의 최적점". 아래에서 prev_best가
    % 현재 r_pat 결과로 덮어써지기 전에 앵커로 따로 보관 (재시도 루프에서도 사용)
    if s > 1
        xPrevRpat = table2array(prev_best);
    else
        xPrevRpat = [];
    end
    fprintf('\n############ r_pat %d/%d : r_pat = %.2f ############\n', s, n_pat, r_pat);

    % r_pat 시작 시 LightTools 클린 상태 보장
    RenewLightTools();
    lt = ltloc.GetLTAPI(ID_swept);
    ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

    if s == 1
        % 첫 r_pat: valid 시드 20개 + 초기 예산 INIT_EVAL_FIRST
        maxEval  = INIT_EVAL_FIRST;
        initX    = array2table(genValidPoints(20, lb, ub), 'VariableNames', varNames);
    else
        % (v4-C) 이후 r_pat: warm-start 상위점은 3개로 줄이고 fresh valid 시드
        % 비중을 늘림 (이전 지역 3 : 새 지역 7). warm-start가 GP를 옛 branch에
        % 과하게 묶어 local trap을 만드는 것을 완화.
        maxEval  = INIT_EVAL_NEXT;
        warmTbl  = [prev_best; topKPoints(prev_results, N_WARM)];
        warmTbl  = unique(warmTbl, 'stable');
        warmTbl  = warmTbl(1:min(height(warmTbl), N_WARM), :);   % 상한 N_WARM개로 제한
        freshTbl = array2table(genValidPoints(N_FRESH, lb, ub), 'VariableNames', varNames);
        initX    = [warmTbl; freshTbl];
    end

    % --- (1) 전역 탐색: BO는 저정밀(빠른) ray로 ---
    ray_nums_current = RAY_SEARCH;

    results = bayesopt(@bo_objective, optVars, ...
        'MaxObjectiveEvaluations', maxEval, ...
        'XConstraintFcn',          @bo_xconstraint, ...
        'IsObjectiveDeterministic', false, ...           % 레이트레이싱 노이즈 모델링
        'AcquisitionFunctionName', 'expected-improvement-plus', ...
        'ExplorationRatio',        EXPLORATION_RATIO, ...  % (v4-B) 탐험 비중 강화
        'Verbose', 1, ...
        'PlotFcn', {}, ...                               % 내부 플롯 끔(속도)
        'InitialX', initX);

    % --- (v2.1) 수렴 판정 루프: 개선이 멈출 때까지 CONV_BLOCK씩 추가 평가 ---
    % 판정 기준: 추정 최적 EQE(-MinEstimatedObjective)의 블록당 상대 개선이
    % CONV_TOL 미만인 상태가 CONV_PATIENCE 블록 연속이면 수렴.
    % (관측 최소가 아니라 GP 추정 최소를 쓰므로 MC 노이즈에 덜 흔들린다)
    noImpCount = 0;
    while noImpCount < CONV_PATIENCE && results.NumObjectiveEvaluations < MAX_EVAL_PER_RPAT
        prevBestEst = -results.MinEstimatedObjective;
        addEval = min(CONV_BLOCK, MAX_EVAL_PER_RPAT - results.NumObjectiveEvaluations);
        results = resume(results, 'MaxObjectiveEvaluations', addEval);  % 이어서 평가 추가
        newBestEst = -results.MinEstimatedObjective;
        relImp = (newBestEst - prevBestEst) / max(abs(newBestEst), eps);
        if relImp < CONV_TOL
            noImpCount = noImpCount + 1;   % 이번 블록에서 유의미한 개선 없음
        else
            noImpCount = 0;                % 개선됨 -> 카운터 리셋, 계속 탐색
        end
        fprintf('[Converge] r_pat=%.2f | evals=%3d | bestEst EQE=%.5f | relImp=%+.4f | noImp %d/%d\n', ...
            r_pat, results.NumObjectiveEvaluations, newBestEst, relImp, noImpCount, CONV_PATIENCE);
    end
    rpat_evals(s) = results.NumObjectiveEvaluations;
    rpat_conv(s)  = (noImpCount >= CONV_PATIENCE);
    if rpat_conv(s)
        fprintf('[Converge] r_pat=%.2f : %d회 평가에서 수렴 판정 통과.\n', r_pat, rpat_evals(s));
    else
        fprintf('[Converge] r_pat=%.2f : 상한 %d회 도달로 종료 (수렴 기준 미충족 - 결과 확인 필요).\n', ...
            r_pat, MAX_EVAL_PER_RPAT);
    end

    % 노이즈가 있으므로 추정 최소(estimated)를 사용
    bestX_BO = results.XAtMinEstimatedObjective;
    x0 = table2array(bestX_BO);

    % --- (2) 국소 정련: BO 추정 최적점에서 patternsearch로 미세 수렴 ---
    % BO는 예산이 작을 때 최적점 "근처"까지만 가는 경우가 많다. 미분 불필요한
    % patternsearch로 소량(POLISH_EVALS) 추가 평가하여 국소 수렴을 마무리.
    fprintf('--- Local polish (patternsearch, %d evals) ---\n', POLISH_EVALS);
    psOpts = optimoptions('patternsearch', ...
        'MaxFunctionEvaluations', POLISH_EVALS, ...
        'InitialMeshSize', 0.1, ...
        'MeshTolerance', 1e-3, ...
        'Cache', 'on', ...
        'Display', 'iter');
    try
        xPol = patternsearch(@polish_objective, x0, [],[],[],[], lb, ub, [], psOpts);
    catch perr
        fprintf('[Warn] patternsearch 실패(%s). BO 결과만 사용합니다.\n', perr.message);
        xPol = x0;
    end

    % --- (3) 최종 검증: 후보들을 고정밀 ray로 반복 평가 후 승자 채택 ---
    % (v4-A) 직전 r_pat의 최적점(prev_best)도 후보로 포함: 직전 최적 형상을
    % 현재 r_pat에서 직접 재평가한 값이 사실상의 하한선 역할을 하므로,
    % 이번 r_pat의 BO가 나쁜 local optimum에 갇혀도 dip이 억제된다.
    if isequal(xPol, x0)
        candX = {x0};
    else
        candX = {x0, xPol};
    end
    if ~isempty(xPrevRpat) && ~any(cellfun(@(cc) isequal(cc, xPrevRpat), candX))
        candX{end+1} = xPrevRpat;
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
        fprintf('  candidate %d : EQE_40_60 = %.5f ± %.5f (N=%d, %d rays)\n', ...
            c, candMean(c), candStd(c), N_FINAL_REP, RAY_FINAL);
    end
    [bestEQE, ci] = max(candMean);
    if ~isfinite(bestEQE)
        % 고정밀 검증이 전부 실패(크래시 등)한 경우 BO 추정값으로 폴백
        bestEQE = -results.MinEstimatedObjective;
        ci = 1;
    end
    bestX = array2table(candX{ci}, 'VariableNames', varNames);

    rpat_trend(s)   = bestEQE;
    rpat_std(s)     = candStd(ci);
    rpat_results{s} = results;
    prev_best       = bestX;             % 다음 r_pat warm-start
    prev_results    = results;

    % ===== (v2.2) 단조성 자기진단 + 재시뮬레이션 피드백 루프 =====
    % r_pat_list(s) >= MONOTONIC_FROM_RPAT 구간에서 rpat_trend(s)가 직전보다
    % (노이즈 감안 허용오차 이상으로) 떨어지면, 같은 r_pat에 대해 추가 예산을
    % 투입해 재탐색 -> 재정련 -> 재검증을 반복한다. 상한 도달 시 강제로 값을
    % 뜯어고치지 않고 rpat_flagged로만 표시한다(아래 주의사항 참고).
    if s > 1 && r_pat_list(s) >= MONOTONIC_FROM_RPAT
        while true
            noiseTol  = MONOTONIC_TOL_K * sqrt(rpat_std(s)^2 + rpat_std(s-1)^2);
            violation = (rpat_trend(s) < rpat_trend(s-1) - noiseTol);
            if ~violation
                break;
            end
            if rpat_retry(s) >= MAX_RETRY_PER_RPAT
                rpat_flagged(s) = true;
                fprintf(['[Monotonic] r_pat=%.2f : %d회 재시도 후에도 단조성 위반 지속 ' ...
                    '(EQE=%.5f < r_pat=%.2f의 %.5f, tol=%.5f). 강제 수정하지 않고 플래그만 ' ...
                    '표시 -> 실제 물리적 비단조성일 수 있으니 수동 확인 필요.\n'], ...
                    r_pat, rpat_retry(s), rpat_trend(s), r_pat_list(s-1), rpat_trend(s-1), noiseTol);
                break;
            end

            rpat_retry(s) = rpat_retry(s) + 1;
            fprintf(['[Monotonic] r_pat=%.2f : EQE=%.5f < 이전 r_pat=%.2f의 %.5f - tol(%.5f). ' ...
                '재시뮬레이션 %d/%d회차 시작...\n'], ...
                r_pat, rpat_trend(s), r_pat_list(s-1), rpat_trend(s-1), noiseTol, ...
                rpat_retry(s), MAX_RETRY_PER_RPAT);

            % 1) (v4-D) 재탐색: fresh 시드 재초기화를 우선 시도.
            %    같은 GP를 resume으로 계속 파면 동일 local optimum만 더 정밀하게
            %    파고들 가능성이 높다. local trap 탈출이 목적이므로 첫 재시도부터
            %    "완전히 새 랜덤 시드 + 직전 r_pat 앵커"로 GP를 새로 만든다.
            %    (마지막 재시도 1회는 기존 GP resume으로 마무리 정련 기회 부여)
            ray_nums_current = RAY_SEARCH;
            if rpat_retry(s) < MAX_RETRY_PER_RPAT
                % fresh 시드 20개 + 직전 r_pat 앵커 -> 새 지역 탐색 유도
                freshMat = genValidPoints(20, lb, ub);
                if ~isempty(xPrevRpat)
                    freshMat = [xPrevRpat; freshMat];
                end
                initX2 = array2table(freshMat, 'VariableNames', varNames);
                results = bayesopt(@bo_objective, optVars, ...
                    'MaxObjectiveEvaluations', INIT_EVAL_NEXT, ...
                    'XConstraintFcn', @bo_xconstraint, ...
                    'IsObjectiveDeterministic', false, ...
                    'AcquisitionFunctionName', 'expected-improvement-plus', ...
                    'ExplorationRatio', EXPLORATION_RATIO, ...
                    'Verbose', 1, 'PlotFcn', {}, 'InitialX', initX2);
                rpat_evals(s) = rpat_evals(s) + results.NumObjectiveEvaluations; % 새 GP: 전량 가산
            else
                % 마지막 재시도: 기존 GP를 이어서 추가 정련
                addEval = min(RETRY_EXTRA_EVAL, MAX_EVAL_PER_RPAT - results.NumObjectiveEvaluations);
                if addEval > 0
                    results = resume(results, 'MaxObjectiveEvaluations', addEval);
                    rpat_evals(s) = rpat_evals(s) + addEval;                     % 추가분만 가산
                end
            end

            % 2) 국소 정련 재수행
            bestX_BO = results.XAtMinEstimatedObjective;
            x0 = table2array(bestX_BO);
            try
                xPol = patternsearch(@polish_objective, x0, [],[],[],[], lb, ub, [], psOpts);
            catch perr
                fprintf('[Warn] patternsearch 실패(%s). BO 결과만 사용합니다.\n', perr.message);
                xPol = x0;
            end

            % 3) 고정밀 재검증 ((v4-A) 직전 r_pat 앵커 포함)
            if isequal(xPol, x0), candX = {x0}; else, candX = {x0, xPol}; end
            if ~isempty(xPrevRpat) && ~any(cellfun(@(cc) isequal(cc, xPrevRpat), candX))
                candX{end+1} = xPrevRpat;
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
            end
            [newEQE, ci] = max(candMean);

            fprintf('  [Monotonic] 재시도 %d 결과: EQE=%.5f (이전 시도 EQE=%.5f)\n', ...
                rpat_retry(s), newEQE, rpat_trend(s));

            % 재시도 결과가 기존보다 나은 경우에만 채택 (더 나쁘면 기존 최선값 유지)
            if isfinite(newEQE) && newEQE > rpat_trend(s)
                rpat_trend(s)   = newEQE;
                rpat_std(s)     = candStd(ci);
                bestX           = array2table(candX{ci}, 'VariableNames', varNames);
                prev_best       = bestX;
                prev_results    = results;
                rpat_results{s} = results;
            end
        end
    end

    save('BO_rpat_trend.mat', 'rpat_trend', 'rpat_std', 'rpat_evals', 'rpat_conv', ...
        'rpat_retry', 'rpat_flagged', 'r_pat_list', 'rpat_results');
    fprintf(['############ r_pat = %.2f Done : Best EQE_40_60 = %.5f ± %.5f ' ...
        '(evals=%d, converged=%d, retries=%d, flagged=%d) ############\n'], ...
        r_pat, bestEQE, rpat_std(s), rpat_evals(s), rpat_conv(s), rpat_retry(s), rpat_flagged(s));
end

disp('=== All r_pat sweeps finished ===');
disp('rpat_trend  ='); disp(rpat_trend);
disp('rpat_std    ='); disp(rpat_std);
disp('rpat_evals  ='); disp(rpat_evals);
disp('rpat_conv   ='); disp(rpat_conv);
disp('rpat_retry  ='); disp(rpat_retry);
disp('rpat_flagged='); disp(rpat_flagged);
if any(rpat_flagged)
    fprintf(['\n[주의] 아래 r_pat 지점들은 재시도 상한(%d회)까지 시도해도 단조성 위반이 ' ...
        '해소되지 않았습니다. 코드가 강제로 값을 조작하지 않았으니, 이 지점들은 ' ...
        '실제로 비단조적인 물리 현상일 가능성을 배제하지 말고 직접 확인하세요:\n'], MAX_RETRY_PER_RPAT);
    disp(r_pat_list(rpat_flagged));
end

figure(2);
errorbar(r_pat_list, rpat_trend, rpat_std, '-o', 'LineWidth', 2);
hold on;
if any(rpat_flagged)
    plot(r_pat_list(rpat_flagged), rpat_trend(rpat_flagged), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
end
hold off;
xlabel('r\_pat'); ylabel('Best EQE\_40\_60');
title('r\_pat Sweep Trend (BO + local polish, high-ray verified)'); grid on;


%% ===== LightTools 1회 평가 공용 래퍼 (v2) =====
% 주기적 재시작 + 크래시 처리. 크래시/기하오류는 NaN 반환:
%  - bayesopt는 NaN을 오류점으로 따로 모델링 -> GP 스케일이 0으로 오염되지 않음
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

%% ===== bayesopt 목적함수 래퍼 (EQE 최대화 -> -EQE 최소화) =====
function obj = bo_objective(Xtbl)
obj = -simulate_EQE(table2array(Xtbl));   % NaN이면 bayesopt가 오류점으로 처리
end

%% ===== patternsearch 정련용 목적함수 =====
function f = polish_objective(x)
if ~isValidPoints(x)
    f = 0;          % 무효 형상: 시뮬 없이 벌점 (valid면 -EQE < 0 이므로 항상 열등)
    return;
end
e = simulate_EQE(x);
if isnan(e), e = 0; end
f = -e;
end

%% ===== bayesopt 제약함수: 기하학적 valid 형상만 통과 =====
function tf = bo_xconstraint(Xtbl)
% bayesopt는 true=feasible(평가 진행)로 해석. isValidPoints도 true=valid.
pts = table2array(Xtbl);   % N x 13 (열 순서 = varNames)
tf  = isValidPoints(pts);  % N x 1 logical
end

%% ===== 직전 bayesopt 결과에서 상위 K개 관측점 추출 (warm-start용) =====
function T = topKPoints(results, K)
X  = results.XTrace;
f  = results.ObjectiveTrace;
ok = isfinite(f);          % 오류(NaN) 평가 제외
X  = X(ok, :);
f  = f(ok);
[~, idx] = sort(f, 'ascend');   % -EQE 오름차순 = EQE 내림차순
K  = min(K, numel(idx));
T  = X(idx(1:K), :);
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


%% Objective Function (v2: ray 수 가변 + .coa fclose 수정 + rng 비오염 파일명)
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

% (v2) 탐색/검증 단계에 따라 ray 수 가변
if isempty(ray_nums_current)
    ray_nums = 50000;
else
    ray_nums = ray_nums_current;
end

List=ltml.LTDbList(lt,'lens_manager[1]','SIMULATIONS');
Key=ltml.LTListByName(lt,List,'ForwardAll');
ltml.LTDbSet(lt,Key,'MaxProgress',ray_nums);
% (v2, 선택) Monte-Carlo 시드 고정(common random numbers): 같은 형상 재평가 시
%   같은 값이 나와 GP가 보는 노이즈가 줄어든다. LightTools GUI의
%   Simulation Input 탭 > Random Number Seed의 데이터액세스 이름을 확인한 뒤
%   아래 주석을 해제해서 사용 (이름이 버전에 따라 다를 수 있음):
% try, ltml.LTDbSet(lt,Key,'RandomNumberSeed',12345); catch, end
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
fclose(fileID);  % (v2 버그수정) LightTools가 읽기 전에 버퍼 플러시 + 파일 잠금 해제

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

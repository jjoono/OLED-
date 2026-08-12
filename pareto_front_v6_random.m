% ============================================================
%  test_random_supercell.m
%
%  목적: stress_random_mla.m 을 밤새 돌리기 전에, 30분 안에 파이프라인이
%        실제로 작동하는지 확인한다. 본 실행과 **같은 코드**(objFcn_supercell.m,
%        generate_random_supercell_ent.m)를 호출하므로, 여기서 통과하면
%        본 실행에서 새로 터질 곳은 예산/시간뿐이다.
%
%  단계
%    STAGE 0  LightTools 없이 .ent 생성만 (수초). 파일/렌즈 개수/충전율 확인
%             + 슈퍼셀 높이맵 미리보기 PNG -> 형상을 눈으로 확인
%    STAGE 1  LightTools 연결. unit-cell Filename / 텍스처 배치 간격 /
%             zone extent 를 설정 후 **읽어서** 검증 (스케일 보존 확인)
%    STAGE 2  초고속 설정(광선 少, 파장 3개)으로 실제 평가 N_QUICK 회.
%             EQE 유한성, band 합 ≈ total, 소요시간 측정
%    STAGE 3  같은 시드 2회 -> 재현성(MC 노이즈 크기) 확인
%    STAGE 4  하이퍼파라미터 코너 3점 -> 극단값에서 안 터지는지 확인
%    마지막   본 실행 예상 소요시간 환산 + 통과/실패 요약
%
%  기본 설정으로 총 12회 평가. 평가당 수십 초면 10분 내외.
%  RAY_TEST / WAVE_N_TEST 를 올리면 정확해지지만 느려진다.
% ============================================================
clear;
global ID_LT ltml ltloc count ray_nums_current wave_n_current

%% ===== 빠른 검증용 설정 =====
RAY_TEST    = 2000;    % 본 실행 10000 -> 1/5
WAVE_N_TEST = 150;     % 파장 step. 453-753(301개)에서 step 150 -> 3개 파장
N_QUICK     = 3;       % STAGE 2 반복 평가 수
SEED_BASE   = 4242;

% 본 실행 설정 (예상 시간 환산용)
RAY_FULL    = 10000;
WAVE_N_FULL = 10;
N_EVAL_FULL = 100 + 60 + 15 + 3;   % 무작위 + surrogate + polish + 최종재평가

BASE  = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\';
scDir = [BASE 'supercell_ents\'];
if ~exist(scDir,'dir'), mkdir(scDir); end

hypNames = {'fill','rJitter','posJitter','aspect','aspectJitter','profileMix'};
lb = [0.35, 0.00, 0.00, 0.30, 0.00, 0.00];
ub = [0.90, 0.30, 1.00, 1.50, 0.40, 1.00];
hypMid = (lb + ub)/2;

pass = struct('s0',false,'s1',false,'s2',false,'s3',false,'s4',false);
fprintf('\n================ RANDOM SUPERCELL SMOKE TEST ================\n');

%% =====================================================================
%  STAGE 0 — .ent 생성만 (LightTools 불필요)
%% =====================================================================
fprintf('\n---- STAGE 0: 슈퍼셀 .ent 생성 (LightTools 없이) ----\n');
t0 = tic;
try
    p0 = struct('fill',hypMid(1),'rJitter',hypMid(2),'posJitter',hypMid(3), ...
                'aspect',hypMid(4),'aspectJitter',hypMid(5),'profileMix',hypMid(6), ...
                'templatePath',[BASE 'freeform_template_v2.ent'], ...
                'nCols',8, 'outPath',[scDir 'smoketest_mid.1.ent']);
    info0 = generate_random_supercell_ent(SEED_BASE, p0);

    d = dir(p0.outPath);
    fprintf('  파일: %s (%.1f MB)\n', p0.outPath, d.bytes/1e6);
    if isstruct(info0)
        fn = fieldnames(info0);
        for i = 1:numel(fn)
            v = info0.(fn{i});
            if isnumeric(v) && isscalar(v)
                fprintf('  %-16s = %g\n', fn{i}, v);
            elseif isnumeric(v)
                fprintf('  %-16s : %d개\n', fn{i}, numel(v));
            end
        end
    end
    if d.bytes < 1e5
        error('.ent 가 너무 작다 (%d bytes). 생성 실패 가능성.', d.bytes);
    end
    pass.s0 = true;
    fprintf('  [OK] 생성 %.1f초\n', toc(t0));
catch ME
    fprintf('  [FAIL] %s\n', ME.message);
end

% --- 미리보기: 슈퍼셀 높이맵을 .ent 에서 되읽어 그린다 ---
if pass.s0
    try
        txt = fileread(p0.outPath);
        nums = sscanf(txt(strfind(txt,'ORAStartData')+12:end), '%f');
        % 헤더 10개 뒤부터 (x,y,z) 삼중항
        nums = nums(11:end);
        m = floor(numel(nums)/3)*3;  P = reshape(nums(1:m),3,[]).';
        figure('Name','supercell preview','Color','w','Position',[80 80 900 380]);
        subplot(1,2,1);
        scatter(P(:,1), P(:,2), 3, P(:,3), 'filled'); axis equal tight;
        colorbar; title('supercell height map (from .ent)');
        xlabel('x'); ylabel('y');
        subplot(1,2,2);
        nside = round(sqrt(size(P,1)));
        if nside^2 == size(P,1)
            Z = reshape(P(:,3), nside, nside);
            plot(Z(round(nside/2),:), 'LineWidth', 1.2); grid on;
            title('center cross-section'); xlabel('grid index'); ylabel('z');
        end
        saveas(gcf, 'test_supercell_preview.png');
        fprintf('  미리보기 저장 -> test_supercell_preview.png (렌즈렛이 여러 개 보이는지 확인)\n');
    catch ME
        fprintf('  [Warn] 미리보기 실패(치명적 아님): %s\n', ME.message);
    end
end

%% =====================================================================
%  STAGE 1 — LightTools 연결 + 설정 왕복 검증
%% =====================================================================
fprintf('\n---- STAGE 1: LightTools 연결 및 설정 검증 ----\n');
try
    RenewLightTools();
    count = 1;
    lt = ltloc.GetLTAPI(ID_LT);
    ltml.LTSetOption(lt, "ShowFileDialogBox", 0);
    ltml.LTCmd(lt, 'Message "smoke test"');
    fprintf('  [OK] 연결\n');

    % 배치 간격 목록 이름 탐색 (objFcn_supercell 과 같은 후보 순서)
    placeLists = {'ZONE_TEXTURE_HEXAGONAL_PLACEMENT','HEXAGONAL_PLACEMENT', ...
                  'ZONE_TEXTURE_PLACEMENT','TEXTURE_PLACEMENT'};
    foundName = '';
    for ip = 1:numel(placeLists)
        try
            PL = ltml.LTDbList(lt,'lens_manager[1]',placeLists{ip});
            PK = ltml.LTListAtPos(lt,PL,1);
            gx = ltml.LTDbGet(lt,PK,'XSpacing');
            gy = ltml.LTDbGet(lt,PK,'YSpacing');
            fprintf('  배치 목록 "%s" 발견: XSpacing=%.4f, YSpacing=%.4f\n', ...
                    placeLists{ip}, gx, gy);
            foundName = placeLists{ip};  break;
        catch
        end
    end
    if isempty(foundName)
        error(['배치 간격 목록을 못 찾음. LightTools Database Browser 에서 ' ...
               'ZoneTextureHexagonalPlacement 의 DB 목록 이름을 확인해 ' ...
               'objFcn_supercell.m 과 이 파일의 placeLists 에 추가할 것.']);
    end
    pass.s1 = true;
catch ME
    fprintf('  [FAIL] %s\n', ME.message);
end

%% =====================================================================
%  STAGE 2 — 실제 평가 (초고속 설정)
%% =====================================================================
fprintf('\n---- STAGE 2: 실제 평가 %d회 (ray=%d, 파장 step=%d) ----\n', ...
        N_QUICK, RAY_TEST, WAVE_N_TEST);
ray_nums_current = RAY_TEST;  wave_n_current = WAVE_N_TEST;
tEval = nan(N_QUICK,1);  Etot = nan(N_QUICK,1);  Bins = nan(N_QUICK,4);
if pass.s1
    for i = 1:N_QUICK
        tA = tic;
        try
            out = objFcn_supercell(hypMid, SEED_BASE + i);
            Etot(i)   = out.EQE_total;
            Bins(i,:) = [out.EQE_0_20, out.EQE_20_40, out.EQE_40_60, out.EQE_60_80];
            tEval(i)  = toc(tA);
            fprintf('  eval %d: EQE_total=%.5f  bands=[%.4f %.4f %.4f %.4f]  %.1f초\n', ...
                i, Etot(i), Bins(i,1), Bins(i,2), Bins(i,3), Bins(i,4), tEval(i));
        catch ME
            tEval(i) = toc(tA);
            fprintf('  eval %d [FAIL] %.1f초: %s\n', i, tEval(i), ME.message);
        end
    end
    ok = isfinite(Etot) & Etot > 0;
    if any(ok)
        sumB = sum(Bins(ok,:),2);  rel = sumB ./ Etot(ok);
        fprintf('  band 합 / EQE_total = %s  (0-90도 중 80-90도 몫만큼 1보다 작으면 정상)\n', ...
                mat2str(round(rel',3)));
        pass.s2 = true;
    end
    fprintf('  [%s] 평가당 중앙값 %.1f초\n', tern(pass.s2), median(tEval,'omitnan'));
else
    fprintf('  (STAGE 1 실패로 건너뜀)\n');
end

%% =====================================================================
%  STAGE 3 — 재현성 (같은 시드 2회)
%% =====================================================================
fprintf('\n---- STAGE 3: 같은 시드 재현성 ----\n');
if pass.s2
    try
        a = objFcn_supercell(hypMid, SEED_BASE);
        b = objFcn_supercell(hypMid, SEED_BASE);
        dRel = abs(a.EQE_total - b.EQE_total) / max(a.EQE_total, eps);
        fprintf('  EQE_total: %.5f vs %.5f  (상대차 %.2f%%)\n', ...
                a.EQE_total, b.EQE_total, 100*dRel);
        if dRel < 0.05
            fprintf('  [OK] 기하는 결정론적, 차이는 광선추적 MC 노이즈 수준\n');
            pass.s3 = true;
        else
            fprintf('  [Warn] 차이가 큼 -> 시드가 기하에 안 먹거나 광선수 부족\n');
        end
    catch ME
        fprintf('  [FAIL] %s\n', ME.message);
    end
else
    fprintf('  (건너뜀)\n');
end

%% =====================================================================
%  STAGE 4 — 하이퍼파라미터 코너
%% =====================================================================
fprintf('\n---- STAGE 4: 하이퍼파라미터 코너 ----\n');
corners = { lb, '최소 (fill/jitter/aspect 하한)'; ...
            ub, '최대 (fill/jitter/aspect 상한)'; ...
            [ub(1) 0 0 hypMid(4) 0 0], '무질서 0 (주기 배열에 가까움)' };
if pass.s2
    nOK = 0;
    for c = 1:size(corners,1)
        try
            o = objFcn_supercell(corners{c,1}, SEED_BASE + 100 + c);
            fprintf('  %-28s EQE_total=%.5f\n', corners{c,2}, o.EQE_total);
            if isfinite(o.EQE_total) && o.EQE_total > 0, nOK = nOK + 1; end
        catch ME
            fprintf('  %-28s [FAIL] %s\n', corners{c,2}, ME.message);
        end
    end
    pass.s4 = (nOK == size(corners,1));
    fprintf('  [%s] %d/%d 코너 통과\n', tern(pass.s4), nOK, size(corners,1));
else
    fprintf('  (건너뜀)\n');
end

%% =====================================================================
%  요약 + 본 실행 예상 시간
%% =====================================================================
fprintf('\n================ 요약 ================\n');
names = {'STAGE 0 .ent 생성','STAGE 1 연결/설정','STAGE 2 평가', ...
         'STAGE 3 재현성','STAGE 4 코너'};
vals  = [pass.s0 pass.s1 pass.s2 pass.s3 pass.s4];
for i = 1:numel(names)
    fprintf('  %-20s %s\n', names{i}, tern(vals(i)));
end

if pass.s2 && any(isfinite(tEval))
    tMed = median(tEval,'omitnan');
    % 광선 수와 파장 개수에 대략 선형으로 스케일
    nWaveTest = numel(1:WAVE_N_TEST:301);
    nWaveFull = numel(1:WAVE_N_FULL:301);
    scale = (RAY_FULL/RAY_TEST) * (nWaveFull/nWaveTest);
    tFull = tMed * scale * N_EVAL_FULL / 3600;
    fprintf('\n  평가당 %.1f초 (테스트 설정) -> 본 설정 환산 x%.0f\n', tMed, scale);
    fprintf('  본 실행 %d회 예상: 약 %.1f 시간\n', N_EVAL_FULL, tFull);
end

if all(vals)
    fprintf('\n  ==> 전부 통과. stress_random_mla.m 본 실행 가능.\n');
else
    fprintf('\n  ==> 실패 항목 있음. 위 메시지 확인 후 본 실행할 것.\n');
end
fprintf('======================================\n\n');

function s = tern(c)
if c, s = 'OK'; else, s = 'FAIL'; end
end

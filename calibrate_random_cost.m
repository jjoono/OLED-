% ============================================================
%  calibrate_random_cost.m
%
%  목적: random supercell 계열의 평가 비용을 어디까지 줄여도 결론이 안 바뀌는지
%        **측정해서** 정한다. (추정으로 두 번 틀렸으므로 재는 편이 빠르다)
%
%  배경: 실측 분해 결과 평가 시간의 대부분이 광선추적이고, 그 비용은
%        (렌즈렛 수) x (메쉬 점수) x (광선 수) x (파장 수) 에 붙는다.
%        슈퍼셀은 렌즈렛 64개 + 201x201 메쉬라 주기 배열보다 훨씬 비싸다.
%
%  방법: 기준 설정 1회(비쌈)를 재고, 값싼 변형들을 각각 1회씩 재서
%        - 시간이 얼마나 줄었는지
%        - **선택성**(= band/total, 논문이 실제로 쓰는 양)이 얼마나 흔들리는지
%        를 비교한다. 총 EQE 절대값보다 선택성 보존이 판정 기준이다.
%
%  소요: 기준 1회 약 15분 + 변형 5회 약 15분 = 30분 내외.
%
%  출력: 표 + calibrate_random_cost.mat + 권장 설정 출력
% ============================================================
clear;
global ID_LT ltml ltloc count ray_nums_current wave_n_current ...
       CACHE_FIXED_STACK RND_NCOLS RND_NGRID

CACHE_FIXED_STACK = true;
SEED = 7777;
hyp  = [0.625, 0.15, 0.5, 0.90, 0.20, 0.50];   % 하이퍼파라미터 중앙값

% 설정: {이름, nCols, nGrid, rays, waveStep}
CFG = { ...
  'reference (현재)',        8, 201, 10000, 10; ...
  'rays 1/2',                8, 201,  5000, 10; ...
  'wave step 20',            8, 201, 10000, 20; ...
  'grid 141',                8, 141, 10000, 10; ...
  'nCols 6',                 6, 201, 10000, 10; ...
  'combo (5000/20/141/6)',   6, 141,  5000, 20 };

nC = size(CFG,1);
tSec = nan(nC,1);  Etot = nan(nC,1);  Sel = nan(nC,4);

RenewLightTools();
count = 1;
fprintf('\n============ RANDOM SUPERCELL COST CALIBRATION ============\n');
fprintf('기준 설정부터 잽니다. 첫 회에는 CPS/코팅 캐시 준비 비용이 포함됩니다.\n');

% 캐시 워밍업 (측정에 캐시 준비 비용이 섞이지 않도록 값싼 설정으로 1회)
RND_NCOLS = 6;  RND_NGRID = 121;
ray_nums_current = 2000;  wave_n_current = 150;
try
    objFcn_supercell(hyp, SEED);
    fprintf('캐시 워밍업 완료.\n\n');
catch ME
    fprintf('[Warn] 워밍업 실패(계속 진행): %s\n\n', ME.message);
end

for k = 1:nC
    RND_NCOLS = CFG{k,2};  RND_NGRID = CFG{k,3};
    ray_nums_current = CFG{k,4};  wave_n_current = CFG{k,5};
    nWave = numel(1:CFG{k,5}:301);
    fprintf('[%d/%d] %-24s nCols=%d grid=%d rays=%d 파장%d개 ... ', ...
            k, nC, CFG{k,1}, CFG{k,2}, CFG{k,3}, CFG{k,4}, nWave);
    tA = tic;
    try
        o = objFcn_supercell(hyp, SEED);       % 같은 시드 = 같은 기하
        tSec(k) = toc(tA);
        Etot(k) = o.EQE_total;
        Sel(k,:) = [o.EQE_0_20, o.EQE_20_40, o.EQE_40_60, o.EQE_60_80] / o.EQE_total;
        fprintf('%.0f초  EQE=%.4f\n', tSec(k), Etot(k));
    catch ME
        tSec(k) = toc(tA);
        fprintf('FAIL (%.0f초): %s\n', tSec(k), ME.message);
    end
end

%% ===== 비교표 =====
fprintf('\n%-24s %8s %8s %8s %10s %10s\n', ...
        '설정','시간(s)','속도배','EQE','ΔEQE(%)','Δ선택성(%p)');
ref = 1;
for k = 1:nC
    if ~isfinite(tSec(k)), continue; end
    spd  = tSec(ref)/tSec(k);
    dE   = 100*(Etot(k)-Etot(ref))/Etot(ref);
    dSel = 100*max(abs(Sel(k,:)-Sel(ref,:)));      % 최대 선택성 편차 (퍼센트포인트)
    fprintf('%-24s %8.0f %8.2f %8.4f %10.1f %10.2f\n', ...
            CFG{k,1}, tSec(k), spd, Etot(k), dE, dSel);
end

fprintf(['\n기준 선택성: %.3f %.3f %.3f %.3f\n'], Sel(ref,:));
fprintf(['판정 기준: 선택성 편차가 0.5%%p 이내면 결론(계열 간 비교)에 영향 없음.\n' ...
         '           총 EQE 는 몇 %% 흔들려도 무방 — 논문이 쓰는 값은 비율이다.\n']);

%% ===== 권장 설정 =====
okIdx = find(isfinite(tSec) & (100*max(abs(Sel-Sel(ref,:)),[],2))' <= 0.5);
if ~isempty(okIdx)
    [~, jj] = max(tSec(ref)./tSec(okIdx));
    best = okIdx(jj);
    N_SEARCH = 50 + 60 + 15;      % 무작위 + surrogate + polish (축소안)
    tEstim = tSec(best) * N_SEARCH / 3600;
    fprintf('\n>>> 권장: %s  (%.1f배 빠름)\n', CFG{best,1}, tSec(ref)/tSec(best));
    fprintf('    stress_random_mla.m 에 넣을 값:\n');
    fprintf('      RND_NCOLS = %d;  RND_NGRID = %d;\n', CFG{best,2}, CFG{best,3});
    fprintf('      RAY_SEARCH = %d;  WAVE_N_SEARCH = %d;\n', CFG{best,4}, CFG{best,5});
    fprintf('      N_RANDOM = 50;   (100 -> 50 축소)\n');
    fprintf('    탐색 %d회 예상: 약 %.1f 시간 (+ 고정밀 최종 재평가 별도)\n', ...
            N_SEARCH, tEstim);
else
    fprintf('\n>>> 선택성 0.5%%p 이내를 만족하는 축소 설정 없음. 기준 설정 유지 필요.\n');
end

save('calibrate_random_cost.mat','CFG','tSec','Etot','Sel','hyp','SEED');
fprintf('\nsaved -> calibrate_random_cost.mat\n');
fprintf('==========================================================\n\n');

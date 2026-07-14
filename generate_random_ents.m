function generate_random_ents(N, outDir, templatePath, ZLO, ZHI, seed)
% GENERATE_RANDOM_ENTS  LightTools import 견고성 점검용 랜덤 freeform .ent 배치 생성.
%
%   generate_random_ents(N, outDir, templatePath, ZLO, ZHI, seed)
%
%   BO 가 실제로 탐색하는 방식(내부 격자점 Z 를 [ZLO,ZHI]에서 무작위, 경계는 0)대로
%   N 개의 .ent 를 만든다. 각 파일을 LightTools 에서 "라이브러리 구성요소 로드"로
%   하나씩 불러 "잘못된 솔리드" 없이 정상 solid 로 열리는지 확인하면, BO 중 매 평가의
%   지오메트리 주입이 안전한지 사전 검증된다.
%
%   입력(모두 선택, 기본값 있음)
%     N            : 생성 개수 (기본 20)
%     outDir       : 출력 폴더 (기본 '.\ent_tests')
%     templatePath : 검증된 템플릿 .ent (기본 '.\freeform_template.ent')
%     ZLO, ZHI     : 내부점 Z 범위 (기본 0.0, 0.8) — 두께 1mm 보다 작게
%     seed         : 난수 시드 (기본 42; 재현성)
%
%   예: generate_random_ents(30)   % 현재 폴더 ent_tests\ 에 30개
%
%   [주의] freeform_grid_info.m, generate_freeform_ent.m 가 같은 경로에 있어야 함.

    if nargin < 1 || isempty(N),            N = 20;                      end
    if nargin < 2 || isempty(outDir),       outDir = fullfile(pwd,'ent_tests'); end
    if nargin < 3 || isempty(templatePath), templatePath = fullfile(pwd,'freeform_template.ent'); end
    if nargin < 4 || isempty(ZLO),          ZLO = 0.0;                   end
    if nargin < 5 || isempty(ZHI),          ZHI = 0.8;                   end
    if nargin < 6 || isempty(seed),         seed = 42;                   end

    if ~exist(outDir,'dir'), mkdir(outDir); end
    [~, innerIdx, Npts] = freeform_grid_info(templatePath);
    nInner = numel(innerIdx);
    rng(seed);

    fprintf('템플릿 격자 %d점, 내부(자유) %d점. Z∈[%.2f,%.2f], %d개 생성 -> %s\n', ...
        Npts, nInner, ZLO, ZHI, N, outDir);

    for k = 1:N
        zin = ZLO + rand(nInner,1) .* (ZHI - ZLO);
        zFull = zeros(Npts,1);
        zFull(innerIdx) = zin;
        outPath = fullfile(outDir, sprintf('ff_rand_%03d.1.ent', k));
        generate_freeform_ent(zFull, templatePath, outPath);
    end
    fprintf('완료. %s 안의 ff_rand_*.1.ent 를 LightTools 에 하나씩 로드해 유효 solid 확인.\n', outDir);
    fprintf('(전부 정상이면 BO 매 평가의 지오메트리 주입이 안전하다는 뜻.)\n');
end

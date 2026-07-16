%% generate_random_freeform_ents.m
% 무작위 비대칭 3D 자유형 렌즈 .ent 를 여러 개 생성.
% (매끈한 비대칭 돔 + 선택적 날카로운 절단면들 -> 절단 있는 것/없는 것/여러 개 섞임)
%
% [핵심] 여기 build_freeform_height() 가 "고정 DOF 벡터 -> 높이장" 매핑이며,
% 나중에 BO 목적함수도 이 함수를 그대로 호출한다. 즉 이 파일 = 형상 파라미터화의
% 단일 소스. 무작위 생성기는 이 벡터를 bounds 안에서 샘플링할 뿐이다.
%
% [DOF] = 1(h) + 2*M(방위각 harmonic) + 3*K(절단면).  기본 M=4,K=3 -> DOF = 18.
%   * MATLAB bayesopt(GP)의 실용 상한 ~20 을 고려한 값. 더 키우려면 M,K 를 올리되
%     DOF>~25 부터는 GP-BO 가 약해지므로 optimizer 교체(CMA-ES/surrogateopt) 권장.
%   * RAM/스레드는 DOF 상한을 못 늘림 - 병렬 평가(여러 LightTools 인스턴스)와 큰
%     초기샘플에만 도움. (BO 단계에서 bayesopt 'UseParallel',true 활용)
%
% [검증된 규칙만 조합]
%   - 사각(모서리 포함) 격자, Ra=1.2139  (원판/극좌표 금지: 자체교차)
%   - 외피는 taper-rim 돔(테두리 특이점 제거 -> 물결 최소)
%   - 날카로운 엣지는 min(높이, 절단평면)  (특이점 없어 해상도만 있으면 날카로움 유지)
%   - RearSurface 평면 그대로, SmoothResample "No", tbase 넉넉(시뮬에서 파묻힘)

%% ===== 설정 =====
templatePath = 'freeform_template.ent';
outDir = 'random_freeform_ents';
Nlens  = 12;     % 생성 개수
seed   = 1;

cfg.M = 4;       % 방위각 harmonic 차수 (m=1..M)
cfg.K = 3;       % 최대 절단면 개수
cfg.Ra   = 1.2139;   % 사각 격자 반폭 (템플릿과 동일, 모서리 포함)
cfg.Rap  = 1.0;      % 광학 조리개 반경
cfg.n    = 81;       % 렌더링 해상도 (DOF 아님)
cfg.taper= 0.03;     % 외피 테두리 완만화 비율
cfg.tbase= 0.30;     % 밑받침 두께 (고정)

[lb, ub] = freeform_bounds(cfg);
DOF = numel(lb);
fprintf('DOF = %d  (h 1 + harmonic %d + 절단면 %d x3)\n', DOF, 2*cfg.M, cfg.K);

if ~exist(outDir,'dir'), mkdir(outDir); end
rng(seed);

% 격자 미리 계산 (모든 렌즈 공통)
g = linspace(-cfg.Ra, cfg.Ra, cfg.n);
[X, Y] = meshgrid(g, g);

for idx = 1:Nlens
    p = lb + rand(1,DOF).*(ub-lb);

    % 다양성: 절단면 일부를 무작위로 비활성(z0 를 돔 위로 올려 no-op) -> "절단 없는" 형상도 생성
    nActive = randi([0, cfg.K]);              % 0~K 개만 실제 절단
    for k = (nActive+1):cfg.K
        p(1 + 2*cfg.M + 3*(k-1) + 1) = 2.0;   % 이 절단면 z0 를 아주 높게 -> min 에서 무효
    end

    H = build_freeform_height(p, X, Y, cfg);
    outPath = fullfile(outDir, sprintf('ff_rand_%03d.1.ent', idx));
    write_freeform_ent(H, X, Y, cfg.n, cfg.tbase, templatePath, outPath);
    fprintf('[%2d] %s  (active cuts=%d, h=%.2f)\n', idx, outPath, nActive, p(1));
end
fprintf('완료. %s 안의 .ent 를 LightTools 로 확인.\n', outDir);


%% ================= 형상 파라미터화 (BO 공용) =================
function H = build_freeform_height(p, X, Y, cfg)
% p : 1 x DOF 파라미터 벡터. 레이아웃:
%   p(1)                 = h (돔 높이)
%   p(2 : 1+2M)          = [a1 b1 a2 b2 ... aM bM] (방위각 harmonic)
%   p(2+2M : end)        = K x [z0 m phi0] (절단면)
    M = cfg.M; K = cfg.K; Ra = cfg.Ra; Rap = cfg.Rap; taper = cfg.taper;
    r   = hypot(X, Y);
    phi = atan2(Y, X);

    h = p(1);

    % --- 외피: taper-rim 반구 반경프로파일 (h=1 기준 후 h 곱) ---
    Pb = sqrt(max(1-(r/Rap).^2, 0));
    r0 = Rap*(1-taper);
    sag0   = sqrt(max(1-(r0/Rap)^2, 0));
    slope0 = -(r0/Rap)/sqrt(max(1-(r0/Rap)^2, 1e-9));
    Plin = max(sag0 + slope0*(r-r0), 0);
    P = Pb;  P(r>=r0) = Plin(r>=r0);

    % --- 방위각 비대칭 harmonic ---
    S = ones(size(r));
    for m = 1:M
        a = p(1 + 2*(m-1) + 1);
        b = p(1 + 2*(m-1) + 2);
        S = S + (a*cos(m*phi) + b*sin(m*phi)) .* (r/Rap).^m;
    end
    S = max(S, 0.05);

    H = h .* P .* S;

    % --- 절단면들: min() 으로 날카로운 엣지 ---
    base = 1 + 2*M;
    for k = 1:K
        z0   = p(base + 3*(k-1) + 1);
        mm   = p(base + 3*(k-1) + 2);
        phi0 = p(base + 3*(k-1) + 3);
        plane = z0 + mm*(X*cos(phi0) + Y*sin(phi0));
        H = min(H, plane);
    end
    H = max(H, 0);
end

%% ================= 파라미터 bounds (BO 공용) =================
function [lb, ub] = freeform_bounds(cfg)
    M = cfg.M; K = cfg.K;
    lb = zeros(1, 1+2*M+3*K);
    ub = zeros(1, 1+2*M+3*K);
    % h
    lb(1) = 0.5;  ub(1) = 1.1;
    % harmonic a_m,b_m  (고차일수록 진폭 축소 -> 고주파 폭주 방지)
    for m = 1:M
        amp = 0.6 / m;
        lb(1+2*(m-1)+1) = -amp;  ub(1+2*(m-1)+1) = amp;
        lb(1+2*(m-1)+2) = -amp;  ub(1+2*(m-1)+2) = amp;
    end
    % 절단면 z0, m, phi0
    base = 1 + 2*M;
    for k = 1:K
        lb(base+3*(k-1)+1) = 0.20;  ub(base+3*(k-1)+1) = 1.40;  % z0 (>~h 이면 no-op)
        lb(base+3*(k-1)+2) = 0.00;  ub(base+3*(k-1)+2) = 1.60;  % m (기울기, 0=수평평면컷)
        lb(base+3*(k-1)+3) = 0.00;  ub(base+3*(k-1)+3) = 2*pi;  % phi0
    end
end

%% ================= .ent writer =================
function write_freeform_ent(H, X, Y, n, tbase, templatePath, outPath)
    Z = -H;   % 음수 = 볼록(양각)
    Xv = X(:); Yv = Y(:); Zv = Z(:);
    N = n*n;
    tpl = fileread(templatePath);
    tok = regexp(tpl, 'ORAStartData;([\s\S]*?)ORAEndData;', 'tokenExtents');
    s0 = tok{1}(1); e0 = tok{1}(2);
    buf = sprintf('0 1 %d %d 0 0 %d 0 0 0', n, n, N);
    for i = 1:N
        buf = [buf sprintf(' %.17g %.17g %.17g', Xv(i), Yv(i), Zv(i))]; %#ok<AGROW>
    end
    buf = [buf ' 0 0 4 CartesianMapper 1 0 0 0 0'];
    newtxt = [tpl(1:s0-1) char(10) buf char(10) tpl(e0+1:end)];
    newtxt = strrep(newtxt, 'setPosition:  { 0. 0. 1.  } ;', ...
        sprintf('setPosition:  { 0. 0. %g  } ;', tbase));
    newtxt = regexprep(newtxt, 'restoreSmoothResample: "Yes"', 'restoreSmoothResample: "No"', 'once');
    fid = fopen(outPath, 'w');
    fwrite(fid, newtxt);
    fclose(fid);
end

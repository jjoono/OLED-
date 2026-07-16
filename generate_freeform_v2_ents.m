%% generate_freeform_v2_ents.m
% 자유도당 표현력을 최대화한 height-field 파라미터화 (RBF 국소요철 + 절단엣지).
% harmonic(전역 방위각) 대신 RBF bump(임의 위치 국소 자유조각)를 써서 같은 DOF로
% 훨씬 다양한 비대칭 형상을 만든다. import 안전 규칙(사각격자/taper-rim/min-cut/
% SmoothResample off/넉넉 tbase)은 그대로.
%
% [DOF] = 1(h) + 4*Nrbf + 3*Ncut.  기본 Nrbf=3, Ncut=2 -> 19.  (GP-bayesopt ~20 상한)
%
% [자유도 더 키우려면 - 중요]
%   파라미터화가 아니라 OPTIMIZER 를 바꿔야 한다. GP-bayesopt 는 ~20 에서 막힘.
%   Nrbf/Ncut 를 크게(예: 40~80 DOF) 하려면 CMA-ES / surrogateopt / TuRBO 로 교체하고
%   bayesopt 대신 그걸 쓴다. 이때 512GB/64thread 는 '병렬 평가 + 큰 집단'으로 활용
%   (여러 LightTools 인스턴스 동시). build_freeform_v2()/freeform_v2_bounds() 는 어느
%   optimizer 든 그대로 재사용 가능.
%
% [render 해상도 n 은 DOF 아님] - n 키우면 매끈해지지만 최적화 변수는 안 늘어남.

%% ===== 설정 =====
templatePath = 'freeform_template_v2.ent';  % 원점(0,0,0)+단위행렬(축정렬) 로 정리된 템플릿
outDir = 'freeform_v2_ents';
Nlens  = 12;  seed = 1;

cfg.Nrbf = 3;        % RBF 국소요철 개수 (개당 4 DOF)
cfg.Ncut = 2;        % 절단평면 개수 (개당 3 DOF)
cfg.Ra    = 1.2139;  cfg.Rap = 1.0;  cfg.n = 81;
cfg.taper = 0.03;    cfg.tbase = 0.30;

[lb, ub] = freeform_v2_bounds(cfg);
DOF = numel(lb);
fprintf('DOF = %d  (h 1 + RBF %d*4 + cut %d*3)\n', DOF, cfg.Nrbf, cfg.Ncut);

if ~exist(outDir,'dir'), mkdir(outDir); end
rng(seed);
g = linspace(-cfg.Ra, cfg.Ra, cfg.n);
[X, Y] = meshgrid(g, g);

for idx = 1:Nlens
    p = lb + rand(1,DOF).*(ub-lb);
    % 다양성: 절단면 일부 무작위 비활성 (z0 를 위로) -> 절단 있는/없는 형상 섞임
    nActive = randi([0, cfg.Ncut]);
    for k = (nActive+1):cfg.Ncut
        p(1 + 4*cfg.Nrbf + 3*(k-1) + 1) = 2.0;
    end
    H = build_freeform_v2(p, X, Y, cfg);
    outPath = fullfile(outDir, sprintf('ff_v2_%03d.1.ent', idx));
    write_freeform_ent(H, X, Y, cfg.n, cfg.tbase, templatePath, outPath);
    fprintf('[%2d] %s  (active cuts=%d, h=%.2f)\n', idx, outPath, nActive, p(1));
end
fprintf('완료. %s 확인.\n', outDir);


%% ============ 형상 파라미터화 (BO/CMA-ES 공용) ============
function H = build_freeform_v2(p, X, Y, cfg)
% p 레이아웃: p(1)=h | 그다음 Nrbf x [x0 y0 amp sigma] | 그다음 Ncut x [z0 m phi0]
    Nrbf=cfg.Nrbf; Ncut=cfg.Ncut; Ra=cfg.Ra; Rap=cfg.Rap; taper=cfg.taper;
    r = hypot(X, Y);
    h = p(1);

    % 외피: taper-rim 돔
    Pb = sqrt(max(1-(r/Rap).^2, 0));
    r0 = Rap*(1-taper);
    sag0 = sqrt(max(1-(r0/Rap)^2,0));  slope0 = -(r0/Rap)/sqrt(max(1-(r0/Rap)^2,1e-9));
    Plin = max(sag0 + slope0*(r-r0), 0);
    P = Pb;  P(r>=r0) = Plin(r>=r0);
    H = h .* P;

    % RBF 국소 요철 (임의 위치 볼록/오목)
    for i = 1:Nrbf
        b = 1 + 4*(i-1);
        x0 = p(b+1); y0 = p(b+2); amp = p(b+3); sig = p(b+4);
        H = H + amp .* exp(-((X-x0).^2 + (Y-y0).^2) ./ (2*sig^2));
    end

    % rim 창: 테두리에서 부드럽게 0 (RBF 누설 차단, 깨끗한 원형 rim)
    W = ones(size(r));
    rw = Rap*0.9;
    trans = (r-rw)/(Rap-rw);
    zoneT = (r>=rw) & (r<Rap);
    W(zoneT) = 1 - (3*trans(zoneT).^2 - 2*trans(zoneT).^3);  % smootherstep 1->0
    W(r>=Rap) = 0;
    H = H .* W;

    % 절단평면 -> 날카로운 엣지
    base = 1 + 4*Nrbf;
    for k = 1:Ncut
        z0 = p(base+3*(k-1)+1); m = p(base+3*(k-1)+2); phi0 = p(base+3*(k-1)+3);
        H = min(H, z0 + m*(X*cos(phi0) + Y*sin(phi0)));
    end
    H = max(H, 0);
end

%% ============ bounds (BO/CMA-ES 공용) ============
function [lb, ub] = freeform_v2_bounds(cfg)
    Nrbf=cfg.Nrbf; Ncut=cfg.Ncut; Rap=cfg.Rap;
    DOF = 1 + 4*Nrbf + 3*Ncut;
    lb = zeros(1,DOF); ub = zeros(1,DOF);
    lb(1)=0.5; ub(1)=1.1;                         % h
    for i=1:Nrbf
        b=1+4*(i-1);
        lb(b+1)=-0.8*Rap; ub(b+1)=0.8*Rap;        % x0
        lb(b+2)=-0.8*Rap; ub(b+2)=0.8*Rap;        % y0
        lb(b+3)=-0.45;    ub(b+3)=0.55;           % amp (오목/볼록)
        lb(b+4)=0.15;     ub(b+4)=0.55;           % sigma (폭)
    end
    base=1+4*Nrbf;
    for k=1:Ncut
        lb(base+3*(k-1)+1)=0.2; ub(base+3*(k-1)+1)=1.4;  % z0
        lb(base+3*(k-1)+2)=0.0; ub(base+3*(k-1)+2)=1.6;  % m
        lb(base+3*(k-1)+3)=0.0; ub(base+3*(k-1)+3)=2*pi; % phi0
    end
end

%% ============ .ent writer ============
function write_freeform_ent(H, X, Y, n, tbase, templatePath, outPath)
    % [방향] +Z 방향으로 볼록 돌출, 플랫 베이스는 아래(-z). z축으로 뒤집은 것:
    %   FrontSurface Z = +H (위로 볼록),  RearSurface(평면)를 z=-tbase 로 이동.
    Z = H;  Xv=X(:); Yv=Y(:); Zv=Z(:);  N=n*n;
    tpl = fileread(templatePath);
    tok = regexp(tpl, 'ORAStartData;([\s\S]*?)ORAEndData;', 'tokenExtents');
    s0 = tok{1}(1); e0 = tok{1}(2);
    buf = sprintf('0 1 %d %d 0 0 %d 0 0 0', n, n, N);
    for i = 1:N
        buf = [buf sprintf(' %.17g %.17g %.17g', Xv(i), Yv(i), Zv(i))]; %#ok<AGROW>
    end
    buf = [buf ' 0 0 4 CartesianMapper 1 0 0 0 0'];
    newtxt = [tpl(1:s0-1) char(10) buf char(10) tpl(e0+1:end)];
    % RearSurface(평면) z 위치 = -tbase  (SecondLensSurfacePrimitive 의 첫 setPosition,
    % 현재 값이 무엇이든 견고하게 치환). FreeformEntity/LensElementPrimitive 는 템플릿에서
    % 이미 원점(0,0,0)+단위행렬 이므로 손대지 않는다.
    newtxt = regexprep(newtxt, ...
        '(CSGLensSurfacePrimitive_1[\s\S]*?setPosition:  \{ 0\. 0\. )[-0-9.eE]+(  \} ;)', ...
        ['$1' num2str(-tbase,'%g') '$2'], 'once');
    newtxt = regexprep(newtxt, 'restoreSmoothResample: "Yes"', 'restoreSmoothResample: "No"', 'once');
    fid = fopen(outPath, 'w');  fwrite(fid, newtxt);  fclose(fid);
end

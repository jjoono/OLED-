%% generate_asym_sharp_ents.m
% 비대칭 + 날카로운 엣지를 갖는 자유형 3D 렌즈 .ent 를 여러 개 생성.
%
% [원리 - 지금까지 검증된 것만 조합]
%  1) 외피(envelope) = taper-rim 돔: 테두리 마지막 taper 비율만 접선으로 완만화해
%     특이점(무한 기울기)을 제거 -> 격자 해상도만 있으면 물결 없이 매끈 (확인됨).
%  2) 비대칭 날카로운 엣지 = min(돔, 절단평면) 을 K 개 연쇄 적용.
%     각 절단면은 (z0,m,phi0) 3개 변수. 절단면끼리 만나는 크리즈는 특이점이 아니라
%     유한 기울기라 해상도만 있으면 바로 날카롭게 재현됨 (확인됨).
%  3) 데이터는 반드시 "모서리 있는 사각(Rectangular) 격자" (템플릿과 동일 반폭
%     Ra=1.2139) 여야 자체교차가 안 남 (확인됨). 원판/극좌표 격자는 쓰지 않는다.
%  4) tbase(밑받침) 는 고정·넉넉(0.3)하게 둔다 - 시뮬레이션에서 기판에 파묻히는
%     값이라 최적화 변수(DOF)가 아니다.
%
% [DOF] = 1(h) + 3*K(절단면) . K=2 -> DOF=7,  K=3 -> DOF=10  (20 DOF 이내 여유)
%
% 사용법: 아래 설정만 바꿔서 그대로 실행. Nlens 개의 .ent 가 저장된다.

%% ===== 설정 =====
templatePath = 'freeform_template.ent';
outDir       = 'asym_sharp_ents';
Nlens        = 6;      % 생성 개수
K            = 2;      % 절단면 개수 (DOF = 1 + 3*K)
seed         = 1;      % 재현성

Ra    = 1.2139;  % 사각 격자 반폭 (템플릿과 동일, 모서리 포함 필수)
Rap   = 1.0;     % 광학 조리개 반경 (원형 트리밍은 소프트웨어가 별도 처리)
n     = 81;      % 렌더링 해상도 (DOF 아님, 물결/날카로움 세밀도만 결정)
taper = 0.03;    % 외피 테두리 완만화 비율 (특이점 제거용, 고정)
tbase = 0.30;    % 밑받침 두께 (고정·넉넉, 시뮬에서 파묻힘 가정)

hRange   = [0.6, 1.0];     % 돔 높이 범위
z0Range  = [0.2, 0.7];     % 절단면 중심높이 범위
mRange   = [0.5, 1.5];     % 절단면 기울기 범위 (클수록 급한 크리즈)
phiRange = [0, 2*pi];      % 절단면 방향 범위

if ~exist(outDir, 'dir'), mkdir(outDir); end
rng(seed);

g = linspace(-Ra, Ra, n);
[X, Y] = meshgrid(g, g);
r = hypot(X, Y);

% 외피(taper-rim 돔)의 기울기 특이점 제거용 상수 (h 에 비례하므로 h 곱해서 사용)
r0 = Rap*(1-taper);
sag0 = sqrt(max(1-(r0/Rap)^2,0));                       % h=1 기준 r0 에서의 높이
slope0 = -(r0/Rap)/sqrt(max(1-(r0/Rap)^2,1e-9));        % h=1 기준 r0 에서의 기울기

fprintf('DOF = %d (h 1개 + 절단면 %d개 x 3)\n', 1+3*K, K);

for idx = 1:Nlens
    h = hRange(1) + rand*(diff(hRange));

    % 1) 외피: taper-rim 돔
    dome = h * sqrt(max(1-(r/Rap).^2,0));
    lin  = max(h*sag0 + h*slope0*(r-r0), 0);
    H = dome; H(r>=r0) = lin(r>=r0);

    % 2) K 개의 절단평면을 순차 적용 -> 비대칭 + 날카로운 크리즈들
    params = zeros(K,3);
    for k = 1:K
        z0   = z0Range(1)   + rand*diff(z0Range);
        m    = mRange(1)    + rand*diff(mRange);
        phi0 = phiRange(1)  + rand*diff(phiRange);
        params(k,:) = [z0, m, phi0];
        plane = z0 + m*(X*cos(phi0) + Y*sin(phi0));
        H = min(H, plane);
    end
    H = max(H, 0);

    outPath = fullfile(outDir, sprintf('ff_asym_sharp_%03d.1.ent', idx));
    write_freeform_ent(H, X, Y, n, tbase, templatePath, outPath);
    fprintf('[%d] %s : h=%.2f | 절단면(z0,m,phi0)=\n', idx, outPath, h);
    disp(params);
end

fprintf('완료. %s 안의 .ent 를 LightTools 에 하나씩 로드해 확인.\n', outDir);


function write_freeform_ent(H, X, Y, n, tbase, templatePath, outPath)
    Z = -H;   % 음수 = 볼록(양각) 규약
    Xv = X(:); Yv = Y(:); Zv = Z(:);
    N = n*n;
    tpl = fileread(templatePath);
    tok = regexp(tpl, 'ORAStartData;([\s\S]*?)ORAEndData;', 'tokenExtents');
    s0 = tok{1}(1); e0 = tok{1}(2);   % FrontSurface 블록만 교체 (RearSurface는 템플릿 그대로 = 평면)
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

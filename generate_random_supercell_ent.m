function [outPath, info] = generate_random_supercell_ent(seed, params)
% GENERATE_RANDOM_SUPERCELL_ENT  비주기(무작위 조립) MLA 의사난수 슈퍼셀 .ent 생성.
%
%   [outPath, info] = generate_random_supercell_ent(seed, params)
%
% [목적] 15x15 mm 어레이를 ~10 um 렌즐릿 수백만 개로 직접 채우는 것은 불가능하므로,
%   표준 pseudo-random supercell 트릭을 쓴다: 하나의 큰 unit cell(= 기존 텍스처
%   패치와 같은 크기) 안에서 렌즐릿의 위치/크기/프로파일을 무작위화하고, 기존
%   텍스처 타일링(x_pattern x y_pattern)이 그 슈퍼셀을 반복하게 한다.
%   상관길이 << 슈퍼셀 크기 이면 통계적으로 무작위 어레이와 동등하다.
%
% [포맷] generate_freeform_v2_ents.m 의 write_freeform_ent() 를 그대로 따른다:
%   - freeform_template_v2.ent (원점정렬 템플릿) 의 첫 ORAStartData 블록만 교체
%   - 사각(모서리 포함) 격자 반폭 Ra = 1.2139 (기존 freeform 렌즈들과 동일 footprint
%     -> 텍스처 unit-cell 크기/타일링을 바꾸지 않고 기존 렌즈 자리에 바로 교체 가능)
%   - Z = +H (위로 볼록), RearSurface 평면을 z = -tbase 로 이동
%   - restoreSmoothResample: "No"
%
% [렌즐릿 프로파일 클래스] pareto_front_freeform.m 의 13-var 클래스 중 형상 부분
%   (10-var: 제어점 x2..x6, y2..y6, 고정 끝점 (0,1),(1,0))과 동일한 규칙을 쓴다:
%   - 동일한 isValidPoints 유효성 규칙 (자기교차/공선/간격/꺾임각) — 아래에 그대로 복사
%   - REQUIRE_MONOTONIC_X 와 동일하게 x2..x6 오름차순 정렬
%   각 렌즐릿은 자기만의 유효 제어점을 뽑고, 반경방향 프로파일 z(rho) 를
%   pchip 보간으로 평가해 회전 대칭(swept)으로 세운다. (swept 모델과 같은 클래스)
%
% [배치] jittered-hexagonal: 육각 격자(피치 p) + 균일 jitter. jitter 는
%   (p - 2*r_max)/2 로 캡 -> 이웃 간 비중첩이 기하학적으로 보장된다.
%
% 입력
%   seed   : 난수 시드. 같은 (seed, params) -> 비트단위 동일한 .ent (결정론적).
%   params : struct, 모든 필드 선택적 (기본값 아래). 하이퍼파라미터 6개 +
%            해상도/경로 설정.
%     --- 무작위 조립 하이퍼파라미터 (stress_random_mla.m 의 최적화 변수) ---
%     .fill         목표 fill factor (렌즐릿 밑면적 / 셀면적)     [기본 0.65]
%     .rJitter      반경 jitter 비율 (r_i = r_mean*(1±rJitter))   [기본 0.15]
%     .posJitter    위치 jitter 강도, 허용 최대치의 비율 0..1      [기본 0.5]
%     .aspect       평균 종횡비 h_i/r_i (StretchZ 역할)           [기본 1.0]
%     .aspectJitter 종횡비 jitter 비율                            [기본 0.15]
%     .profileMix   0 = 전부 반구 프로파일, 1 = 전부 무작위 유효 제어점,
%                   중간값 = 반구와 무작위의 볼록결합               [기본 0.5]
%     --- 구조/해상도 (최적화 변수 아님) ---
%     .Ra     슈퍼셀 반폭 [정규 단위, 템플릿과 동일]               [기본 1.2139]
%     .nCols  육각 격자 열 수 (피치 p = 2*Ra/nCols 결정)          [기본 8]
%     .n      렌더링 격자 해상도 (n x n)                          [기본 201]
%     .tbase  밑받침 두께                                         [기본 0.30]
%     .templatePath                                               [기본 'freeform_template_v2.ent']
%     .outPath                                                    [기본 'supercell_<seed>.1.ent']
%
% 출력
%   outPath : 쓰인 .ent 경로
%   info    : struct — nLens, pitch, rMean, rList, centers, fillRealized, H (미리보기용)
%
% [주의] 렌즐릿이 슈퍼셀 경계를 넘지 않도록 중심을 안쪽으로 제한한다(경계에서
%   잘린 렌즐릿은 타일 이음매 불연속을 만든다). 그 대가로 실현 fill 이 목표보다
%   약간 낮다 (info.fillRealized 로 확인).

%% ===== 기본값 =====
if nargin < 2, params = struct(); end
def.fill = 0.65;  def.rJitter = 0.15;  def.posJitter = 0.5;
def.aspect = 1.0; def.aspectJitter = 0.15; def.profileMix = 0.5;
def.Ra = 1.2139;  def.nCols = 8;  def.n = 201;  def.tbase = 0.30;
def.templatePath = 'freeform_template_v2.ent';
def.outPath = '';
fn = fieldnames(def);
for k = 1:numel(fn)
    if ~isfield(params, fn{k}) || isempty(params.(fn{k}))
        params.(fn{k}) = def.(fn{k});
    end
end
if isempty(params.outPath)
    params.outPath = sprintf('supercell_%d.1.ent', seed);
end
outPath = params.outPath;

rng(seed, 'twister');    % 결정론: 같은 seed -> 같은 슈퍼셀

Ra = params.Ra;  n = params.n;

%% ===== 1) jittered-hexagonal 배치 =====
p = 2*Ra / params.nCols;                       % 육각 피치
% 목표 fill 로부터 평균 반경: 육각 셀 면적 = sqrt(3)/2 * p^2
r_mean = p * sqrt(params.fill * sqrt(3) / (2*pi));
r_mean = min(r_mean, 0.48*p);                  % 비중첩 절대 상한 (기저 격자 기준)
r_max_draw = r_mean * (1 + params.rJitter);
r_max_draw = min(r_max_draw, 0.495*p);

% 육각 격자 중심 (짝수/홀수 행 offset)
dy = p * sqrt(3)/2;
ys = 0:dy:(Ra - r_max_draw);   ys = unique([-ys, ys]);
centers = [];
for iy = 1:numel(ys)
    rowIdx = round(ys(iy)/dy);
    xoff = 0.5*p * mod(rowIdx, 2);
    xs = (-Ra):p:(Ra);  xs = xs + xoff;
    for ix = 1:numel(xs)
        if abs(xs(ix)) <= Ra - r_max_draw && abs(ys(iy)) <= Ra - r_max_draw
            centers(end+1,:) = [xs(ix), ys(iy)]; %#ok<AGROW>
        end
    end
end
nLens = size(centers,1);
if nLens == 0
    error('generate_random_supercell_ent:layout', ...
        '렌즐릿이 하나도 안 들어감: nCols=%d, fill=%.2f 확인.', params.nCols, params.fill);
end

% 반경 추첨 (균일 ±rJitter)
rList = r_mean * (1 + params.rJitter * (2*rand(nLens,1) - 1));

% 위치 jitter: 이웃 비중첩 보장 캡 = (p - 2*max r)/2, 여기에 posJitter 비율 적용
jmax = max(0, (p - 2*max(rList))/2) * params.posJitter;
theta = 2*pi*rand(nLens,1);  rad = jmax * sqrt(rand(nLens,1));   % 원판 내 균일
centers = centers + [rad.*cos(theta), rad.*sin(theta)];
% 경계 재클램프 (jitter 후에도 셀 내부 유지)
for i = 1:nLens
    lim = Ra - rList(i);
    centers(i,:) = max(min(centers(i,:), lim), -lim);
end

% 종횡비 추첨
aspList = params.aspect * (1 + params.aspectJitter * (2*rand(nLens,1) - 1));

%% ===== 2) 렌즐릿별 프로파일 (13-var 클래스의 10-var 형상부와 동일 규칙) =====
% 반구 canonical 제어점: 사분원 위의 등각 5점 (x=sin, y=cos)
th_c = (1:5) * (pi/2) / 6;
hemiX = sin(th_c);  hemiY = cos(th_c);

profiles = cell(nLens,1);
for i = 1:nLens
    mix = params.profileMix;
    ok = false;
    for attempt = 1:200
        xr = sort(rand(1,5));              % REQUIRE_MONOTONIC_X 와 동일: 오름차순
        yr = 1.5 * rand(1,5);              % pareto lb/ub 와 동일: y ∈ [0, 1.5]
        cx = (1-mix)*hemiX + mix*xr;
        cy = (1-mix)*hemiY + mix*yr;
        cand = [cx, cy];
        if isValidPoints(cand)
            ok = true; break;
        end
    end
    if ~ok
        cand = [hemiX, hemiY];             % fallback: 반구 (항상 유효)
    end
    profiles{i} = cand;                    % [x2..x6, y2..y6]
end

%% ===== 3) 높이장 합성 =====
g = linspace(-Ra, Ra, n);
[X, Y] = meshgrid(g, g);
H = zeros(n, n);
rho_s = linspace(0, 1, 101);               % 프로파일 평가용 반경 샘플
for i = 1:nLens
    c = profiles{i};
    px = [0, c(1:5), 1];                   % 끝점 고정: (0,1) apex, (1,0) rim
    py = [1, c(6:10), 0];
    % 중복 x 방지 (pchip 은 strictly increasing 필요)
    px = px + (0:6)*1e-9;
    zprof = pchip(px, py, rho_s);
    zprof = max(zprof, 0);  zprof(end) = 0;
    h_i = aspList(i) * rList(i);           % 높이 = 종횡비 * 반경 (반구: aspect=1)
    R = hypot(X - centers(i,1), Y - centers(i,2)) / rList(i);
    Hi = zeros(n, n);
    inside = R <= 1;
    Hi(inside) = h_i * interp1(rho_s, zprof, R(inside), 'linear', 0);
    H = max(H, Hi);                        % 비중첩이 보장되지만 max 로 안전 합성
end

%% ===== 4) .ent 쓰기 (generate_freeform_v2_ents.m 의 writer 와 동일) =====
write_freeform_ent(H, X, Y, n, params.tbase, params.templatePath, outPath);

%% ===== info =====
info = struct();
info.nLens = nLens;  info.pitch = p;  info.rMean = r_mean;
info.rList = rList;  info.aspList = aspList;  info.centers = centers;
info.fillRealized = sum(pi*rList.^2) / (2*Ra)^2;
info.jmax = jmax;    info.H = H;  info.grid = g;
end


%% ================= .ent writer =================
% generate_freeform_v2_ents.m 의 write_freeform_ent() 를 그대로 복사 (검증된 포맷).
function write_freeform_ent(H, X, Y, n, tbase, templatePath, outPath)
    % [방향] +Z 방향으로 볼록 돌출, 플랫 베이스는 아래(-z):
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
    % RearSurface(평면) z 위치 = -tbase
    newtxt = regexprep(newtxt, ...
        '(CSGLensSurfacePrimitive_1[\s\S]*?setPosition:  \{ 0\. 0\. )[-0-9.eE]+(  \} ;)', ...
        ['$1' num2str(-tbase,'%g') '$2'], 'once');
    newtxt = regexprep(newtxt, 'restoreSmoothResample: "Yes"', 'restoreSmoothResample: "No"', 'once');
    fid = fopen(outPath, 'w');
    if fid < 0, error('generate_random_supercell_ent:io', '출력 파일 열기 실패: %s', outPath); end
    fwrite(fid, newtxt);  fclose(fid);
end


%% ===== Spline 제약 (pareto_front_freeform.m 의 isValidPoints 와 동일) =====
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

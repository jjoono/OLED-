%% sharp_cut_freeform_ent.m
% 돔을 평면으로 실제로 잘라내(min 연산) 진짜 크리즈(꺾임)를 만드는 예시.
% DOF = 5 : h, z0, m, phi0, tbase  (20 DOF 이내)
%
% [중요/원인분석] 이전 극좌표(원판) 격자는 "자체교차" 에러가 났다. 이유:
%   템플릿의 표면은 setElementShape:"Rectangular", setUseSurfaceBoundary:"Yes" 로
%   정의되어 있어 데이터가 "모서리(코너) 있는 사각 도메인"이어야 한다. 원형 조리개는
%   LensElement 의 별도 설정(ElementShape:"Circular")이 사각 데이터를 사후에 잘라내는
%   것이지, 점 데이터 자체를 원판으로 만들 필요가 없다(오히려 만들면 깨짐). 그래서
%   여기서는 정사각 Cartesian 격자(코너 포함)를 쓰고, 그 위에서 min(돔,평면)만 적용한다.

h     = 1.0;   % 돔 높이
z0    = 0.45;  % 절단평면의 중심 높이 (낮을수록 잘리는 면적이 커짐)
m     = 1.2;   % 절단평면 기울기 (클수록 크리즈가 급함)
phi0  = 0;     % 절단평면 방향 [rad]
tbase = 0.02;  % base 두께 (밑받침 collar 두께. 0에 가까울수록 안 보이지만 완전히 0은
               % 두께0 degenerate solid라 불가 - 이 정도가 실용적 최소)
Ra    = 1.2139;  % 템플릿과 동일한 사각 격자 반폭(코너 포함); 광학 조리개는 별도로 원형 1.0

% [주의] n 은 렌더링 해상도일 뿐 DOF(자유변수)가 아니다. DOF는 위 5개(h,z0,m,phi0,tbase)
% 로 고정이며, n 을 키워도 BO 변수 수는 그대로다. 원형 트리밍 경계의 각짐(물결)은 격자와
% 원의 근본적 불일치에서 오는 한계지만, n 을 키울수록 계속 작아진다(공짜로 완화 가능).
n = 71;   % 격자 해상도(자유변수 아님). 더 매끈하게 하려면 더 키워도 됨(파일 크기만 증가)
templatePath = 'freeform_template.ent';
outPath = 'ff_sharp_cut.1.ent';

g = linspace(-Ra, Ra, n);
[X, Y] = meshgrid(g, g);
r = hypot(X, Y);

dome  = h * sqrt(max(1 - (r/1.0).^2, 0));         % 조리개 반경 1.0 기준 돔
plane = z0 + m * (X*cos(phi0) + Y*sin(phi0));
H = max(min(dome, plane), 0);
Z = -H;   % 음수 = 볼록(양각) 규약

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

% FrontSurface(첫 번째)만 SmoothResample 끔 -> 내부 스무딩이 크리즈를 뭉개는 것 방지.
% (파일 안 첫 occurrence = Front. Rear 는 어차피 평면이라 영향 없음.)
newtxt = regexprep(newtxt, 'restoreSmoothResample: "Yes"', 'restoreSmoothResample: "No"', 'once');

fid = fopen(outPath, 'w');
fwrite(fid, newtxt);
fclose(fid);

fprintf('저장됨: %s (사각 격자 %dx%d, 코너 포함, DOF=5)\n', outPath, n, n);

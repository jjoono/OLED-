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
z0    = 0.75;  % 절단평면의 중심 높이
m     = 0.6;   % 절단평면 기울기
phi0  = 0;     % 절단평면 방향 [rad]
tbase = 0.1;   % base 두께
Ra    = 1.2139;  % 템플릿과 동일한 사각 격자 반폭(코너 포함); 광학 조리개는 별도로 원형 1.0

n = 15;   % 격자 해상도(자유변수 아님)
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

fid = fopen(outPath, 'w');
fwrite(fid, newtxt);
fclose(fid);

fprintf('저장됨: %s (사각 격자 %dx%d, 코너 포함, DOF=5)\n', outPath, n, n);

%% generate_halfball_ent.m
% Half-ball(반구) 렌즈 .ent 생성. v2 규약: 원점 템플릿, +Z 돌출, taper-rim(물결 억제).
% 대조군(hemisphere baseline)으로 사용.
%
% [기하 - 명시]
%   조리개 반경 Rap = 1.0,  돔 높이 h = 1.0  -> 진짜 반구(높이=반경)
%   +Z 방향 볼록: 반구가 z=0(rim) 에서 z=+1.0(정점) 까지 돌출
%   밑받침(base collar) 두께 tbase = 0.30  -> RearSurface 평면이 z = -0.30 에 위치
%     => 전체 solid 높이 = h + tbase = 1.30, 그중 0.30 이 평평한 base (시뮬에서 파묻힘)
%   테두리 taper = 0.03 (바깥 3% 반경만 접선으로 완만화 -> rim 특이점/물결 제거)

templatePath = 'freeform_template_v2.ent';
outPath      = 'ff_halfball.1.ent';

Rap   = 1.0;      % 조리개 반경 = 반구 반경
h     = 1.0;      % 돔 높이 = Rap -> 반구
tbase = 0.30;     % 밑받침 두께 (RearSurface 평면 z = -tbase)
taper = 0.03;     % rim 완만화 비율
Ra    = 1.2139;   % 사각 격자 반폭 (템플릿과 동일, 모서리 포함)
n     = 81;       % 렌더 해상도 (DOF 아님)

g = linspace(-Ra, Ra, n);
[X, Y] = meshgrid(g, g);
r = hypot(X, Y);

% taper-rim 반구
Pb = sqrt(max(1-(r/Rap).^2, 0));
r0 = Rap*(1-taper);
sag0   = sqrt(max(1-(r0/Rap)^2, 0));
slope0 = -(r0/Rap)/sqrt(max(1-(r0/Rap)^2, 1e-9));
Plin = max(sag0 + slope0*(r-r0), 0);
P = Pb;  P(r>=r0) = Plin(r>=r0);
H = h .* P;

fprintf('Half-ball: Rap=%.2f, h=%.2f (+Z 돌출), tbase=%.2f (base plane z=%.2f), taper=%.2f\n', ...
    Rap, h, tbase, -tbase, taper);
fprintf('  -> solid: z=-%.2f(flat base) ~ z=0(rim) ~ z=+%.2f(apex), 전체높이 %.2f\n', ...
    tbase, h, h+tbase);

%% ---- .ent 쓰기 (v2 규약: Z=+H, rear=-tbase, 원점 템플릿) ----
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
newtxt = regexprep(newtxt, ...
    '(CSGLensSurfacePrimitive_1[\s\S]*?setPosition:  \{ 0\. 0\. )[-0-9.eE]+(  \} ;)', ...
    ['$1' num2str(-tbase,'%g') '$2'], 'once');
newtxt = regexprep(newtxt, 'restoreSmoothResample: "Yes"', 'restoreSmoothResample: "No"', 'once');
fid = fopen(outPath, 'w');  fwrite(fid, newtxt);  fclose(fid);
fprintf('저장됨: %s\n', outPath);

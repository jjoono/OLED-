%% sharp_cut_freeform_ent.m
% 돔을 평면으로 실제로 잘라내(min 연산) 진짜 크리즈(꺾임)를 만드는 예시.
% 다항식/스플라인 근사가 아니라 기하학적 절단이라 "칼로 자른 단면"이 정확히 생김.
% 자유변수(DOF) = 5 : h, z0, m, phi0, tbase  (20 DOF 이내)
% 극좌표 격자(각도x반지름) 사용 -> X,Y 배치가 안 꼬여 자체교차 안전 (polar_freeform_ent.m과 동일 원리).

h     = 1.0;   % 돔 높이
z0    = 0.75;  % 절단평면의 중심 높이
m     = 0.6;   % 절단평면 기울기
phi0  = 0;     % 절단평면 방향 [rad]
tbase = 0.1;   % base 두께

nr = 24; nt = 48;   % 렌더링 해상도(자유변수 아님, 크리즈를 촘촘히 표현하기 위함)
Ra = 1.0;

templatePath = 'freeform_template.ent';
outPath = 'ff_sharp_cut.1.ent';

r  = linspace(Ra/(nr+1), Ra, nr);
th = linspace(0, 2*pi, nt+1); th(end) = [];

X = zeros(nt, nr); Y = zeros(nt, nr); H = zeros(nt, nr);
for iv = 1:nr
    xr = r(iv) * cos(th)';  yr = r(iv) * sin(th)';
    dome  = h * sqrt(max(1 - (r(iv)/Ra)^2, 0));
    plane = z0 + m * (xr*cos(phi0) + yr*sin(phi0));
    hcut  = max(min(dome, plane), 0);
    X(:,iv) = xr; Y(:,iv) = yr; H(:,iv) = hcut;
end
Z = -H;   % 음수 = 볼록(양각) 규약

Xv = X(:); Yv = Y(:); Zv = Z(:);
N = nt * nr;

tpl = fileread(templatePath);
tok = regexp(tpl, 'ORAStartData;([\s\S]*?)ORAEndData;', 'tokenExtents');
s0 = tok{1}(1); e0 = tok{1}(2);

buf = sprintf('0 1 %d %d 0 0 %d 0 0 0', nt, nr, N);
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

fprintf('저장됨: %s (DOF=5: h=%.2f,z0=%.2f,m=%.2f,phi0=%.2f,tbase=%.2f | 격자 %dx%d)\n', ...
    outPath, h, z0, m, phi0, tbase, nt, nr);

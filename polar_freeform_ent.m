%% polar_freeform_ent.m
% 극좌표(반지름 x 각도) 격자로 FrontSurface 를 만들어 .ent 로 저장.
% 원리: U=각도(고정 V에서 등분) -> 경계가 자동으로 원. V=반지름(고정 U에서 단조증가)
%       -> 인덱스 이웃 = 공간 이웃이 항상 성립 -> 자체교차 불가능.
% 중심 특이점 방지를 위해 r_min>0 부터 시작(중앙에 작은 평평한 원판만 남음).

nr = 8;      % 반지름 링 수 (V)
nt = 16;     % 각도 분할 수 (U)
Ra = 1.0;    % 조리개 반지름
r_min = Ra/(nr+1);   % 중심 특이점(다중좌표 중복) 방지
h = 1.0;     % 돔 높이
tbase = 0.1; % base 두께

templatePath = 'freeform_template.ent';
outPath = 'ff_polar_hemi.1.ent';

r = linspace(r_min, Ra, nr);
th = linspace(0, 2*pi, nt+1); th(end) = [];   % 0..2pi, 중복 각도 제외

X = zeros(nt, nr); Y = zeros(nt, nr); Z = zeros(nt, nr);
for iv = 1:nr
    X(:,iv) = r(iv) * cos(th)';
    Y(:,iv) = r(iv) * sin(th)';
    Z(:,iv) = -h * sqrt(max(1 - (r(iv)/Ra)^2, 0));   % 음수 Z = 볼록(양각), rim에서 0
end

% U(각도, 빠름) x V(반지름, 느림) 순서로 평탄화 = 기존 템플릿과 동일한 순서 규칙
Xv = X(:); Yv = Y(:); Zv = Z(:);
N = nt * nr;

tpl = fileread(templatePath);
tok = regexp(tpl, 'ORAStartData;([\s\S]*?)ORAEndData;', 'tokenExtents');
s0 = tok{1}(1); e0 = tok{1}(2);   % FrontSurface 블록만 교체

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

fprintf('저장됨: %s (극좌표 격자 %dx%d, r_min=%.3f, base=%.2f)\n', outPath, nt, nr, r_min, tbase);

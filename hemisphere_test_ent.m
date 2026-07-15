%% hemisphere_test_ent.m
% 반구형 렌즈 .ent 2종 생성 (밑받침은 넉넉하게: 시뮬 코드에서 오프셋 보정 가정).
%   (1) ff_hemi_true.1.ent      : 수학적으로 순수한 반구. 테두리(적도) 기울기가
%       무한대(특이점)라, 격자를 아무리 촘촘히 해도 테두리 근처 물결은 "0"이 아니라
%       느리게(제곱근 특이점 -> O(1/sqrt(n)) 수준) 수렴만 한다.
%   (2) ff_hemi_taperrim.1.ent  : 바깥쪽 3% 반경 구간만 그 지점 접선으로 완만화해
%       기울기를 유한하게 만든 버전. 광학적으로 거의 동일하지만(실제 제작 렌즈도
%       테두리가 완벽한 수직은 아님) 격자 수렴이 훨씬 빠를 것으로 기대.
%
% [주의] tbase(밑받침 두께)와 테두리 물결은 서로 독립적인 변수다. tbase 를 늘려도
% 줄여도 물결에는 영향 없음 - 물결은 오직 격자 해상도 n 과, (1)의 경우 특이점 자체가
% 결정한다. 어느 쪽이 실제로 물결이 적은지는 LightTools 에 직접 로드해 비교할 것
% (내부 보간 방식을 여기서 정확히 재현할 수 없어 이론 추정만으로는 단정 불가).

Ra = 1.2139; Rap = 1.0; h = 1.0; tbase = 0.3; n = 81;
templatePath = 'freeform_template.ent';

g = linspace(-Ra, Ra, n);
[X, Y] = meshgrid(g, g);
r = hypot(X, Y);

% (1) 순수 반구
H1 = h * sqrt(max(1 - (r/Rap).^2, 0));
write_hemi_ent(H1, X, Y, n, tbase, templatePath, 'ff_hemi_true.1.ent');

% (2) 테두리 3% 완만화(접선 연장, 유한 기울기)
taper = 0.03; r0 = Rap*(1-taper);
H0 = h*sqrt(max(1-(r0/Rap)^2,0));
slope0 = -h*(r0/Rap)/sqrt(max(1-(r0/Rap)^2,1e-9));
lin = max(H0 + slope0*(r-r0), 0);
H2 = H1; H2(r>=r0) = lin(r>=r0);
write_hemi_ent(H2, X, Y, n, tbase, templatePath, 'ff_hemi_taperrim.1.ent');

fprintf('두 파일 생성 완료. LightTools 로 각각 로드해 테두리 물결을 직접 비교.\n');


function write_hemi_ent(H, X, Y, n, tbase, templatePath, outPath)
    Z = -H;
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
    fprintf('저장됨: %s\n', outPath);
end

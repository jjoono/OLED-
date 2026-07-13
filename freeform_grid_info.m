function [XY, innerIdx, N] = freeform_grid_info(templatePath)
% FREEFORM_GRID_INFO  freeform 템플릿(.ent)의 FrontSurface 격자 정보를 파싱.
%
%   [XY, innerIdx, N] = freeform_grid_info(templatePath)
%
%   출력
%     XY       : N x 2, 각 격자점의 (X,Y) [정규 좌표]
%     innerIdx : 내부(경계가 아닌) 점들의 인덱스. 경계점은 Z=0 고정(기판에 맞물림),
%                내부점만 BO 변수로 삼는다. (규칙 5x5 -> 내부 3x3 = 9개)
%     N        : 총 격자점 수 (예: 25)
%
%   경계 판정: 바깥 링(|X| 또는 |Y| 가 최대값)에 있는 점 = 경계.

    tpl = fileread(templatePath);
    tok = regexp(tpl, 'ORAStartData;([\s\S]*?)ORAEndData;', 'tokens', 'once');
    if isempty(tok)
        error('freeform_grid_info:tpl', '템플릿에서 데이터 블록을 못 찾음: %s', templatePath);
    end
    toks = strsplit(strtrim(tok{1}));
    N = round(str2double(toks{7}));
    XY = zeros(N, 2);
    for i = 1:N
        XY(i,1) = str2double(toks{10 + 3*(i-1) + 1});
        XY(i,2) = str2double(toks{10 + 3*(i-1) + 2});
    end

    xmax = max(abs(XY(:,1)));  ymax = max(abs(XY(:,2)));
    tol = 1e-6 * max(xmax, 1);
    onBoundary = (abs(abs(XY(:,1)) - xmax) < tol) | (abs(abs(XY(:,2)) - ymax) < tol);
    innerIdx = find(~onBoundary);
end

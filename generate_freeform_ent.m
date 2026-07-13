function [outPath, XY] = generate_freeform_ent(zFront, templatePath, outPath)
% GENERATE_FREEFORM_ENT  자유형(Freeform) LensElement 라이브러리 파일(.ent)을 직접 생성.
%
%   [outPath, XY] = generate_freeform_ent(zFront, templatePath, outPath)
%
%   [핵심] LightTools FreeformEntity(.ent) 는 표면 격자점을 X,Y,Z 텍스트로 그대로
%   담는다(restorePoints 블록). 따라서 STL import / SaveLibrary 없이 MATLAB 이 이
%   .ent 를 직접 써서, 기존 코드가 swept_XXX.1.ent 를 texture unit-cell 에 물리던 것과
%   똑같이 물리면 된다.
%
%   사용자가 GUI 에서 만든 실제 .ent(templatePath)를 "템플릿"으로 두고, 첫 번째
%   자유형 표면(FrontSurface) 격자의 Z 값만 zFront 로 치환한다. X,Y 및 나머지 구조
%   (헤더/mapper/RearSurface 평면/재질/토폴로지)는 검증된 실제 파일 그대로 → 포맷
%   오류 위험 최소. (FrontSurface=자유형 곡면, RearSurface=평면 → 바닥이 기판에
%   파묻힌 plano-freeform 렌즈.)
%
%   입력
%     zFront       : 길이 N(=FrontSurface 격자점 수, 예 4x4=16) Z 벡터 [정규 단위]
%     templatePath : 기준 .ent (freeform_template.ent)
%     outPath      : 출력 .ent 경로
%   출력
%     outPath : 쓰인 경로
%     XY      : Nx2, 각 격자점 (X,Y) [정규 좌표] — report/미리보기/유효성용

    tpl = fileread(templatePath);

    % 첫 번째 ORAStartData;...ORAEndData; = FrontSurface 데이터 블록의 본문 범위
    te = regexp(tpl, 'ORAStartData;([\s\S]*?)ORAEndData;', 'tokenExtents', 'once');
    if isempty(te)
        error('generate_freeform_ent:tpl', '템플릿에서 데이터 블록을 못 찾음: %s', templatePath);
    end
    body = tpl(te(1):te(2));
    toks = strsplit(strtrim(body));

    N = round(str2double(toks{7}));     % 헤더 7번째 토큰 = 점 개수
    hdr    = toks(1:10);
    pts    = toks(11:10+3*N);
    footer = toks(11+3*N:end);          % 예: '0 0 4 CartesianMapper 1 0 0 0 0'

    XY = zeros(N,2);
    for i = 1:N
        XY(i,1) = str2double(pts{3*(i-1)+1});
        XY(i,2) = str2double(pts{3*(i-1)+2});
    end

    if numel(zFront) ~= N
        error('generate_freeform_ent:N', 'zFront 길이(%d) != 템플릿 격자점 수(%d)', numel(zFront), N);
    end

    % 새 본문: 헤더 + (X Y Znew) x N + footer  (모두 공백 구분)
    buf = strjoin(hdr, ' ');
    for i = 1:N
        buf = [buf sprintf(' %.17g %.17g %.17g', XY(i,1), XY(i,2), zFront(i))]; %#ok<AGROW>
    end
    buf = [buf ' ' strjoin(footer, ' ')];
    newbody = [char(10) buf char(10)];   % 앞뒤 개행 (원본과 동일한 배치)

    newtxt = [tpl(1:te(1)-1), newbody, tpl(te(2)+1:end)];

    fid = fopen(outPath, 'w');
    if fid < 0, error('generate_freeform_ent:io', '출력 파일 열기 실패: %s', outPath); end
    fwrite(fid, newtxt);
    fclose(fid);
end

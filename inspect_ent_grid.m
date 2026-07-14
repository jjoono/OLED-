%% inspect_ent_grid.m
% .ent 파일의 자유형 표면 격자를 읽어 X, Y, Z 를 각각 U x V 행렬로 저장.
% 실행하면 파일 선택창이 열림. FrontSurface = 블록1, RearSurface = 블록2.
%
% 결과 변수: X1,Y1,Z1 (FrontSurface, U1 x V1), X2,Y2,Z2 (RearSurface, U2 x V2)

[fname, fpath] = uigetfile({'*.ent', 'LightTools Library (*.ent)'; '*.*', 'All Files'}, ...
    '.ent 파일 선택');
if isequal(fname, 0)
    disp('선택 취소됨.');
    return;
end
entPath = fullfile(fpath, fname);

txt = fileread(entPath);
toks = regexp(txt, 'ORAStartData;([\s\S]*?)ORAEndData;', 'tokens');
if isempty(toks)
    error('데이터 블록을 찾지 못함: %s', entPath);
end

blockNames = {'FrontSurface', 'RearSurface'};

for b = 1:numel(toks)
    t = strsplit(strtrim(toks{b}{1}));
    U = round(str2double(t{3}));
    V = round(str2double(t{4}));
    N = round(str2double(t{7}));

    Xv = zeros(N,1); Yv = zeros(N,1); Zv = zeros(N,1);
    for i = 1:N
        Xv(i) = str2double(t{10 + 3*(i-1) + 1});
        Yv(i) = str2double(t{10 + 3*(i-1) + 2});
        Zv(i) = str2double(t{10 + 3*(i-1) + 3});
    end

    % 점 순서는 U가 안쪽 루프(행 우선)라고 가정 -> U x V 행렬로 reshape
    X = reshape(Xv, U, V).';
    Y = reshape(Yv, U, V).';
    Z = reshape(Zv, U, V).';

    if b <= numel(blockNames)
        name = blockNames{b};
    else
        name = sprintf('Block%d', b);
    end
    fprintf('\n=== %s : %d x %d ===\n', name, U, V);
    disp('X ='); disp(X);
    disp('Y ='); disp(Y);
    disp('Z ='); disp(Z);

    assignin('base', sprintf('X%d', b), X);
    assignin('base', sprintf('Y%d', b), Y);
    assignin('base', sprintf('Z%d', b), Z);
end

fprintf('\n워크스페이스 변수: X1,Y1,Z1 (FrontSurface), X2,Y2,Z2 (RearSurface)\n');

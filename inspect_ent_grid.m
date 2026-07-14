function inspect_ent_grid(entPath)
% INSPECT_ENT_GRID  .ent 파일의 자유형 표면 격자 X,Y,Z 값을 표로 출력.
%
%   inspect_ent_grid()               % 인자 없이 실행 -> 파일 선택창이 열림
%   inspect_ent_grid('freeform_template.ent')   % 경로 직접 지정
%
%   .ent 안의 ORAStartData...ORAEndData 블록(보통 2개: Front, Rear)을 찾아
%   각 블록의 U,V 격자 크기와 (X,Y,Z) 점들을 표로 보여준다. 디버깅용.

    if nargin < 1 || isempty(entPath)
        [fname, fpath] = uigetfile({'*.ent', 'LightTools Library (*.ent)'; '*.*', 'All Files'}, ...
            '.ent 파일 선택');
        if isequal(fname, 0)
            disp('선택 취소됨.');
            return;
        end
        entPath = fullfile(fpath, fname);
    end

    txt = fileread(entPath);
    toks = regexp(txt, 'ORAStartData;([\s\S]*?)ORAEndData;', 'tokens');
    if isempty(toks)
        error('inspect_ent_grid:none', '데이터 블록을 찾지 못함: %s', entPath);
    end

    names = {'FrontSurface(1번째 블록)', 'RearSurface(2번째 블록)', '블록3', '블록4'};
    for b = 1:numel(toks)
        t = strsplit(strtrim(toks{b}{1}));
        U = round(str2double(t{3}));
        V = round(str2double(t{4}));
        N = round(str2double(t{7}));
        fprintf('\n=== %s : U=%d V=%d N=%d ===\n', names{min(b,numel(names))}, U, V, N);

        X = zeros(N,1); Y = zeros(N,1); Z = zeros(N,1);
        for i = 1:N
            X(i) = str2double(t{10 + 3*(i-1) + 1});
            Y(i) = str2double(t{10 + 3*(i-1) + 2});
            Z(i) = str2double(t{10 + 3*(i-1) + 3});
        end

        T = table((1:N)', X, Y, Z, 'VariableNames', {'idx','X','Y','Z'});
        disp(T);
        fprintf('X range [%.4f, %.4f] | Y range [%.4f, %.4f] | Z range [%.4f, %.4f]\n', ...
            min(X), max(X), min(Y), max(Y), min(Z), max(Z));
    end
end

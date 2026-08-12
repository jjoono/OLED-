% ============================================================
%  probe_texture_keys.m
%
%  목적: LightTools 3D texture 의 양각/음각(bump vs hole) 토글이
%        데이터베이스에서 어떤 객체의 어떤 속성 이름으로 노출되는지 찾아낸다.
%
%  배경: 검증된 스크립트들은 텍스처 관련해서 아래 두 가지만 사용해 왔다.
%          LIBRARY_ELEMENT_UNIT_CELL / 'LibraryElement' -> 'Filename'
%          TEXTURE_PARAMETER         / 'StretchZ'       -> 'Value'
%        bump/hole 토글의 DB 이름 선례가 레포에 없으므로, 실행 전에 이 스크립트로
%        실제 이름을 확인한 뒤 stress_inverted_mla.m 상단 상수에 적어 넣는다.
%
%  사용법: LightTools 를 평소처럼 (렌즈 모델이 열린 상태로) 띄워 두고 실행.
%          출력에서 Bump / Hole / Relief / Protrusion / Depression / Subtract /
%          Boolean / CellType 등으로 보이는 항목을 찾으면 된다.
%          찾은 이름을 stress_inverted_mla.m 의 TEX_RELIEF_PARAM /
%          TEX_RELIEF_VALUE_HOLE 에 적는다.
%
%  주의: 이 스크립트는 읽기 전용이다. 아무 값도 바꾸지 않는다.
% ============================================================
clear;
global ID_swept ID_LT ltml ltloc

RenewLightTools();
try
    ltml.LTCmd(ltml.GetLTAPI(ID_LT), 'Message "Check Connection"');
catch
    ltml = actxserver('ltcom64.LTAPI2');
    ltloc = actxserver('ltlocator.Locator');
end
lt2 = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt2, "ShowFileDialogBox", 0);

% 텍스처와 관련될 수 있는 DB 목록들. 존재하지 않는 것은 조용히 건너뛴다.
LIST_TYPES = { ...
    'TEXTURE_PARAMETER', ...
    'LIBRARY_ELEMENT_UNIT_CELL', ...
    'TEXTURE', ...
    'TEXTURE_ZONE', ...
    'TEXTURE_ZONE_EXTENT', ...
    'SURFACE_TEXTURE', ...
    'TEXTURE_UNIT_CELL' };

% 각 항목에서 시도해 볼 속성 이름 후보. 읽히는 것만 출력된다.
PROP_CANDIDATES = { ...
    'Name', 'Value', 'Filename', ...
    'Bump', 'Hole', 'BumpOrHole', 'BumpHole', ...
    'Relief', 'ReliefType', 'CellType', 'TextureType', 'UnitCellType', ...
    'Protrusion', 'Depression', 'Polarity', 'Sense', 'Invert', 'Inverted', ...
    'Subtract', 'SubtractFromSurface', 'BooleanOperation', 'Operation', ...
    'AddOrSubtract', 'Direction', 'Sign' };

MAX_ITEMS = 40;   % 목록당 훑어볼 최대 개수

fprintf('\n================ LightTools texture DB probe ================\n');
for t = 1:numel(LIST_TYPES)
    typeName = LIST_TYPES{t};
    try
        List = ltml.LTDbList(lt2, 'LENS_MANAGER[1]', typeName);
    catch
        continue;   % 이 모델에 없는 목록 타입
    end
    if isempty(List) || (isnumeric(List) && List == 0)
        continue;
    end

    fprintf('\n--- %s ---\n', typeName);
    nFound = 0;
    for i = 1:MAX_ITEMS
        try
            Key = ltml.LTListAtPos(lt2, List, i);
        catch
            break;      % 목록 끝
        end
        if isempty(Key) || (isnumeric(Key) && Key == 0)
            break;
        end
        nFound = nFound + 1;

        % 이 항목에서 읽히는 속성만 모아 한 줄로 출력
        parts = {};
        for p = 1:numel(PROP_CANDIDATES)
            prop = PROP_CANDIDATES{p};
            try
                v = ltml.LTDbGet(lt2, Key, prop);
            catch
                continue;   % 이 항목에 없는 속성
            end
            if isempty(v)
                continue;
            end
            if isnumeric(v) || islogical(v)
                parts{end+1} = sprintf('%s=%g', prop, double(v(1))); %#ok<SAGROW>
            else
                parts{end+1} = sprintf('%s=%s', prop, char(string(v))); %#ok<SAGROW>
            end
        end

        if isempty(parts)
            fprintf('  [%2d] (읽히는 후보 속성 없음)\n', i);
        else
            fprintf('  [%2d] %s\n', i, strjoin(parts, '  |  '));
        end
    end
    if nFound == 0
        fprintf('  (항목 없음)\n');
    end
end

fprintf('\n============================================================\n');
fprintf('찾는 것: bump/hole (양각/음각) 를 나타내는 속성.\n');
fprintf('위 출력에서 해당 이름과 값(음각일 때의 값)을 확인한 뒤\n');
fprintf('stress_inverted_mla.m 상단의\n');
fprintf('    TEX_RELIEF_PARAM      (TEXTURE_PARAMETER 항목 이름)\n');
fprintf('    TEX_RELIEF_VALUE_HOLE (음각에 해당하는 값)\n');
fprintf('에 적어 넣으면 된다.\n');
fprintf('후보 목록에 없는 이름이면 LightTools GUI 의 Database Browser 에서\n');
fprintf('해당 텍스처 속성을 직접 확인하는 것이 가장 빠르다.\n');
fprintf('============================================================\n\n');

%% ===== RenewLightTools (검증된 원본 그대로) =====
function RenewLightTools()
global ID_swept ID_LT ltml ltloc
ltml = actxserver('ltcom64.LTAPI2');
ltloc = actxserver('ltlocator.Locator');
[~, cmdout] = system('tasklist /FI "IMAGENAME eq lt.exe" /FO CSV /NH');
lines = strsplit(strtrim(cmdout), '\n');
pids = [];
for i = 1:numel(lines)
    tok = strsplit(lines{i}, ',');
    if numel(tok) >= 2
        p = str2double(erase(tok{2}, '"'));
        if ~isnan(p), pids(end+1) = p; end %#ok<AGROW>
    end
end
if numel(pids) < 2
    error(['lt.exe 프로세스가 2개 필요하다 (swept 편집용 + 어레이 모델용). ' ...
           '현재 %d개. 평소 최적화 스크립트와 동일하게 띄워 둘 것.'], numel(pids));
end
ID_swept = pids(1);
ID_LT    = pids(2);
end
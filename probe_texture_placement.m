% ============================================================
%  probe_texture_placement.m
%
%  목적: 텍스처 배치(ZoneTextureHexagonalPlacement)의 **실제 DB 속성 이름**을
%        찾아낸다. .lts 파일에는 setXSpacing / setYSpacing 으로 저장되지만,
%        이 모델은 파일 토큰과 DB 이름이 다른 전례가 있다
%        (zone: 파일 setZoneWidth  <->  DB Geometry_1).
%
%  사용법: LightTools 를 평소처럼 띄운 상태에서 실행. **읽기 전용**이 아니라
%        값을 썼다 되돌리므로(원복함), 시뮬 중이 아닐 때 돌릴 것.
%
%  출력에서 "쓰기→되읽기 일치" 로 표시된 (목록, 위치, 속성) 조합을 찾으면 된다.
%  그 조합을 set_texture_spacing.m 의 listNames / propPairs 맨 앞에 넣으면
%  이후 실행이 곧바로 그 경로를 쓴다.
% ============================================================
clear;
global ID_LT ltml ltloc

RenewLightTools();
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

listNames = {'ZONE_TEXTURE_HEXAGONAL_PLACEMENT','HEXAGONAL_PLACEMENT', ...
             'ZONE_TEXTURE_PLACEMENT','TEXTURE_PLACEMENT','PLACEMENT', ...
             'VARIABLE_SPACED_TEXTURE','TEXTURE_ZONE_EXTENT','TEXTURE_PARAMETER'};

props = {'Name','XSpacing','YSpacing','XOffset','YOffset', ...
         'Geometry_1','Geometry_2','Geometry_3','Geometry_4', ...
         'X_Spacing','Y_Spacing','SpacingX','SpacingY', ...
         'XPitch','YPitch','Pitch_1','Pitch_2','Spacing_1','Spacing_2', ...
         'Value','Value_1','Value_2','DeltaX','DeltaY','Type','PlacementType'};

TESTVAL = 0.123456;   % 쓰기 테스트용 값 (성공 시 원복)

% --- 부모 후보: 최상위 + zone 키 ---
parents = {[], 'zone'};

fprintf('\n================ TEXTURE PLACEMENT PROBE ================\n');
for ip = 1:numel(parents)
    if isempty(parents{ip})
        pdesc = 'lens_manager[1] (top-level)';  PK = [];
    else
        pdesc = 'zone key (child)';
        try
            ZL = ltml.LTDbList(lt,'lens_manager[1]','TEXTURE_ZONE_EXTENT');
            PK = ltml.LTListByName(lt, ZL, 'zone');
        catch
            fprintf('\n--- %s: zone 키 획득 실패, 건너뜀 ---\n', pdesc);
            continue;
        end
    end
    fprintf('\n########## 부모: %s ##########\n', pdesc);

    for il = 1:numel(listNames)
        try
            if isempty(PK)
                L = ltml.LTDbList(lt, 'lens_manager[1]', listNames{il});
            else
                L = ltml.LTDbList(lt, PK, listNames{il});
            end
        catch
            continue;
        end
        if isempty(L) || (isnumeric(L) && L == 0), continue; end

        printedHeader = false;
        for pos = 1:12   % Scale/StretchX/Y/Z 는 6번 뒤에 온다
            try
                K = ltml.LTListAtPos(lt, L, pos);
            catch
                break;
            end
            if isempty(K) || (isnumeric(K) && K == 0), break; end

            if ~printedHeader
                fprintf('\n--- 목록 %s ---\n', listNames{il});
                printedHeader = true;
            end

            % 1) 읽히는 속성 나열
            readable = {};
            for q = 1:numel(props)
                try
                    v = ltml.LTDbGet(lt, K, props{q});
                catch
                    continue;
                end
                if isempty(v), continue; end
                if isnumeric(v) || islogical(v)
                    if ~isfinite(double(v(1))), continue; end     % NaN = 없는 속성
                    readable{end+1} = sprintf('%s=%g', props{q}, double(v(1))); %#ok<SAGROW>
                else
                    readable{end+1} = sprintf('%s=%s', props{q}, char(string(v))); %#ok<SAGROW>
                end
            end
            if isempty(readable)
                fprintf('  [pos %d] 읽히는 속성 없음\n', pos);
            else
                fprintf('  [pos %d] %s\n', pos, strjoin(readable, '  |  '));
            end

            % 2) 쓰기→되읽기 검사 (수치 속성만)
            for q = 1:numel(props)
                nm = props{q};
                if any(strcmp(nm, {'Name','Type','PlacementType'})), continue; end
                try
                    old = ltml.LTDbGet(lt, K, nm);
                    if ~isnumeric(old) || ~isfinite(double(old)), continue; end
                    ltml.LTDbSet(lt, K, nm, TESTVAL);
                    chk = ltml.LTDbGet(lt, K, nm);
                    ltml.LTDbSet(lt, K, nm, old);        % 원복
                    if isfinite(chk) && abs(chk - TESTVAL) < 1e-9
                        fprintf('      >>> 쓰기→되읽기 일치: %s (원래값 %g)\n', nm, double(old));
                    end
                catch
                end
            end
        end
    end
end

fprintf('\n=========================================================\n');
fprintf('찾을 것: 배치 간격에 해당하는 두 속성 (원래값 0.0866 / 0.1 근처).\n');
fprintf('그 (목록, 속성쌍) 을 set_texture_spacing.m 의 listNames / propPairs\n');
fprintf('맨 앞에 넣으면 본 실행이 곧바로 그 경로를 사용한다.\n');
fprintf('=========================================================\n\n');

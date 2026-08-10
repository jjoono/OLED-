function [ok, desc] = set_texture_spacing(lt, tgtX, tgtY, verbose)
% SET_TEXTURE_SPACING  텍스처 배치 간격(XSpacing/YSpacing 상당)을 설정하고 검증.
%
%   [ok, desc] = set_texture_spacing(lt, tgtX, tgtY, verbose)
%
% [왜 필요한가] 슈퍼셀 하나에 nCols x nCols 개의 렌즈렛이 들어가므로, 배치 간격을
%   같은 배수로 키우지 않으면 렌즈렛 하나의 물리 크기가 1/nCols 로 줄어든다.
%   그러면 "무작위성의 효과" 가 아니라 "렌즈렛이 작아진 효과" 를 재게 된다.
%
% [왜 이렇게 복잡한가] 이 모델은 .lts 파일의 토큰과 DB 속성 이름이 다르다.
%   예: zone 은 파일에서 setZoneWidth/setZoneHeight 지만 DB 에서는
%   Geometry_1/Geometry_2 로 접근한다(검증된 기존 코드가 그렇게 쓴다).
%   placement 도 마찬가지일 수 있으므로, (목록 x 위치 x 속성명) 조합을 훑고
%   **되읽어서** 확인된 조합만 성공으로 친다.
%   LTDbGet 은 잘못된 키/이름에 예외 대신 NaN 을 돌려주므로 NaN 은 실패로 본다.
%
% 반환: ok=true 면 desc 에 성공한 조합 설명이 담긴다. 실패 시 ok=false.

global ltml
if nargin < 4, verbose = false; end
ok = false;  desc = '';
TOL = 1e-9;

% 후보 목록: (a) lens_manager 최상위, (b) zone 키의 자식
listNames = {'ZONE_TEXTURE_HEXAGONAL_PLACEMENT','HEXAGONAL_PLACEMENT', ...
             'ZONE_TEXTURE_PLACEMENT','TEXTURE_PLACEMENT','PLACEMENT'};

% 후보 속성 쌍 (파일 토큰 / Geometry_n 관례 / 흔한 변형)
propPairs = { {'XSpacing','YSpacing'}, {'Geometry_1','Geometry_2'}, ...
              {'X_Spacing','Y_Spacing'}, {'SpacingX','SpacingY'}, ...
              {'XPitch','YPitch'}, {'Pitch_1','Pitch_2'}, ...
              {'Spacing_1','Spacing_2'}, {'Value_1','Value_2'}, ...
              {'DeltaX','DeltaY'} };

% --- 부모 키 후보 수집: [] = 최상위, 그 외 = zone 키 ---
parents = {[]};
try
    ZL = ltml.LTDbList(lt,'lens_manager[1]','TEXTURE_ZONE_EXTENT');
    ZK = ltml.LTListByName(lt, ZL, 'zone');
    if ~isempty(ZK) && ~(isnumeric(ZK) && ZK == 0)
        parents{end+1} = ZK;
    end
catch
end

for ip = 1:numel(parents)
    for il = 1:numel(listNames)
        try
            if isempty(parents{ip})
                L = ltml.LTDbList(lt, 'lens_manager[1]', listNames{il});
            else
                L = ltml.LTDbList(lt, parents{ip}, listNames{il});
            end
        catch
            continue;
        end
        if isempty(L) || (isnumeric(L) && L == 0), continue; end

        for pos = 1:5
            try
                K = ltml.LTListAtPos(lt, L, pos);
            catch
                break;
            end
            if isempty(K) || (isnumeric(K) && K == 0), break; end

            for pp = 1:numel(propPairs)
                nx = propPairs{pp}{1};  ny = propPairs{pp}{2};
                try
                    ltml.LTDbSet(lt, K, nx, tgtX);
                    ltml.LTDbSet(lt, K, ny, tgtY);
                    gx = ltml.LTDbGet(lt, K, nx);
                    gy = ltml.LTDbGet(lt, K, ny);
                catch
                    continue;
                end
                if isfinite(gx) && isfinite(gy) && ...
                   abs(gx - tgtX) < TOL && abs(gy - tgtY) < TOL
                    ok = true;
                    desc = sprintf('%s[pos %d]%s .%s/.%s', listNames{il}, pos, ...
                        tern_parent(parents{ip}), nx, ny);
                    if verbose
                        fprintf('  [spacing] 설정 확인: %s = %.4f / %.4f\n', desc, gx, gy);
                    end
                    return;
                end
            end
        end
    end
end

function s = tern_parent(pk)
if isempty(pk), s = ' (top-level)'; else, s = ' (child of zone)'; end
end
end

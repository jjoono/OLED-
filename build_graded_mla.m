%% build_graded_mla.m
% 위치별로 다른 .ent(편심 grading) 를 LightTools 배열에 적용한다.
% .lts 는 직접 못 쓰므로(바이너리), 기존 파이프라인처럼 COM 으로 모델을 구성한다.
%
% [설계] 배열을 Nzone 개의 동심 반경 zone 으로 나누고, zone k 의 렌즐릿에는
%   ff_grad_z{k}.1.ent (편심 dx_k) 를 물린다. 중심->가장자리로 편심 증가 =
%   큰 렌즈를 프레넬화(집단 collimation) -> 전역 에텐듀 harvest 시도.
%
% [두 가지 배치 방법] 아래 METHOD 로 선택. 사용자 .lts 가 어느 쪽을 지원하는지에 따라.
%   'zones'   : 텍스처를 Nzone 개 두고 각 zone 의 unit-cell Filename 을 다르게.
%               (모델에 zone 별 텍스처가 미리 있어야 함. 이름 규약 @@VERIFY)
%   'entities': 렌즐릿을 (x,y) 마다 개별 FreeformEntity 로 복제 배치.
%               (엔티티 복제/이동 명령 @@VERIFY. 렌즐릿 수 많으면 느림)
%
% [!] @@VERIFY 로 표시된 곳: 사용자 .lts 의 실제 오브젝트 이름/명령. 나머지 매핑
%     로직(위치->zone->파일)은 그대로 사용 가능.

clear; global ID_LT ltml ltloc
RenewLightTools();
lt = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt, "ShowFileDialogBox", 0);

%% ===== 배열/소스 지오메트리 (사용자 값) =====
BASE      = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\';  % @@VERIFY .ent 가 있는 경로
arraySize = 15.0;    % [mm] 어레이 한 변
pitch     = 1.0;     % [mm] 렌즐릿 피치 @@VERIFY (모델 unit-cell 크기)
gap       = 1.4;     % [mm] 렌즈-기판 간격
Nzone     = 6;       % zone(=.ent) 개수 (generate_graded_offset_ents.m 와 일치)
entName   = @(k) sprintf('%sff_grad_z%02d.1.ent', BASE, k);   % zone k -> 파일
METHOD    = 'zones'; % 'zones' | 'entities'

%% ===== 렌즐릿 격자 위치 + zone 매핑 =====
c  = (arraySize - pitch)/2;
xs = -c:pitch:c;  ys = -c:pitch:c;
[GX, GY] = meshgrid(xs, ys);
Rg = hypot(GX, GY);
Rmax = max(Rg(:));
% zone 경계: 반경을 Nzone 등분 (필요시 조명반경 기준으로 교체)
zoneOf = min(floor(Rg / (Rmax/Nzone)), Nzone-1);   % 0..Nzone-1
fprintf('렌즐릿 %d개, zone 분포:', numel(GX));
for k=0:Nzone-1, fprintf(' z%d=%d', k, sum(zoneOf(:)==k)); end; fprintf('\n');

%% ===== 방법별 적용 =====
switch METHOD
case 'zones'
    % 모델에 zone 별 텍스처(예: 'zone0'..'zone{N-1}') 가 있고 각각 LibraryElement 를
    % 가진다고 가정. 각 텍스처의 unit-cell Filename 을 해당 zone .ent 로 지정.
    for k = 0:Nzone-1
        texName = sprintf('zone%d', k);          % @@VERIFY 모델 텍스처 이름 규약
        try
            List = ltml.LTDbList(lt,'LENS_MANAGER[1]','LIBRARY_ELEMENT_UNIT_CELL');
            Key  = ltml.LTListByName(lt, List, texName);   % @@VERIFY 이름으로 찾기
            ltml.LTDbSet(lt, Key, 'Filename', entName(k));
            fprintf('  [zones] %s <- %s\n', texName, entName(k));
        catch me
            fprintf('  [zones] %s 실패: %s\n', texName, me.message);
        end
    end
    % 각 텍스처 zone 의 공간범위(TEXTURE_ZONE_EXTENT)를 동심 링으로 맞추는 것은
    % 모델에서 설정(반경 링 마스크). @@VERIFY

case 'entities'
    % 각 렌즐릿을 개별 엔티티로 배치. 기준 엔티티 하나를 복제→위치·파일 지정.
    % @@VERIFY: 아래 복제/배치 명령은 사용자 모델의 엔티티 이름/명령에 맞춰야 함.
    baseEnt = 'FreeformLens_base';   % @@VERIFY 복제 원본 엔티티 이름
    idx = 0;
    for i = 1:numel(GX)
        k = zoneOf(i);
        name = sprintf('FF_%03d', i);
        try
            % 예시(개념): 원본 복제 -> 위치 이동 -> 이 엔티티의 .ent 지정
            ltml.LTCmd(lt, sprintf('\\O"LENS_MANAGER[1].COMPONENTS[Components].SOLID[%s]" Copy=', baseEnt)); %#ok<*NASGU> @@VERIFY
            % 위치: 소스 위 gap, 셀 중심 (GX,GY)
            % ltml.LTDbSet(lt, thisKey, 'X', GX(i)); ...
            % ltml.LTDbSet(lt, freeformKey, 'Filename', entName(k));
            idx = idx + 1;
        catch me
            fprintf('  [entities] %s 실패: %s\n', name, me.message);
        end
    end
    fprintf('  [entities] 배치 시도 %d/%d (명령 @@VERIFY 필요)\n', idx, numel(GX));
end

fprintf('\n구성 완료(가정 하). 시뮬 실행 전 텍스처/엔티티 배치를 GUI 에서 육안 확인 권장.\n');

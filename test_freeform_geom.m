% ============================================================
%  test_freeform_geom.m  —  지오메트리 파이프라인 단독 점검용
%
%  전체 BO 를 돌리기 전에, 아래 두 @@VERIFY 지점만 빠르게 반복 확인한다:
%    (A) freeform STL 생성 -> LightTools import -> RepairEntities -> SaveLibrary(.ent)
%        가 실제로 .ent 파일을 만드는가 (import 명령이 맞는가)
%    (B) far-field INTENSITY_MESH 의 유효 인덱스 범위(longitude x latitude) 확인
%
%  사용법:
%    1) 아래 LT_IMPORT_CMD 를 매크로 레코더로 얻은 실제 명령으로 채운다.
%       - LightTools 에서 Console/Macro 레코더 시작 -> File>Import 로 .stl 하나 수동
%         임포트 -> 레코더 정지 -> 기록된 명령을 복사.
%       - 그 명령의 파일경로 부분을 '%s' 로 바꿔 LT_IMPORT_CMD 에 넣는다.
%       - 임포트된 solid 이름(레코더의 DefaultSelect 줄)을 LT_IMPORT_SOLIDNAME 에.
%    2) 이 스크립트를 실행. "[OK] .ent 생성됨" 이 나오면 (A) 통과 -> 본 코드
%       BO_Freeform3D_asym.m 상단의 LT_IMPORT_CMD/LT_IMPORT_SOLIDNAME 에 그대로 복사.
%    3) 출력되는 mesh 유효 인덱스 범위를 BO 코드 상단 MESH_NLONG/MESH_NLAT 에 반영.
% ============================================================
clear;
global ID_swept ID_LT ltml ltloc

% ---- 여기만 채우세요 (@@VERIFY) ----
LT_IMPORT_CMD       = '';               % 예: 'ImportFile "%s" Format=STL'  ('%s'=STL경로)
LT_IMPORT_SOLIDNAME = 'freeform_lens';  % import 후 solid 이름
% ------------------------------------

base = 'C:\Users\jhkim\Desktop\Green_CE_Calculation\';

%% LightTools 연결
RenewLightTools();     % BO_Freeform3D_asym.m 의 것과 동일 함수. 경로는 그 안에서 설정.
lt  = ltloc.GetLTAPI(ID_swept);
lt2 = ltloc.GetLTAPI(ID_LT);
ltml.LTSetOption(lt,  "ShowFileDialogBox", 0);
ltml.LTSetOption(lt2, "ShowFileDialogBox", 0);

%% (A) 지오메트리 파이프라인 점검
% 고정 예시 4x4 격자(대칭 아님) 하나로 STL 생성
C = [0.20 0.35 0.30 0.15;
     0.35 0.70 0.55 0.25;
     0.30 0.55 0.45 0.20;
     0.15 0.25 0.20 0.10];
hfun = @(X,Y) freeform_grid_height(X, Y, C, 1.0);

tagc = char(datetime('now','Format','yyyyMMddHHmmssSSS'));
stlPath       = [base 'freeform_' tagc '.stl'];
entPath_local = [base 'swept_'    tagc '.ent'];

info = generate_freeform_mesh(hfun, stlPath, ...
    struct('nr',40,'nt',120,'Rfoot',1,'solidName',LT_IMPORT_SOLIDNAME));
fprintf('[Mesh] STL 저장: %s\n', stlPath);
fprintf('[Mesh] zmax=%.4f mm, 부피=%.4f, centroid=(%+.4f,%+.4f)\n', ...
    info.zmax, info.volume, info.centroidXY(1), info.centroidXY(2));

if isempty(LT_IMPORT_CMD)
    fprintf(2, ['\n[!] LT_IMPORT_CMD 가 비어 있습니다. STL 은 만들어졌으니, 지금 이 STL 을\n' ...
        '    LightTools GUI 에서 매크로 레코더로 File>Import 해 보고 기록된 명령을\n' ...
        '    LT_IMPORT_CMD 에 넣은 뒤 다시 실행하세요.\n    STL: %s\n'], stlPath);
    return;
end

if exist(entPath_local, 'file'), delete(entPath_local); end
ltml.LTCmd(lt, sprintf(LT_IMPORT_CMD, stlPath));
ltml.LTCmd(lt, sprintf('DefaultSelect "%s"', LT_IMPORT_SOLIDNAME));
ltml.LTCmd(lt, 'RepairEntities');
ltml.LTCmd(lt, sprintf('SaveLibrary XYZ 0,0,0 "%s"', entPath_local));
pause(0.5);

if exist(entPath_local, 'file')
    fprintf('\n[OK] .ent 생성됨: %s\n', entPath_local);
    fprintf('     -> import 파이프라인 통과. 이 LT_IMPORT_CMD/SOLIDNAME 을 본 코드에 복사하세요.\n');
else
    fprintf(2, ['\n[FAIL] .ent 미생성. LightTools 콘솔의 오류 메시지를 확인하세요.\n' ...
        '       (Unknown variable ... 이면 import 명령 문자열이 틀린 것)\n']);
end

%% (B) far-field mesh 유효 인덱스 범위 스캔
% CellValue_UI(iLong, iLat) 를 키워가며 반환값이 유효한 최대 인덱스를 찾는다.
% (LTCmd 와 달리 LTDbGet 는 잘못된 인덱스에서 콘솔 오류를 내지만 MATLAB 값은 보통 0.
%  따라서 여기서는 "연속으로 0 이 시작되는 지점" 을 대략적 상한 후보로 출력만 한다.
%  정확한 값은 GUI 의 mesh 셀 수와 대조할 것.)
try
    List = ltml.LTDbList(lt2,'lens_manager[1]','INTENSITY_MESH');
    Key  = ltml.LTListAtPos(lt2, List, 3);   % MESH_POS 와 동일하게
    fprintf('\n[Mesh probe] CellValue_UI(iLong, iLat) 스캔 (참고용):\n');
    maxProbe = 400;
    % longitude 방향 (iLat=1 고정)
    lastLong = 0;
    for k = 1:maxProbe
        v = ltml.LTDbGet(lt2, Key, 'CellValue_UI', k, 1);
        if isempty(v) || ~isfinite(v), break; end
        lastLong = k;
    end
    % latitude 방향 (iLong=1 고정)
    lastLat = 0;
    for k = 1:maxProbe
        v = ltml.LTDbGet(lt2, Key, 'CellValue_UI', 1, k);
        if isempty(v) || ~isfinite(v), break; end
        lastLat = k;
    end
    fprintf('  longitude(iLong) 최대 유효 근처 = %d,  latitude(iLat) 최대 유효 근처 = %d\n', ...
        lastLong, lastLat);
    fprintf('  -> 콘솔 "잘못된 인덱스(N,·)" 의 N-1 이 실제 longitude 셀 수. GUI 와 대조 후\n');
    fprintf('     BO_Freeform3D_asym.m 의 MESH_NLONG / MESH_NLAT 에 정확히 반영하세요.\n');
catch me
    fprintf(2, '[Mesh probe] 실패: %s\n', me.message);
end

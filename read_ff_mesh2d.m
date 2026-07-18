function [Igrid, thC, phC] = read_ff_mesh2d(ltml, lt, meshPos, cfg)
% READ_FF_MESH2D  far-field INTENSITY_MESH 를 2D [nLat x nLong] 로 읽는다.
%
%   [Igrid, thC, phC] = read_ff_mesh2d(ltml, lt, meshPos, cfg)
%     meshPos : INTENSITY_MESH 리스트 내 위치 (기존 코드의 far-field mesh = 3)
%     cfg     : .nLong .nLat .longMin .longMax .latMin .latMax
%               (export 데이터가 theta 2도 x phi 10도였으므로 기본 45x36 가정)
%               + (선택) .exportObjPath, .exportTmpDir  -> 아래 [속도] 참고
%     Igrid   : nLat x nLong  (행=theta, 열=phi)
%     thC,phC : 각 bin 중심각 [deg]
%
% [주의] CellValue_UI(iLong, iLat) : 첫 인자 longitude, 둘째 latitude.
%        mesh 는 Symmetry = "No Symmetry" 여야 phi 분해능이 살아있다. @@VERIFY
% [주의] 잘못된 인덱스 UI 에러가 나면 nLong/nLat 이 모델과 다른 것. 모델 mesh
%        속성창의 (longitude, latitude) bin 수를 그대로 넣을 것.
%
% [속도] 기본 경로는 셀 하나하나를 COM 으로 가져온다(nLong*nLat 회 왕복 ->
% 90x36 이면 3240회, 매 평가·매 파장마다 반복돼 매우 느리다). cfg.exportObjPath
% 를 지정하면 LightTools 콘솔 명령 ExportMeshToFile 로 mesh 전체를 파일 1회
% 저장 후 벡터화된 파일 읽기로 대체한다 (COM 호출 1회로 감소).
%   cfg.exportObjPath : ExportMeshToFile 에 넘길 mesh 의 전체 오브젝트 경로
%     예) 'ILLUM_MANAGER[Illumination Manager].RECEIVERS[Receiver List].FARFIELD_RECEIVER[farFieldReceiver_21].FORWARD_SIM_FUNCTION[Forward Simulation].INTENSITY_MESH[Intensity Mesh]'
%     @@VERIFY: 리시버 인스턴스 이름(farFieldReceiver_21)은 모델마다 다를 수 있음.
%   cfg.exportTmpDir  : 임시 export 파일을 쓸 디렉토리 (없으면 현재 폴더)
%
% [!] @@VERIFY 파일 포맷: parse_exported_mesh_file() 은 사용자가 이전에 공유한
%     xlsx export(첫 행=phi 헤더, 첫 열=theta, 나머지=intensity, tab 구분)와
%     동일하다고 가정한 초안이다. mesh2.1.txt 실제 내용 몇 줄을 확인해 다르면
%     그 함수만 고치면 된다(이 함수의 나머지 로직은 그대로 재사용 가능).

if nargin < 4 || isempty(cfg)
    cfg = struct('nLong',36,'nLat',45,'longMin',-180,'longMax',180, ...
                 'latMin',0,'latMax',90);
end

if isfield(cfg,'exportObjPath') && ~isempty(cfg.exportObjPath)
    % ---- 빠른 경로: ExportMeshToFile (COM 호출 1회) ----
    if isfield(cfg,'exportTmpDir') && ~isempty(cfg.exportTmpDir)
        tmpDir = cfg.exportTmpDir;
    else
        tmpDir = '';
    end
    rng('shuffle');
    charSet = ['a':'z' 'A':'Z' '0':'9'];
    tag = charSet(randi(numel(charSet), 1, 10));
    tmpFile = [tmpDir 'ffmesh_' tag '.txt'];   % 문자열 접합(sprintf 아님) -> 경로의 '\' 안전

    ltml.LTCmd(lt, ['ExportMeshToFile "' cfg.exportObjPath '" "' tmpFile '"']);
    Igrid = parse_exported_mesh_file(tmpFile, cfg.nLat, cfg.nLong);
    delete(tmpFile);
else
    % ---- 기존 경로: 셀 단위 COM 읽기 (폴백) ----
    List = ltml.LTDbList(lt,'lens_manager[1]','INTENSITY_MESH');
    Key  = ltml.LTListAtPos(lt,List,meshPos);
    Igrid = zeros(cfg.nLat, cfg.nLong);
    for iL = 1:cfg.nLong
        for iT = 1:cfg.nLat
            v = ltml.LTDbGet(lt, Key, 'CellValue_UI', iL, iT);
            if isempty(v) || ~isfinite(v), v = 0; end
            Igrid(iT, iL) = v;
        end
    end
end

thC = cfg.latMin  + (cfg.latMax -cfg.latMin) /cfg.nLat  * ((1:cfg.nLat)  - 0.5);
phC = cfg.longMin + (cfg.longMax-cfg.longMin)/cfg.nLong * ((1:cfg.nLong) - 0.5);
end


function Igrid = parse_exported_mesh_file(filePath, nLat, nLong)
% PARSE_EXPORTED_MESH_FILE  ExportMeshToFile 출력을 [nLat x nLong] 로 파싱.
%
% [!] @@VERIFY: 아래는 이전에 공유된 xlsx(45x36, 탭 구분, 1행=phi 헤더,
%     1열=theta) 와 동일 포맷이라는 가정의 초안이다. mesh2.1.txt 앞 3~5줄을
%     보내주면 실제 구분자/헤더줄수/행-열 방향을 맞춰 이 함수만 수정하면 된다.
raw = readmatrix(filePath, 'FileType','text');
Igrid = raw(2:end, 2:end);
if ~isequal(size(Igrid), [nLat, nLong])
    error('parse_exported_mesh_file:size', ...
        'export 파일 크기(%dx%d)가 기대(%dx%d)와 다름 - 포맷/헤더 가정을 확인할 것.', ...
        size(Igrid,1), size(Igrid,2), nLat, nLong);
end
end

function [Igrid, thC, phC] = read_ff_mesh2d(ltml, lt, meshPos, cfg)
% READ_FF_MESH2D  far-field INTENSITY_MESH 를 2D [nLat x nLong] 로 읽는다.
%
%   [Igrid, thC, phC] = read_ff_mesh2d(ltml, lt, meshPos, cfg)
%     meshPos : INTENSITY_MESH 리스트 내 위치 (기존 코드의 far-field mesh = 3)
%     cfg     : .nLong .nLat .longMin .longMax .latMin .latMax
%               (export 데이터가 theta 2도 x phi 10도였으므로 기본 45x36 가정)
%     Igrid   : nLat x nLong  (행=theta, 열=phi)
%     thC,phC : 각 bin 중심각 [deg]
%
% [주의] CellValue_UI(iLong, iLat) : 첫 인자 longitude, 둘째 latitude.
%        mesh 는 Symmetry = "No Symmetry" 여야 phi 분해능이 살아있다. @@VERIFY
% [주의] 잘못된 인덱스 UI 에러가 나면 nLong/nLat 이 모델과 다른 것. 모델 mesh
%        속성창의 (longitude, latitude) bin 수를 그대로 넣을 것.

if nargin < 4 || isempty(cfg)
    cfg = struct('nLong',36,'nLat',45,'longMin',-180,'longMax',180, ...
                 'latMin',0,'latMax',90);
end
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
thC = cfg.latMin  + (cfg.latMax -cfg.latMin) /cfg.nLat  * ((1:cfg.nLat)  - 0.5);
phC = cfg.longMin + (cfg.longMax-cfg.longMin)/cfg.nLong * ((1:cfg.nLong) - 0.5);
end

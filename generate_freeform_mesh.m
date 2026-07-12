function out = generate_freeform_mesh(H, pCtrl, harm, stlPath, opts)
% GENERATE_FREEFORM_MESH  비대칭 freeform 렌즈를 watertight 삼각 메쉬로 만들어
% STL(ASCII)로 저장한다. LightTools 로 import -> (모델에서) SaveLibrary 로 .ent(ACIS)
% 변환하는 파이프라인의 앞단이다.
%
%   out = generate_freeform_mesh(H, pCtrl, harm, stlPath, opts)
%
%   solid 구성 (plano-convex, 평평한 바닥):
%     - 윗면  : freeform cap  z = freeform_height(x,y,...)  (rim rho=1 에서 z=0)
%     - 아랫면: 평면 disk z=0 (같은 footprint)
%     - rim(rho=1)에서 윗/아랫면이 z=0 으로 만나 side wall 없이 닫힌 solid.
%       (P(1)=0 이므로 cap 이 자연히 바닥으로 내려와 닫힘)
%
%   입력
%     H, pCtrl, harm : freeform_height 파라미터
%     stlPath        : 저장할 .stl 경로(문자열). 비우면 저장 생략(메쉬만 반환).
%     opts           : (선택) 구조체
%                        .nr        반경 방향 링 수     (기본 40)
%                        .nt        방위각 분할 수      (기본 120)
%                        .Rfoot     footprint 물리 반경 [mm] (기본 1; 좌표 스케일)
%                        .solidName STL solid 이름      (기본 'freeform_lens')
%
%   출력 out 구조체
%     .V         Nx3 정점 [x y z] (mm)
%     .F         Mx3 삼각형 (1-based 인덱스)
%     .centroidXY  [cx cy]  높이가중 무게중심 (비대칭/tilt 방향 진단; 대칭이면 ~0)
%     .volume    근사 부피
%     .zmax      최대 높이
%     .stlPath   저장 경로(저장 시)
%
%   주의: 좌표는 opts.Rfoot 로 스케일된 물리 단위(mm). 기존 파이프라인이 프로파일을
%   max<=1 로 정규화하던 것과 달리 여기서는 footprint 반경을 명시적으로 준다.

    if nargin < 5 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'nr'),        opts.nr = 40;              end
    if ~isfield(opts, 'nt'),        opts.nt = 120;             end
    if ~isfield(opts, 'Rfoot'),     opts.Rfoot = 1;            end
    if ~isfield(opts, 'solidName'), opts.solidName = 'freeform_lens'; end

    nr = opts.nr;  nt = opts.nt;  R = opts.Rfoot;

    % --- 극좌표 그리드 (정점) : 정점(중심) + nr 개 링 x nt 방위각 ---
    rings = linspace(0, 1, nr + 1);         % 0 .. 1 (정규 반경)
    rings = rings(2:end);                    % 링(중심 제외): nr 개
    ang   = linspace(0, 2*pi, nt + 1);
    ang   = ang(1:end-1);                    % nt 개(중복 제거)

    % 정규 좌표에서 높이 계산 후 물리 좌표(R 배)로 확대
    % top 정점
    topApex = [0, 0, freeform_height(0, 0, H, pCtrl, harm)];  % (0,0,H)
    Vtop = zeros(nr*nt, 3);
    for i = 1:nr
        r = rings(i);
        xr = r*cos(ang);  yr = r*sin(ang);
        zr = freeform_height(xr, yr, H, pCtrl, harm);
        zr(~isfinite(zr)) = 0;
        idx = (i-1)*nt + (1:nt);
        Vtop(idx, :) = [R*xr(:), R*yr(:), zr(:)];
    end

    % bottom 정점 (z=0). rim 링(i=nr, z=0)은 top 과 좌표가 동일하므로 공유한다.
    botCenter = [0, 0, 0];
    Vbot = zeros((nr-1)*nt, 3);              % rim 제외한 내부 링만 (마지막 링 공유)
    for i = 1:nr-1
        r = rings(i);
        xr = r*cos(ang);  yr = r*sin(ang);
        idx = (i-1)*nt + (1:nt);
        Vbot(idx, :) = [R*xr(:), R*yr(:), zeros(nt,1)];
    end

    % --- 정점 배열 조립 & 인덱스 헬퍼 ---
    % 순서: [topApex ; Vtop(ring1..nr) ; botCenter ; Vbot(ring1..nr-1)]
    V = [topApex; Vtop; botCenter; Vbot];
    iApex = 1;
    topRing = @(i,j) 1 + (i-1)*nt + mod(j-1, nt) + 1;      % top 링 i(1..nr), 방위 j
    iBotCenter = 1 + nr*nt + 1;
    botRing = @(i,j) iBotCenter + (i-1)*nt + mod(j-1, nt) + 1; % bot 링 i(1..nr-1) (+1: 중심 정점 다음부터)
    rimRing = @(j) topRing(nr, j);                         % rim = top 마지막 링(z=0) 공유

    F = zeros(0, 3);

    % --- 윗면(cap) 삼각형 ---
    % 정점 팬 (apex -> 링1)
    for j = 1:nt
        F(end+1, :) = [iApex, topRing(1,j), topRing(1,j+1)]; %#ok<AGROW>
    end
    % 링 사이 quad(2 삼각형)
    for i = 1:nr-1
        for j = 1:nt
            a = topRing(i,   j);   b = topRing(i,   j+1);
            c = topRing(i+1, j);   d = topRing(i+1, j+1);
            F(end+1, :) = [a, c, d]; %#ok<AGROW>
            F(end+1, :) = [a, d, b]; %#ok<AGROW>
        end
    end

    % --- 아랫면(flat) 삼각형 (법선 -Z, 감김 반대로) ---
    for j = 1:nt
        F(end+1, :) = [iBotCenter, botRing(1,j+1), botRing(1,j)]; %#ok<AGROW>
    end
    for i = 1:nr-2
        for j = 1:nt
            a = botRing(i,   j);   b = botRing(i,   j+1);
            c = botRing(i+1, j);   d = botRing(i+1, j+1);
            F(end+1, :) = [a, d, c]; %#ok<AGROW>
            F(end+1, :) = [a, b, d]; %#ok<AGROW>
        end
    end
    % 마지막 내부 bot 링(nr-1) -> rim(top 마지막 링, 공유) 연결
    for j = 1:nt
        a = botRing(nr-1, j);   b = botRing(nr-1, j+1);
        c = rimRing(j);         d = rimRing(j+1);
        F(end+1, :) = [a, d, c]; %#ok<AGROW>
        F(end+1, :) = [a, b, d]; %#ok<AGROW>
    end

    % --- 진단값 ---
    zt = Vtop(:,3);
    wsum = sum(zt);
    if wsum > 0
        cx = sum(Vtop(:,1).*zt)/wsum;
        cy = sum(Vtop(:,2).*zt)/wsum;
    else
        cx = 0; cy = 0;
    end
    out.V = V;  out.F = F;
    out.centroidXY = [cx, cy];
    out.zmax = max(zt);
    out.volume = tri_mesh_volume(V, F);

    % --- STL(ASCII) 저장 ---
    if nargin >= 4 && ~isempty(stlPath)
        write_ascii_stl(stlPath, V, F, opts.solidName);
        out.stlPath = stlPath;
    end
end

% ------------------------------------------------------------------------
function write_ascii_stl(path, V, F, name)
    fid = fopen(path, 'w');
    if fid < 0, error('generate_freeform_mesh:stl', 'STL 파일 열기 실패: %s', path); end
    fprintf(fid, 'solid %s\n', name);
    for k = 1:size(F,1)
        p1 = V(F(k,1),:);  p2 = V(F(k,2),:);  p3 = V(F(k,3),:);
        n  = cross(p2-p1, p3-p1);
        nn = norm(n);
        if nn > 0, n = n/nn; else, n = [0 0 0]; end
        fprintf(fid, '  facet normal %.7e %.7e %.7e\n', n(1), n(2), n(3));
        fprintf(fid, '    outer loop\n');
        fprintf(fid, '      vertex %.7e %.7e %.7e\n', p1(1), p1(2), p1(3));
        fprintf(fid, '      vertex %.7e %.7e %.7e\n', p2(1), p2(2), p2(3));
        fprintf(fid, '      vertex %.7e %.7e %.7e\n', p3(1), p3(2), p3(3));
        fprintf(fid, '    endloop\n  endfacet\n');
    end
    fprintf(fid, 'endsolid %s\n', name);
    fclose(fid);   % LightTools 가 읽기 전 반드시 flush/unlock
end

% ------------------------------------------------------------------------
function vol = tri_mesh_volume(V, F)
% 닫힌 삼각 메쉬의 부호있는 부피(발산정리). watertight 이면 양수여야 함.
    v1 = V(F(:,1),:);  v2 = V(F(:,2),:);  v3 = V(F(:,3),:);
    vol = sum(dot(v1, cross(v2, v3, 2), 2)) / 6;
    vol = abs(vol);
end

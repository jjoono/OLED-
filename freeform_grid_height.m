function Z = freeform_grid_height(X, Y, C, winPow)
% FREEFORM_GRID_HEIGHT  2D B-spline 제어격자 기반 "임의 3D" freeform 높이장 z(x,y)
%
%   Z = freeform_grid_height(X, Y, C, winPow)
%
%   디스크(반경 1) 위에 정의된 NgxNg 제어높이 격자 C 를 bicubic spline(interp2)으로
%   보간하고, rim(rho=1)에서 0 이 되도록 윈도우 (1-rho^2)^winPow 를 곱한다.
%
%       z(x,y) = max( interp2_spline(C)(x,y), 0 ) * (1 - rho^2)^winPow ,   rho<=1
%       z      = NaN,                                                      rho> 1
%
%   [왜 격자인가] 회전대칭 프로파일 + 저차 harmonic(정점 tilt) 은 "약간 틀어진"
%   비대칭만 만든다. NgxNg 제어높이는 국소 bump/dimple 까지 포함한 임의의 3D 형상을
%   표현한다(높이장이므로 single-valued -> 기하학적으로 항상 유효). Ng=4 -> 16 DOF.
%
%   [rim 윈도우] winPow>=1 이면 rim 기울기가 유한 -> 제조/레이트레이싱/메쉬 품질 안정.
%   winPow=1 이 기본. rho=1 에서 정확히 0 이므로 generate_freeform_mesh 의 rim ring 이
%   바닥과 만나 watertight solid 가 유지된다.
%
%   입력
%     X, Y   : 같은 크기 좌표 배열. footprint 반경 1 로 정규화된 좌표.
%     C      : Ng x Ng 제어높이 [mm]. C(i,j) = 격자점 (gx(j), gy(i)) 의 높이.
%     winPow : (선택) rim 윈도우 지수. 기본 1.0.
%
%   출력
%     Z      : X,Y 와 같은 크기. rho>1 은 NaN.
%
%   * generate_freeform_mesh / isValidFreeformGrid / BO 목적함수가 공유한다.
%     동일 수식을 (검증용) Python 프로토타입과 일치시켜 유지할 것.

    if nargin < 4 || isempty(winPow), winPow = 1.0; end

    Ng = size(C, 1);
    gx = linspace(-1, 1, Ng);
    [GX, GY] = meshgrid(gx, gx);

    Zint = interp2(GX, GY, C, X, Y, 'spline');   % 전역 bicubic spline (질의점은 [-1,1] 내)

    r2  = X.^2 + Y.^2;
    win = max(1 - r2, 0) .^ winPow;              % rim(rho=1) -> 0

    Z = max(Zint, 0) .* win;                     % 음수 높이 방지 후 윈도우
    Z(r2 > 1) = NaN;                             % footprint 밖
end

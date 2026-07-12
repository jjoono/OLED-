function Z = freeform_height(X, Y, H, pCtrl, harm)
% FREEFORM_HEIGHT  비대칭 freeform MLA 렌즈의 높이장 z(x,y)
%
%   Z = freeform_height(X, Y, H, pCtrl, harm)
%
%   단일 렌즈의 표면 높이를 극좌표(rho, phi)로 정의한다:
%
%       z(rho,phi) = H * P(rho) * S(rho,phi)
%
%     - P(rho)      : 반경 방향 base 프로파일 (정점 P(0)=1, rim P(1)=0).
%                     제어높이 pCtrl 를 PCHIP(단조 보존) 보간.
%     - S(rho,phi)  : 방위각 비대칭 인자
%                       S = 1 + sum_m rho^m ( c_m cos(m phi) + s_m sin(m phi) )
%                     m=1 항(c1,s1)이 "정점 tilt" = 방향성 비대칭의 핵심 DOF.
%                     rho^m 가중 -> 정점 근처는 대칭 유지, rim으로 갈수록 비대칭
%                     심화 -> 렌즈 footprint(=rim)는 원형으로 유지되어 배열/제조 용이.
%
%   [핵심] 모든 개별 렌즈가 "동일한 비대칭 형상"을 갖고 같은 방향으로 정렬되면
%   array 전체가 방향성 비대칭 발광을 만든다. hemisphere는 rotational symmetry가
%   기하학적으로 강제되어(어떤 파라미터로도 S=1로 고정) 이 자유도가 원리적으로 없다.
%
%   입력
%     X, Y  : 같은 크기의 좌표 배열. 단위 footprint(반경 1)로 정규화된 좌표.
%     H     : 렌즈 전체 높이(aspect) 스칼라.
%     pCtrl : 1xK 반경 base 제어높이 (rho = 1/(K+1) .. K/(K+1) 위치). 정점/ rim 제외.
%     harm  : Mx3 행렬, 각 행 [m, c_m, s_m].  비면 회전대칭.
%
%   출력
%     Z     : X,Y 와 같은 크기. rho>1 (footprint 밖)은 NaN.
%
%   * 이 함수는 generate_freeform_mesh / isValidFreeform / BO 목적함수가 공유한다.
%     MATLAB 과 (검증용) Python 프로토타입이 동일 수식을 계산하도록 유지할 것.

    rho = hypot(X, Y);
    phi = atan2(Y, X);

    P = radial_profile(rho, pCtrl);
    S = asym_factor(rho, phi, harm);

    Z = H .* P .* S;
    Z(rho > 1) = NaN;
end

% ------------------------------------------------------------------------
function P = radial_profile(rho, pCtrl)
% 정점(rho=0)=1, rim(rho=1)=0 고정 + 내부 제어높이 pCtrl 를 PCHIP 보간.
    K  = numel(pCtrl);
    rc = linspace(0, 1, K + 2);              % 0, 내부 K개, 1
    hc = [1, pCtrl(:).', 0];
    P  = pchip(rc, hc, min(max(rho, 0), 1)); % 단조 보존 -> 물결/오버슛 억제
    P  = max(P, 0);                          % 음수 높이 방지
end

% ------------------------------------------------------------------------
function S = asym_factor(rho, phi, harm)
% 방위각 비대칭 인자. harm 각 행 [m c_m s_m]. rho^m 가중.
    S = ones(size(rho));
    for k = 1:size(harm, 1)
        m  = harm(k, 1);
        cm = harm(k, 2);
        sm = harm(k, 3);
        S  = S + (rho.^m) .* (cm .* cos(m .* phi) + sm .* sin(m .* phi));
    end
    S = max(S, 0.05);   % 양수 유지(단일값 표면/제조성). 위반 여부는 isValidFreeform에서 별도 판정.
end

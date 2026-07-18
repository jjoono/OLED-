function R = detect_phi_window(Wacc, thC, phC, thLo, thHi, phiWidth)
% DETECT_PHI_WINDOW  theta 밴드는 고정, phi 창 위치는 자동 검출.
%
%   R = detect_phi_window(Wacc, thC, phC, thLo, thHi, phiWidth)
%     Wacc     : nLat x nLong 파워가중 그리드 (= I .* sin(theta), 파장가중 누적본)
%     thLo,thHi: 응용에 따라 고정하는 theta 밴드 [deg] (예: 고각 40~60)
%     phiWidth : phi 창 전체폭 [deg] (예: 80 -> 중심 ±40)
%
%   출력 R:
%     .phiC       검출된 phi 창 중심 [deg]  (빛이 실제로 몰린 방향)
%     .fracWin    창내 파워 / 전체 파워      (EQE_region 분율)
%     .PWin       창내 파워 (Wacc 합, 절대값)
%     .Ptot       전체 파워
%     .contrast   창내 평균 / 반대편(phiC+180 중심 동일폭) 평균  (phi-대비비)
%     .fracIso    phi 균일 가정시 같은 창의 분율 (= 밴드분율 x phiWidth/360)
%                 -> fracWin/fracIso > 1 이면 phi 집중이 실재
%
% [원리] 어느 phi 로 몰릴지 모르므로, 고정폭 창을 wrap-around 로 슬라이드하며
% 창내 파워가 최대가 되는 중심을 찾는다 (원형 이동합 = 최적 창 위치).

tm = (thC >= thLo) & (thC <= thHi);
band = sum(Wacc(tm, :), 1);          % 1 x nLong : theta 밴드 내 phi 분포
Ptot = sum(Wacc(:));
Pband = sum(band);

% 원형(wrap) 이동합: 각 bin 을 중심으로 phiWidth 창의 파워
nL   = numel(phC);
dphi = abs(phC(2) - phC(1));
half = phiWidth/2;
winP = zeros(1, nL);
for c = 1:nL
    d = abs(mod(phC - phC(c) + 180, 360) - 180);   % wrap 각거리
    winP(c) = sum(band(d <= half));
end
[PWin, ic] = max(winP);
phiC = phC(ic);

% 반대편 동일폭 창 (대비비)
dOpp = abs(mod(phC - (phiC+180) + 180, 360) - 180);
POpp = sum(band(dOpp <= half));
nWin = sum(abs(mod(phC - phiC + 180, 360) - 180) <= half);   % 창 bin 수(동일)

R = struct();
R.phiC     = phiC;
R.PWin     = PWin;
R.Ptot     = Ptot;
R.fracWin  = PWin / max(Ptot, eps);
R.contrast = (PWin/max(nWin,1)) / max(POpp/max(nWin,1), eps);
R.fracIso  = (Pband/max(Ptot,eps)) * (phiWidth/360);   % phi 균일이면 이 값
end

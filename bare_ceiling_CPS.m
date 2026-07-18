% ============================================================
%  bare_ceiling_CPS.m  -- 렌즈 없이 "CPS만"으로 EQE_region 상한 (LightTools 불필요)
%
%  [원리] 렌즈가 없으면 substrate(n=1.51)->air 는 "평평한 Fresnel 계면"이다.
%    -> ray tracing 없이 해석적으로 계산 가능:
%       CPS 로 I_sub(θ_sub, λ) 계산 -> 평면계면 Snell+Fresnel 로 air 로 투과
%       -> 목표 (θ_air 밴드 x φ 창) 로 나가는 EQE 를 적분.
%    목표 밴드 θ_air∈[40,60] 는 θ_sub∈[25.2°,35.0°](임계각 41.5° 이내)에서만
%    오므로 전부 "투과광" -> TIR/재활용 무관, 해석식이 정확.
%    bare 는 φ-대칭이라 φ 창은 (창폭/360) 만 곱하면 됨.
%
%  [변수] dETL,dHTL,dAg  (Ag 스택 = v4 와 동일). 렌즈 변수 없음.
%  [출력] EQE_region(상한) + 최적 캐비티 + θ_air 분포. 빠르므로 surrogateopt.
%
%  [정합성] LightTools 목적함수와 동일 정의:
%    EQE_region = (φ창/360) * Σ_λ EQE_sub_matrix(λ) * w_air_band(λ)/w_sub(λ)
%    (w_air = I_sub 를 Fresnel 투과시킨 air 파워, w_sub = Σ I_sub sinθ = weight_factor)
% ============================================================
clear;
global TH_LO TH_HI PHI_W N_SUB
TH_LO = 40;  TH_HI = 60;  PHI_W = 40;   % 목표 밴드/창 (v4 와 동일)
N_SUB = 1.51;                            % 기판 굴절률 (스택 마지막층)

lb = [ 10,  10,  0];    % dETL, dHTL, dAg
ub = [150, 150, 50];
varNames = {'dETL','dHTL','dAg'};

opts = optimoptions('surrogateopt','MaxFunctionEvaluations',300, ...
    'UseParallel',false,'PlotFcn',[],'Display','iter');
fprintf('bare ceiling (CPS only): θ_air∈[%d,%d], φ창=%d°, LightTools 불필요\n', TH_LO,TH_HI,PHI_W);
[xBest,fBest] = surrogateopt(@(x) -eqe_region_cps(x(1),x(2),x(3)), lb, ub, opts);

[Ereg, Etot, Eth, dist] = eqe_region_cps(xBest(1),xBest(2),xBest(3));
save('bare_ceiling_CPS_result.mat','xBest','Ereg','Etot','Eth','dist','TH_LO','TH_HI','PHI_W');

fprintf('\n######## Bare ceiling (CPS only) ########\n');
fprintf('  EQE_region(상한) = %.5g   [= θ밴드 EQE %.5g x φ창/360(%.3f)]\n', Ereg, Eth, PHI_W/360);
fprintf('  EQE_total(전방위 air) = %.5g\n', Etot);
fprintf('  최적 캐비티: dETL=%.1f  dHTL=%.1f  dAg=%.1f\n', xBest(1),xBest(2),xBest(3));
fprintf('  --> v4 렌즈 EQE_region 이 %.5g 를 넘는 만큼이 순수 φ-fold 렌즈 이득.\n', Ereg);

% θ_air 분포 플롯
figure('Name','bare 최적 캐비티 θ_air 분포','Color','w');
plot(dist.th_air, dist.I_air, 'LineWidth',1.8); grid on; hold on;
xline(TH_LO,'--'); xline(TH_HI,'--');
xlabel('\theta_{air} (deg)'); ylabel('air-side intensity (a.u.)');
title(sprintf('bare 최적 (dETL=%.0f,dHTL=%.0f,dAg=%.0f) | 회색선=목표밴드',xBest(1),xBest(2),xBest(3)));


%% =====================================================================
function [EQE_region, EQE_total, EQE_theta, dist] = eqe_region_cps(dETL, dHTL, dAg)
% CPS -> 평면 substrate->air Fresnel -> 목표 밴드 EQE
global TH_LO TH_HI PHI_W N_SUB
% nk_JH33.mat 에 material, spectrum 구조체가 들어있음(원본 코드와 동일 가정).
% 매 CPS 마다 재로드해도 되지만, 반복 최적화라 persistent 로 1회만.
persistent material spectrum
if isempty(material)
    D=load('nk_JH33.mat');  material=D.material;  spectrum=D.spectrum;
end

wavelength_start=580; wavelength_end=590;
wavelength=(wavelength_start:wavelength_end).';  wavelength_num=length(wavelength);
emission_spectrum=spectrum.l_I_Irdmppyph2tmd(wavelength_start-399:wavelength_end-399,:);
eta_rad=0.98; horizontal_dipole_ratio=0.865;
bottom_air_refractive_index=ones(wavelength_num,1);

no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_Ag_McPeak 1.51*ones(401,1)];
ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_Ag_McPeak 1.51*ones(401,1)];
no_bar=no_bar(wavelength_start-399:wavelength_end-399,:);
ne_bar=ne_bar(wavelength_start-399:wavelength_end-399,:);
thickness=[100 dETL 25 10 dHTL dAg];

CPS_result=CPS_for_Isub(no_bar,ne_bar,thickness,emission_spectrum,eta_rad, ...
    horizontal_dipole_ratio,bottom_air_refractive_index,4,12.5,499,3,wavelength);
EQE_sub_CPS=CPS_result.EQE_sub;
I_sub_s=CPS_result.I_sub_s;  I_sub_p=CPS_result.I_sub_p;   % [wavelength_num x 90], θ_sub=0..89
EQE_sub_matrix=CPS_result.EQE_sub_matrix;                  % [wavelength_num x 1] (or per-λ)
EQE_sub_matrix=EQE_sub_matrix(:)/sum(EQE_sub_matrix)*EQE_sub_CPS;

% --- 평면 substrate(N_SUB) -> air(1) Fresnel + Snell ---
th_sub = (0:89);  n1=N_SUB; n2=1;
s1=sind(th_sub);  sin2=n1/n2*s1;                     % Snell
esc = sin2<1;                                        % 임계각 이내만 투과
c1=cosd(th_sub);  c2=sqrt(max(1-sin2.^2,0));
Rs=((n1*c1-n2*c2)./(n1*c1+n2*c2)).^2;  Ts=1-Rs;
Rp=((n1*c2-n2*c1)./(n1*c2+n2*c1)).^2;  Tp=1-Rp;
Ts(~esc)=0;  Tp(~esc)=0;
th_air = nan(1,90);  th_air(esc)=asind(sin2(esc));   % 대응 air 각

bandAir = (th_air>=TH_LO) & (th_air<=TH_HI);          % 목표 θ_air 밴드
sinsub = sind(th_sub);

EQE_theta=0; EQE_total=0; I_air_accum=zeros(1,90);
for wi=1:wavelength_num
    Is=I_sub_s(wi,:); Ip=I_sub_p(wi,:);
    w_sub = sum( 0.5*(Is+Ip).*sinsub );                       % = weight_factor(λ)
    trans = Is.*Ts + Ip.*Tp;                                  % air 로 투과된 강도
    w_air_total = sum( trans.*sinsub );
    w_air_band  = sum( trans(bandAir).*sinsub(bandAir) );
    if w_sub>0
        EQE_theta = EQE_theta + EQE_sub_matrix(wi) * w_air_band /w_sub;
        EQE_total = EQE_total + EQE_sub_matrix(wi) * w_air_total/w_sub;
    end
    I_air_accum = I_air_accum + EQE_sub_matrix(wi)*trans;      % 분포(가중)
end
EQE_region = (PHI_W/360) * EQE_theta;   % bare φ-대칭 -> 창폭/360

% 플롯용 분포 (air 각 순으로 정렬)
[tha_sorted, idx] = sort(th_air);
dist.th_air = tha_sorted(~isnan(tha_sorted));
Ia = I_air_accum(idx);  dist.I_air = Ia(~isnan(tha_sorted));
end

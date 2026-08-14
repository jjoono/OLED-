%% 특정 두께(dETL, dHTL, dCap)에 대한 OLED 효율 계산
% PSO 최적화를 제거하고, 두께를 직접 지정해서 EQE / 모드 분포 / 각도별 효율을
% 확인하는 스크립트입니다.
%
% 사용법
%   - 아래 dETL_list / dHTL_list / dCap_list 에 값을 넣습니다.
%   - 스칼라 하나면 단일 조건 계산 (상세 결과를 화면에 출력)
%   - 벡터를 넣으면 모든 조합을 스윕해서 results 테이블로 정리
%       예) dETL_list = 30:5:80;  dHTL_list = 150:10:220;  dCap_list = 60;

clear; clc;

%% ===== 1. 계산할 두께 [nm] =====
dETL_list = 51;     % ETL (B3PyMPM) 두께
dHTL_list = 197;    % HTL (TAPC) 두께
dCap_list = 65;     % Capping layer 두께

%% ===== 2. 계산 파장 [nm] =====
% 원본 n/k, 스펙트럼 데이터는 400~800 nm, 1 nm 간격 401 포인트로 저장되어 있습니다.
% 전체 대역(400~800)이 물리적으로 올바른 EQE를 줍니다.
% 빠르게 경향만 볼 때는 595~600 처럼 좁혀도 되지만, EQE/CE/FWHM 절대값은 왜곡됩니다.
LAMBDA_MIN = 400;
LAMBDA_MAX = 800;

%% ===== 3. 스윕 실행 =====
[G_ETL, G_HTL, G_CAP] = ndgrid(dETL_list, dHTL_list, dCap_list);
n_case = numel(G_ETL);

varnames = {'dETL','dHTL','dCap','EQE_air','EQE_sub','EQE_sub_confined', ...
    'EQE_wg','EQE_spp','EQE_abs','EQE_0_20','EQE_20_40','EQE_40_60','EQE_60_80', ...
    'CE','FWHM','EQE_factor_air','Purcell_mean','obj_4060'};
R = zeros(n_case, numel(varnames));

fprintf('총 %d 개 조건 계산 (파장 %d ~ %d nm)\n', n_case, LAMBDA_MIN, LAMBDA_MAX);
t_all = tic;

for c = 1:n_case

    res = calc_OLED_efficiency(G_ETL(c), G_HTL(c), G_CAP(c), LAMBDA_MIN, LAMBDA_MAX);

    R(c,:) = [G_ETL(c), G_HTL(c), G_CAP(c), ...
        res.EQE_air, res.EQE_sub, res.EQE_sub_confined, ...
        res.EQE_wg, res.EQE_spp, res.EQE_abs, ...
        res.EQE_0_20, res.EQE_20_40, res.EQE_40_60, res.EQE_60_80, ...
        res.CE, res.FWHM, res.EQE_factor_air, mean(res.Purcell_factor), res.obj_4060];

    if n_case > 1 && mod(c, max(1,round(n_case/20))) == 0
        fprintf('  진행률 %5.1f%% (%d/%d)\n', 100*c/n_case, c, n_case);
    end

end

fprintf('계산 완료 (%.1f 초)\n\n', toc(t_all));

results = array2table(R, 'VariableNames', varnames);

%% ===== 4. 결과 출력 =====
if n_case == 1

    fprintf('====== 두께 조건 ======\n');
    fprintf('  dETL = %.2f nm\n', G_ETL(1));
    fprintf('  dHTL = %.2f nm\n', G_HTL(1));
    fprintf('  dCap = %.2f nm\n', G_CAP(1));
    fprintf('  전체 층 구조 [nm] : %s\n', mat2str(res.thickness));
    fprintf('\n====== 효율 (외부 방출 기준) ======\n');
    fprintf('  EQE (air, out-coupled) : %8.4f %%\n', 100*res.EQE_air);
    fprintf('  EQE (substrate)        : %8.4f %%\n', 100*res.EQE_sub);
    fprintf('  Current efficacy (CE)  : %8.4f cd/A  (per unit IQE)\n', res.CE);
    fprintf('  정면 스펙트럼 FWHM     : %8.1f nm\n', res.FWHM);
    fprintf('\n====== 광학 모드 분포 (power fraction) ======\n');
    fprintf('  Air (out-coupled)      : %8.4f %%\n', 100*res.Power_ratio_air);
    fprintf('  Substrate confined     : %8.4f %%\n', 100*res.Power_ratio_sub_confined);
    fprintf('  Waveguide              : %8.4f %%\n', 100*res.Power_ratio_wg);
    fprintf('  SPP                    : %8.4f %%\n', 100*res.Power_ratio_spp);
    fprintf('  Absorption             : %8.4f %%\n', 100*res.Power_ratio_abs);
    fprintf('\n====== 각도 구간별 EQE 기여 ======\n');
    fprintf('   0~20 deg : %8.4f %%\n', 100*res.EQE_0_20);
    fprintf('  20~40 deg : %8.4f %%\n', 100*res.EQE_20_40);
    fprintf('  40~60 deg : %8.4f %%\n', 100*res.EQE_40_60);
    fprintf('  60~80 deg : %8.4f %%\n', 100*res.EQE_60_80);
    fprintf('\n  40~60 deg 비중 x EQE (기존 목적함수) : %.6f\n', -res.obj_4060);
    fprintf('  평균 Purcell factor : %.4f\n', mean(res.Purcell_factor));
    fprintf('=====================================\n');

    % 각도 분포 / 정면 스펙트럼 확인용 그래프
    figure('Name','Angular & spectral profile');
    subplot(1,2,1);
    plot(0:89, res.I_air_total, 'LineWidth', 1.5); grid on;
    xlabel('Viewing angle (deg)'); ylabel('Normalized intensity');
    title('Angular emission profile (air)');
    subplot(1,2,2);
    plot(res.wavelength, res.I_air_front/max(res.I_air_front), 'LineWidth', 1.5); grid on;
    xlabel('Wavelength (nm)'); ylabel('Normalized intensity');
    title('Front (0 deg) spectrum');

else

    disp(results);

    [best_eqe, i_best] = max(R(:, strcmp(varnames,'EQE_air')));
    fprintf('\n최대 EQE(air) 조건: dETL = %.2f, dHTL = %.2f, dCap = %.2f  ->  EQE = %.4f %%\n', ...
        R(i_best,1), R(i_best,2), R(i_best,3), 100*best_eqe);

end

% 필요하면 저장
% writetable(results, 'EQE_sweep_results.csv');
% save('EQE_sweep_results.mat', 'results');


%% ================= 효율 계산 함수 =================
function res = calc_OLED_efficiency(dETL, dHTL, dCap, lambda_min, lambda_max)

if nargin < 4 || isempty(lambda_min), lambda_min = 400; end
if nargin < 5 || isempty(lambda_max), lambda_max = 800; end

% mat 파일은 한 번만 로드
persistent material spectrum V_401
if isempty(material)
    S_nk       = load('nk_JH_total.mat');       % material, spectrum 포함
    material   = S_nk.material;
    spectrum   = S_nk.spectrum;
    S_photopic = load('Photopic_400_800.mat');  % V_401 포함
    V_401      = S_photopic.V_401;
end

% Copyright ⓒ All Rights Reserved.

%% ===== 파장 및 데이터 인덱스 =====
% 원본 데이터는 400 nm 가 1행 -> lambda 의 행 인덱스 = lambda - 399
idx = (lambda_min - 399):(lambda_max - 399);

wavelength     = (lambda_min:lambda_max).';
wavelength_num = length(wavelength);

V_sel = V_401(idx, 1);   % 시감도 (CE 계산용)

%% ===== 발광 스펙트럼 =====
emission_spectrum = spectrum.l_I_Irdmppyph2tmd; %% JOSong, eta_rad=0.98, hdr=0.865
emission_spectrum = emission_spectrum(idx, 1);
emission_spectrum = emission_spectrum / sum(emission_spectrum);

eta_rad = 1;
horizontal_dipole_ratio = 0.95;
bottom_air_refractive_index = ones(wavelength_num, 1);

%% JO 구조 %%
no_bar = [ones(wavelength_num,1) material.l_Ag_McPeak(idx,1) material.l_B3_o_JO(idx,1) material.l_TCTA_B3_o_JO(idx,1) material.l_TCTA_o_JO(idx,1) material.l_TAPC_o_JO(idx,1) material.l_Ag_McPeak(idx,1) 2.3*ones(wavelength_num,1) ones(wavelength_num,1)];
ne_bar = [ones(wavelength_num,1) material.l_Ag_McPeak(idx,1) material.l_B3_e_JO(idx,1) material.l_TCTA_B3_e_JO(idx,1) material.l_TCTA_e_JO(idx,1) material.l_TAPC_e_JO(idx,1) material.l_Ag_McPeak(idx,1) 2.3*ones(wavelength_num,1) ones(wavelength_num,1)];

%% ===== 층 구조 =====
% [Ag(cathode) ETL EML(TCTA:B3) TCTA HTL(TAPC) Ag(anode, semi-transparent) Cap]
thickness = [100 dETL 25 10 dHTL 12 dCap];

EML_position = 4;

z0 = 12.5;          % emitting position in EML
u_data_num = 485;
max_u = 3;

%%

sin089 = sind(0:89);

layer_num = size(no_bar,2);

u = [(0:u_data_num-1)/u_data_num (u_data_num+1:u_data_num*max_u)/u_data_num];

u_num = length(u);

TMF_bottom = TMF_birefringence_whole(no_bar(:,EML_position:layer_num),ne_bar(:,EML_position:layer_num),[thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0],u,wavelength);
TMF_top    = TMF_birefringence_whole(no_bar(:,EML_position:-1:1),ne_bar(:,EML_position:-1:1),[z0 thickness(EML_position-2:-1:1) 0],u,wavelength);

K_p_v = 3/4*real(ne_bar(:,EML_position)./no_bar(:,EML_position)*(u.^2./sqrt(1-u.^2)).*(1+TMF_bottom.r_p).*(1+TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
K_p_h = real(3./(6*(no_bar(:,EML_position)./ne_bar(:,EML_position)).^2+2)*sqrt(1-u.^2).*(1-TMF_bottom.r_p).*(1-TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
K_s_h = real(3./(2*(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2+6)*(1./sqrt(1-u.^2)).*(1+TMF_bottom.r_s).*(1+TMF_top.r_s)./(1-TMF_bottom.r_s.*TMF_top.r_s));

K_p_v2 = K_p_v;
K_p_h2 = K_p_h;
K_s_h2 = K_s_h;

K_p_v2_total = K_p_v;
K_p_h2_total = K_p_h;
K_s_h2_total = K_s_h;

K_p_v3 = K_p_v;
K_p_h3 = K_p_h;
K_s_h3 = K_s_h;

u_sub_max_p = zeros(wavelength_num,1);
u_sub_max_s = zeros(wavelength_num,1);

u_air_max_p = zeros(wavelength_num,1);
u_air_max_s = zeros(wavelength_num,1);

for i = 1:wavelength_num

    if ne_bar(i,layer_num) > ne_bar(i,EML_position)
        u_sub_max_p(i) = ceil(u_data_num*ne_bar(i,layer_num)/ne_bar(i,EML_position))-1;
    else
        u_sub_max_p(i) = ceil(u_data_num*ne_bar(i,layer_num)/ne_bar(i,EML_position));
    end

    exp_phase = ones(1,u_sub_max_p(i));

    if u_sub_max_p(i) > u_data_num
        exp_phase(u_data_num+1:u_sub_max_p(i)) = exp((-4*pi*no_bar(i,EML_position)*sqrt(u(u_data_num+1:u_sub_max_p(i)).^2-1)*(thickness(EML_position-1)-z0))/wavelength(i));
    end

    K_p_v2(i,1:u_sub_max_p(i)) = 3/8*ne_bar(i,EML_position)*no_bar(i,layer_num)/no_bar(i,EML_position)^2*sqrt(1-(ne_bar(i,EML_position)*u(1:u_sub_max_p(i))/ne_bar(i,layer_num)).^2).*exp_phase.*u(1:u_sub_max_p(i)).^2.*abs((1+TMF_top.r_p(i,1:u_sub_max_p(i))).*TMF_bottom.t_p(i,1:u_sub_max_p(i))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2./abs(1-u(1:u_sub_max_p(i)).^2);
    K_p_h2(i,1:u_sub_max_p(i)) = 3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2*(1-(ne_bar(i,EML_position)*u(1:u_sub_max_p(i))/ne_bar(i,layer_num)).^2)).*exp_phase.*abs((1-TMF_top.r_p(i,1:u_sub_max_p(i))).*TMF_bottom.t_p(i,1:u_sub_max_p(i))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2/(12*(no_bar(i,EML_position)/ne_bar(i,EML_position))^2+4);

    if no_bar(i,layer_num) > no_bar(i,EML_position)
        u_sub_max_s(i) = ceil(u_data_num*no_bar(i,layer_num)/no_bar(i,EML_position))-1;
    else
        u_sub_max_s(i) = ceil(u_data_num*no_bar(i,layer_num)/no_bar(i,EML_position));
    end

    exp_phase = ones(1,u_sub_max_s(i));

    if u_sub_max_s(i) > u_data_num
        exp_phase(u_data_num+1:u_sub_max_s(i)) = exp((-4*pi*no_bar(i,EML_position)*sqrt(u(u_data_num+1:u_sub_max_s(i)).^2-1)*(thickness(EML_position-1)-z0))/wavelength(i));
    end

    K_s_h2(i,1:u_sub_max_s(i)) = 3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2-u(1:u_sub_max_s(i)).^2).*exp_phase.*abs((1+TMF_top.r_s(i,1:u_sub_max_s(i))).*TMF_bottom.t_s(i,1:u_sub_max_s(i))./(1-TMF_bottom.r_s(i,1:u_sub_max_s(i)).*TMF_top.r_s(i,1:u_sub_max_s(i)))).^2./((4*(ne_bar(i,EML_position)/no_bar(i,EML_position))^2+12)*abs(1-u(1:u_sub_max_s(i)).^2));

    K_p_v2_total(i,1:u_sub_max_p(i)) = 3/8*(u(1:u_sub_max_p(i)).^2).*real((1+TMF_bottom.r_p(i,1:u_sub_max_p(i))).*(1-conj(TMF_bottom.r_p(i,1:u_sub_max_p(i))))./sqrt(1-u(1:u_sub_max_p(i)).^2)).*abs((1+TMF_top.r_p(i,1:u_sub_max_p(i)))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2;
    K_p_h2_total(i,1:u_sub_max_p(i)) = 3*real((1-TMF_bottom.r_p(i,1:u_sub_max_p(i))).*(1+conj(TMF_bottom.r_p(i,1:u_sub_max_p(i)))).*sqrt(1-u(1:u_sub_max_p(i)).^2)).*abs((1-TMF_top.r_p(i,1:u_sub_max_p(i)))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2/(12*(no_bar(i,EML_position)/ne_bar(i,EML_position))^2+4);
    K_s_h2_total(i,1:u_sub_max_s(i)) = 3*real((1+TMF_bottom.r_s(i,1:u_sub_max_s(i))).*(1-conj(TMF_bottom.r_s(i,1:u_sub_max_s(i))))./sqrt(1-u(1:u_sub_max_s(i)).^2)).*abs((1+TMF_top.r_s(i,1:u_sub_max_s(i)))./(1-TMF_bottom.r_s(i,1:u_sub_max_s(i)).*TMF_top.r_s(i,1:u_sub_max_s(i)))).^2/(4*(ne_bar(i,EML_position)/no_bar(i,EML_position))^2+12);

    K_p_v3(i,:) = K_p_v2(i,:);
    K_p_h3(i,:) = K_p_h2(i,:);
    K_s_h3(i,:) = K_s_h2(i,:);

    if bottom_air_refractive_index(i) > ne_bar(i,EML_position)
        u_air_max_p(i) = min(u_sub_max_p(i),ceil(bottom_air_refractive_index(i)*u_data_num/ne_bar(i,EML_position))-1);
    else
        u_air_max_p(i) = min(u_sub_max_p(i),ceil(bottom_air_refractive_index(i)*u_data_num/ne_bar(i,EML_position)));
    end

    TMF_OLED_bottom_p = TMF_birefringence_whole_p(no_bar(i,layer_num:-1:1),ne_bar(i,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(i,EML_position)*u(1:u_air_max_p(i)),wavelength(i));

    R_p_bottom = abs(TMF_OLED_bottom_p.r_p).^2;

    cos_theta_sub = sqrt(1-(ne_bar(i,EML_position)*u(1:u_air_max_p(i))/ne_bar(i,layer_num)).^2);
    cos_theta_air = sqrt(1-(ne_bar(i,EML_position)*u(1:u_air_max_p(i))/bottom_air_refractive_index(i)).^2);

    r_p = (bottom_air_refractive_index(i)*cos_theta_sub-no_bar(i,layer_num)*cos_theta_air)./(bottom_air_refractive_index(i)*cos_theta_sub+no_bar(i,layer_num)*cos_theta_air);

    R_sub_air_bottom_p = abs(r_p).^2;
    T_sub_air_bottom_p = 1-R_sub_air_bottom_p;

    if bottom_air_refractive_index(i) > no_bar(i,EML_position)
        u_air_max_s(i) = min(u_sub_max_s(i),ceil(bottom_air_refractive_index(i)*u_data_num/no_bar(i,EML_position))-1);
    else
        u_air_max_s(i) = min(u_sub_max_s(i),ceil(bottom_air_refractive_index(i)*u_data_num/no_bar(i,EML_position)));
    end

    TMF_OLED_bottom_s = TMF_birefringence_whole_s(no_bar(i,layer_num:-1:1),ne_bar(i,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(i,EML_position)*u(1:u_air_max_s(i)),wavelength(i));

    R_s_bottom = abs(TMF_OLED_bottom_s.r_s).^2;

    cos_theta_sub = sqrt(1-(no_bar(i,EML_position)*u(1:u_air_max_s(i))/no_bar(i,layer_num)).^2);
    cos_theta_air = sqrt(1-(no_bar(i,EML_position)*u(1:u_air_max_s(i))/bottom_air_refractive_index(i)).^2);

    r_s = (no_bar(i,layer_num)*cos_theta_sub-bottom_air_refractive_index(i)*cos_theta_air)./(no_bar(i,layer_num)*cos_theta_sub+bottom_air_refractive_index(i)*cos_theta_air);

    R_sub_air_bottom_s = abs(r_s).^2;
    T_sub_air_bottom_s = 1-R_sub_air_bottom_s;

    K_p_v3(i,1:u_air_max_p(i)) = K_p_v2(i,1:u_air_max_p(i)).*T_sub_air_bottom_p./(1-R_p_bottom.*R_sub_air_bottom_p);
    K_p_h3(i,1:u_air_max_p(i)) = K_p_h2(i,1:u_air_max_p(i)).*T_sub_air_bottom_p./(1-R_p_bottom.*R_sub_air_bottom_p);
    K_s_h3(i,1:u_air_max_s(i)) = K_s_h2(i,1:u_air_max_s(i)).*T_sub_air_bottom_s./(1-R_s_bottom.*R_sub_air_bottom_s);

end

const  = (1-horizontal_dipole_ratio)*ne_bar(:,EML_position)./(wavelength.^4);
const2 = horizontal_dipole_ratio*no_bar(:,EML_position).*(3+(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2)./(4.*wavelength.^4);
const3 = const*u;
const4 = const2*u;

U_tot = 2*(const3.*K_p_v+const4.*(K_p_h+K_s_h));

U_bottom_transmit_p  = 2*(const3.*K_p_v2+const4.*K_p_h2);
U_bottom_transmit_ph = 2*(const4.*K_p_h2);
U_bottom_transmit_pv = 2*(const3.*K_p_v2);
U_bottom_transmit_s  = 2*const4.*K_s_h2;

U_bottom_transmit_total_p = 2*(const3.*K_p_v2_total+const4.*K_p_h2_total);
U_bottom_transmit_total_s = 2*const4.*K_s_h2_total;

U_bottom_transmit_thick_p  = 2*(const3.*K_p_v3+const4.*K_p_h3);
U_bottom_transmit_thick_ph = 2*(const4.*K_p_h3);
U_bottom_transmit_thick_pv = 2*(const3.*K_p_v3);
U_bottom_transmit_thick_s  = 2*const4.*K_s_h3;

const3 = repmat(const,1,u_num);
const4 = repmat(const2,1,u_num);

K_bottom_transmit_p = const3.*K_p_v2+const4.*K_p_h2;
K_bottom_transmit_s = const4.*K_s_h2;

K_bottom_transmit_thick_p = const3.*K_p_v3+const4.*K_p_h3;
K_bottom_transmit_thick_s = const4.*K_s_h3;

Power_ratio_air_matrix  = zeros(wavelength_num,1);
Power_ratio_air2_matrix = zeros(wavelength_num,1);
Power_ratio_sub_matrix  = zeros(wavelength_num,1);
Power_ratio_abs_matrix  = zeros(wavelength_num,1);
Power_ratio_wg_matrix   = zeros(wavelength_num,1);
Power_ratio_spp_matrix  = zeros(wavelength_num,1);

EQE_air_matrix  = zeros(wavelength_num,1);
EQE_air2_matrix = zeros(wavelength_num,1);
EQE_sub_matrix  = zeros(wavelength_num,1);

EQE_air_matrix_TE  = zeros(wavelength_num,1);
EQE_air_matrix_TMh = zeros(wavelength_num,1);
EQE_air_matrix_TMv = zeros(wavelength_num,1);
EQE_sub_matrix_TE  = zeros(wavelength_num,1);
EQE_sub_matrix_TMh = zeros(wavelength_num,1);
EQE_sub_matrix_TMv = zeros(wavelength_num,1);

EQE_abs_matrix = zeros(wavelength_num,1);
EQE_wg_matrix  = zeros(wavelength_num,1);
EQE_spp_matrix = zeros(wavelength_num,1);

sumUtot = sum(U_tot,2);

Purcell_factor = sumUtot./((const+const2)*u_data_num);
eta_eff = eta_rad*Purcell_factor./(1-eta_rad+eta_rad*Purcell_factor);

emissioneta_sumUtot       = emission_spectrum.*eta_eff./sumUtot;
lambdaemissioneta_sumUtot = wavelength.*emissioneta_sumUtot/sum(wavelength.*emission_spectrum);

P_air_p = zeros(wavelength_num,90);
P_air_s = zeros(wavelength_num,90);

P_sub_p = zeros(wavelength_num,90);
P_sub_s = zeros(wavelength_num,90);

const3 = pi*(const+const2);
const  = bottom_air_refractive_index./ne_bar(:,EML_position);
const2 = bottom_air_refractive_index./no_bar(:,EML_position);

for i = 1:wavelength_num

    Power_ratio_air_matrix(i)  = emissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_p(i,1:u_air_max_p(i)))+sum(U_bottom_transmit_thick_s(i,1:u_air_max_s(i))));
    Power_ratio_air2_matrix(i) = emissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,1:u_air_max_p(i)))+sum(U_bottom_transmit_s(i,1:u_air_max_s(i))));
    Power_ratio_sub_matrix(i)  = emissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_s(i,1:u_sub_max_s(i))));
    Power_ratio_abs_matrix(i)  = emissioneta_sumUtot(i)*(sum(U_bottom_transmit_total_p(i,1:u_sub_max_p(i))-U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_total_s(i,1:u_sub_max_s(i))-U_bottom_transmit_s(i,1:u_sub_max_s(i)))); % bottom 방향 abs만 계산
    Power_ratio_wg_matrix(i)   = emissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,u_sub_max_p(i)+1:u_data_num))+sum(U_bottom_transmit_s(i,u_sub_max_s(i)+1:u_data_num)));
    Power_ratio_spp_matrix(i)  = emissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,max(u_sub_max_p(i),u_data_num)+1:end))+sum(U_bottom_transmit_s(i,max(u_sub_max_s(i),u_data_num)+1:end)));

    EQE_air_matrix(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_p(i,1:u_air_max_p(i)))+sum(U_bottom_transmit_thick_s(i,1:u_air_max_s(i))));

    EQE_air_matrix_TE(i)  = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_s(i,1:u_air_max_s(i))));
    EQE_air_matrix_TMh(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_ph(i,1:u_air_max_p(i))));
    EQE_air_matrix_TMv(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_pv(i,1:u_air_max_p(i))));

    EQE_air2_matrix(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,1:u_air_max_p(i)))+sum(U_bottom_transmit_s(i,1:u_air_max_s(i))));
    EQE_sub_matrix(i)  = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_s(i,1:u_sub_max_s(i))));

    EQE_sub_matrix_TE(i)  = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_s(i,1:u_sub_max_s(i))));
    EQE_sub_matrix_TMh(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_ph(i,1:u_sub_max_p(i))));
    EQE_sub_matrix_TMv(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_pv(i,1:u_sub_max_p(i))));

    EQE_abs_matrix(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_total_p(i,1:u_sub_max_p(i))-U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_total_s(i,1:u_sub_max_s(i))-U_bottom_transmit_s(i,1:u_sub_max_s(i)))); % bottom 방향 abs만 계산
    EQE_wg_matrix(i)  = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,u_sub_max_p(i)+1:u_data_num))+sum(U_bottom_transmit_s(i,u_sub_max_s(i)+1:u_data_num)));
    EQE_spp_matrix(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,max(u_sub_max_p(i),u_data_num)+1:end))+sum(U_bottom_transmit_s(i,max(u_sub_max_s(i),u_data_num)+1:end)));

    P_air_p(i,:) = const(i)*spline(ne_bar(i,EML_position)*u(1:u_air_max_p(i)),sqrt(const(i)^2-u(1:u_air_max_p(i)).^2).*K_bottom_transmit_thick_p(i,1:u_air_max_p(i)),bottom_air_refractive_index(i)*sin089)/const3(i);
    P_air_s(i,:) = const2(i)*spline(no_bar(i,EML_position)*u(1:u_air_max_s(i)),sqrt(const2(i)^2-u(1:u_air_max_s(i)).^2).*K_bottom_transmit_thick_s(i,1:u_air_max_s(i)),bottom_air_refractive_index(i)*sin089)/const3(i);

    if bottom_air_refractive_index(i) > ne_bar(i,layer_num)
        P_air_p(i,ceil(asind(ne_bar(i,layer_num)/bottom_air_refractive_index(i)))+1:90) = 0;
    end

    if bottom_air_refractive_index(i) > no_bar(i,layer_num)
        P_air_s(i,ceil(asind(no_bar(i,layer_num)/bottom_air_refractive_index(i)))+1:90) = 0;
    end

    P_sub_p(i,:) = (ne_bar(i,layer_num)/ne_bar(i,EML_position))*spline(ne_bar(i,EML_position)*u(1:u_sub_max_p(i)),sqrt((ne_bar(i,layer_num)/ne_bar(i,EML_position))^2-u(1:u_sub_max_p(i)).^2).*K_bottom_transmit_p(i,1:u_sub_max_p(i)),ne_bar(i,layer_num)*sin089)/const3(i);
    P_sub_s(i,:) = (no_bar(i,layer_num)/no_bar(i,EML_position))*spline(no_bar(i,EML_position)*u(1:u_sub_max_s(i)),sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2-u(1:u_sub_max_s(i)).^2).*K_bottom_transmit_s(i,1:u_sub_max_s(i)),no_bar(i,layer_num)*sin089)/const3(i);

end

Power_ratio_air  = sum(Power_ratio_air_matrix);
Power_ratio_air2 = sum(Power_ratio_air2_matrix);
Power_ratio_sub  = sum(Power_ratio_sub_matrix);
Power_ratio_sub_confined = Power_ratio_sub-Power_ratio_air;
Power_ratio_abs  = sum(Power_ratio_abs_matrix);
Power_ratio_wg   = sum(Power_ratio_wg_matrix);
Power_ratio_spp  = sum(Power_ratio_spp_matrix);

EQE_air  = sum(EQE_air_matrix);
EQE_air2 = sum(EQE_air2_matrix);
EQE_sub  = sum(EQE_sub_matrix);
EQE_sub_confined = EQE_sub-EQE_air;
EQE_abs  = sum(EQE_abs_matrix);
EQE_wg   = sum(EQE_wg_matrix);
EQE_spp  = sum(EQE_spp_matrix);

P_air = P_air_p+P_air_s;
P_sub = P_sub_p+P_sub_s;

I_air = P_air.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);
I_sub = P_sub.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);

I_air_total = sum(I_air);
I_air_total = I_air_total/I_air_total(1);

I_sub_total  = sum(I_sub);
I_sub_sum    = sum(I_sub);
I_sub_sum_30 = sum(I_sub_sum(1,1:31).*sin089(1,1:31));
I_sub_total  = I_sub_total/I_sub_total(1);

EQE_factor_air = pi*sum(I_air_total.*sin089)/90;
EQE_factor_sub = pi*sum(I_sub_total.*sin089)/90;

spec_lambda = emission_spectrum./(wavelength*1e-9);

CE = 683*6.626*10^(-34)*(3*10^8)/(1.6*10^-(19))*sum(V_sel.*spec_lambda.*eta_eff.*P_air(:,1)./Purcell_factor);

I_FWHM = I_air(:,1);
I_FWHM = I_FWHM/max(I_FWHM);   % normalized 정면 spectrum
FWHM = sum(I_FWHM>=0.5);       % 정면 spectrum의 반치폭 (파장 구간을 좁히면 의미 제한적)

LEE_out_TE  = EQE_air_matrix_TE./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_out_TMh = EQE_air_matrix_TMh./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_out_TMv = EQE_air_matrix_TMv./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_sub_TE  = EQE_sub_matrix_TE./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_sub_TMh = EQE_sub_matrix_TMh./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_sub_TMv = EQE_sub_matrix_TMv./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_out = EQE_air_matrix./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_sub = EQE_sub_matrix./(lambdaemissioneta_sumUtot)./sumUtot;

%% Angular-range EQE
% Index convention: sin089 = sind(0:89), index i -> (i-1) deg
% Range [theta1, theta2): indices (theta1+1):(theta2)
aa = I_air_total.*sin089;
aa_total  = sum(aa(1:90));
EQE_0_20  = sum(aa(1:20))  / aa_total * EQE_air;   %  0~19 deg
EQE_20_40 = sum(aa(21:40)) / aa_total * EQE_air;   % 20~39 deg
EQE_40_60 = sum(aa(41:60)) / aa_total * EQE_air;   % 40~59 deg
EQE_60_80 = sum(aa(61:80)) / aa_total * EQE_air;   % 60~79 deg

% 기존 PSO 목적함수 값 (참고용, 최소화 기준이라 음수)
obj_4060 = -sum(aa(41:60))/aa_total*EQE_air;

%% ===== 결과 구조체 =====
res.thickness = thickness;
res.wavelength = wavelength;

res.EQE_air  = EQE_air;
res.EQE_air2 = EQE_air2;
res.EQE_sub  = EQE_sub;
res.EQE_sub_confined = EQE_sub_confined;
res.EQE_wg   = EQE_wg;
res.EQE_spp  = EQE_spp;
res.EQE_abs  = EQE_abs;

res.Power_ratio_air  = Power_ratio_air;
res.Power_ratio_air2 = Power_ratio_air2;
res.Power_ratio_sub  = Power_ratio_sub;
res.Power_ratio_sub_confined = Power_ratio_sub_confined;
res.Power_ratio_wg   = Power_ratio_wg;
res.Power_ratio_spp  = Power_ratio_spp;
res.Power_ratio_abs  = Power_ratio_abs;

res.EQE_0_20  = EQE_0_20;
res.EQE_20_40 = EQE_20_40;
res.EQE_40_60 = EQE_40_60;
res.EQE_60_80 = EQE_60_80;
res.obj_4060  = obj_4060;

res.CE   = CE;
res.FWHM = FWHM;
res.EQE_factor_air = EQE_factor_air;
res.EQE_factor_sub = EQE_factor_sub;
res.I_sub_sum_30   = I_sub_sum_30;

res.Purcell_factor = Purcell_factor;
res.eta_eff        = eta_eff;
res.emission_spectrum = emission_spectrum;

res.I_air_total = I_air_total;      % 각도별 정규화 세기 (0~89 deg)
res.I_sub_total = I_sub_total;
res.I_air_front = I_air(:,1);       % 정면 스펙트럼
res.I_air = I_air;
res.I_sub = I_sub;

res.LEE_out = LEE_out;
res.LEE_sub = LEE_sub;
res.LEE_out_TE  = LEE_out_TE;
res.LEE_out_TMh = LEE_out_TMh;
res.LEE_out_TMv = LEE_out_TMv;
res.LEE_sub_TE  = LEE_sub_TE;
res.LEE_sub_TMh = LEE_sub_TMh;
res.LEE_sub_TMv = LEE_sub_TMv;

end

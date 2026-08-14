%% PSO 기반 OLED 두께 최적화 (파장 595~600 nm 축소 버전)
% 원본 대비 변경점:
%   1) 계산 파장을 400~800 nm (401 pts) -> 595~600 nm (6 pts) 로 축소
%   2) n/k 데이터, 발광 스펙트럼, 시감도(V) 를 모두 동일한 인덱스로 슬라이싱
%   3) mat 파일을 worker 당 1회만 로드하도록 persistent 캐싱
%   4) 결과 출력부를 실제 변수 개수(nvars=3)에 맞게 수정

clear; clc;

%% --- 파장 설정 (여기만 바꾸면 전체 계산 범위가 바뀝니다) ---
% 원본 nk / 스펙트럼 데이터는 400~800 nm 를 1 nm 간격 401 포인트로 저장하고 있으므로
% 파장 lambda 에 대응하는 행 인덱스는 (lambda - 399) 입니다.
LAMBDA_MIN = 595;
LAMBDA_MAX = 600;

nvars = 3;       % 최적화할 변수 개수
lb = [ 10, 10,  5];
ub = [250,250, 30];

% --- 2. 병렬 처리 설정 ---
numWorkers = 4;
if isempty(gcp('nocreate'))
    parpool('local', numWorkers);
end

options = optimoptions('particleswarm', ...
    'SwarmSize', 30, ...       % 파티클(입자) 개수 (기본값은 10*nvars)
    'MaxIterations', 5, ...    % 최대 반복 횟수
    'PlotFcn', @pswplotbestf, ...
    'UseParallel', true, ...
    'HybridFcn', @fmincon, ...
    'OutputFcn', @realtime_pso_monitor); % 실시간으로 최적값 그래프 표시

% --- 3. PSO 실행 ---
fprintf('PSO 최적화를 시작합니다... (계산 파장: %d ~ %d nm)\n', LAMBDA_MIN, LAMBDA_MAX);

objfun = @(x) evaluateEQE(x, LAMBDA_MIN, LAMBDA_MAX);

[x_optimal, fval] = particleswarm(objfun, nvars, lb, ub, options);

% --- 4. 결과 출력 ---
fprintf('------------------------------------\n');
fprintf('최적화가 완료되었습니다.\n');
fprintf('  > 최적의 dETL: %.4f\n', x_optimal(1));
fprintf('  > 최적의 dHTL: %.4f\n', x_optimal(2));
fprintf('  > 최적의 dCap: %.4f\n', x_optimal(3));
fprintf('  > 최대 Totalpower (fval의 음수): %.4f\n', -fval);
fprintf('------------------------------------\n');


%% --- 목적 함수 (Objective Function) ---

function output = evaluateEQE(point, lambda_min, lambda_max)

if nargin < 2 || isempty(lambda_min), lambda_min = 595; end
if nargin < 3 || isempty(lambda_max), lambda_max = 600; end

% mat 파일은 worker 당 한 번만 로드 (반복 load 제거로 속도 향상)
persistent material spectrum V_401
if isempty(material)
    S_nk       = load('nk_JH_total.mat');       % material, spectrum 포함
    material   = S_nk.material;
    spectrum   = S_nk.spectrum;
    S_photopic = load('Photopic_400_800.mat');  % V_401 포함
    V_401      = S_photopic.V_401;
end

dETL = point(1);
dHTL = point(2);
dCap = point(3);

% Copyright ⓒ All Rights Reserved.

%% ===== 파장 축소: 데이터 행 인덱스 계산 =====
% 원본 데이터는 400 nm 가 1행 -> lambda 의 행 인덱스 = lambda - 399
idx = (lambda_min - 399):(lambda_max - 399);

wavelength     = (lambda_min:lambda_max).';
wavelength_num = length(wavelength);

V_sel = V_401(idx, 1);   % 축소된 시감도 (CE 계산용)

%% ===== 발광 스펙트럼 =====
emission_spectrum = spectrum.l_I_Irdmppyph2tmd; %% JOSong, eta_rad=0.98, hdr=0.865
emission_spectrum = emission_spectrum(idx, 1);
emission_spectrum = emission_spectrum / sum(emission_spectrum);

eta_rad = 1;
horizontal_dipole_ratio = 0.95;
bottom_air_refractive_index = ones(wavelength_num, 1);

%% JO 구조 %%
% 모든 물질 n/k 도 동일한 idx 로 슬라이싱
no_bar = [ones(wavelength_num,1) material.l_Ag_McPeak(idx,1) material.l_B3_o_JO(idx,1) material.l_TCTA_B3_o_JO(idx,1) material.l_TCTA_o_JO(idx,1) material.l_TAPC_o_JO(idx,1) material.l_Ag_McPeak(idx,1) 2.3*ones(wavelength_num,1) ones(wavelength_num,1)];
ne_bar = [ones(wavelength_num,1) material.l_Ag_McPeak(idx,1) material.l_B3_e_JO(idx,1) material.l_TCTA_B3_e_JO(idx,1) material.l_TCTA_e_JO(idx,1) material.l_TAPC_e_JO(idx,1) material.l_Ag_McPeak(idx,1) 2.3*ones(wavelength_num,1) ones(wavelength_num,1)];

d1 = 70;   %Anode TAPC
Nd1 = length(d1);
d2 = 197;  %HTL TPBi
Nd2 = length(d2);
data_matrix = zeros(Nd1*Nd2, 12);

for k1 = 1:length(d1)
    for k2 = 1:length(d2)

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
        U_bottom_transmit    = U_bottom_transmit_p+U_bottom_transmit_s;

        U_bottom_transmit_total_p = 2*(const3.*K_p_v2_total+const4.*K_p_h2_total);
        U_bottom_transmit_total_s = 2*const4.*K_s_h2_total;
        U_bottom_transmit_total   = U_bottom_transmit_total_p+U_bottom_transmit_total_s;

        U_bottom_transmit_thick_p  = 2*(const3.*K_p_v3+const4.*K_p_h3);
        U_bottom_transmit_thick_ph = 2*(const4.*K_p_h3);
        U_bottom_transmit_thick_pv = 2*(const3.*K_p_v3);
        U_bottom_transmit_thick_s  = 2*const4.*K_s_h3;
        U_bottom_transmit_thick    = U_bottom_transmit_thick_p+U_bottom_transmit_thick_s;

        const3 = repmat(const,1,u_num);
        const4 = repmat(const2,1,u_num);

        K_bottom_transmit_p = const3.*K_p_v2+const4.*K_p_h2;
        K_bottom_transmit_s = const4.*K_s_h2;
        K_bottom_transmit   = K_bottom_transmit_p+K_bottom_transmit_s;

        K_bottom_transmit_thick_p = const3.*K_p_v3+const4.*K_p_h3;
        K_bottom_transmit_thick_s = const4.*K_s_h3;
        K_bottom_transmit_thick   = K_bottom_transmit_thick_p+K_bottom_transmit_thick_s;

        Power_ratio_air_matrix  = zeros(wavelength_num,1);
        Power_ratio_air2_matrix = zeros(wavelength_num,1);
        Power_ratio_sub_matrix  = zeros(wavelength_num,1);
        Power_ratio_abs_matrix  = zeros(wavelength_num,1);
        Power_ratio_wg_matrix   = zeros(wavelength_num,1);
        Power_ratio_spp_matrix  = zeros(wavelength_num,1);

        EQE_air_matrix  = zeros(wavelength_num,1);
        EQE_air2_matrix = zeros(wavelength_num,1);
        EQE_sub_matrix  = zeros(wavelength_num,1);
        %%
        EQE_air_matrix_TE  = zeros(wavelength_num,1);
        EQE_air_matrix_TMh = zeros(wavelength_num,1);
        EQE_air_matrix_TMv = zeros(wavelength_num,1);
        EQE_sub_matrix_TE  = zeros(wavelength_num,1);
        EQE_sub_matrix_TMh = zeros(wavelength_num,1);
        EQE_sub_matrix_TMv = zeros(wavelength_num,1);
        %%
        EQE_abs_matrix = zeros(wavelength_num,1);
        EQE_wg_matrix  = zeros(wavelength_num,1);
        EQE_spp_matrix = zeros(wavelength_num,1);

        sumUtot = sum(U_tot,2);

        Purcell_factor = sumUtot./((const+const2)*u_data_num);
        eta_eff = eta_rad*Purcell_factor./(1-eta_rad+eta_rad*Purcell_factor);

        emission_spectrum = emission_spectrum/sum(emission_spectrum);

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
            %%
            EQE_air_matrix_TE(i)  = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_s(i,1:u_air_max_s(i))));
            EQE_air_matrix_TMh(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_ph(i,1:u_air_max_p(i))));
            EQE_air_matrix_TMv(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_pv(i,1:u_air_max_p(i))));
            %%
            EQE_air2_matrix(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,1:u_air_max_p(i)))+sum(U_bottom_transmit_s(i,1:u_air_max_s(i))));
            EQE_sub_matrix(i)  = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_s(i,1:u_sub_max_s(i))));
            %%
            EQE_sub_matrix_TE(i)  = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_s(i,1:u_sub_max_s(i))));
            EQE_sub_matrix_TMh(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_ph(i,1:u_sub_max_p(i))));
            EQE_sub_matrix_TMv(i) = lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_pv(i,1:u_sub_max_p(i))));
            %%
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

        I_air   = P_air.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);
        I_air_p = P_air_p.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);

        I_air_total = sum(I_air);
        I_air_sum   = sum(I_air);
        I_air_total = I_air_total/I_air_total(1);

        I_sub   = P_sub.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);
        I_sub_p = P_sub_p.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);
        I_sub_total  = sum(I_sub);
        I_sub_sum    = sum(I_sub);
        I_sub_sum_30 = sum(I_sub_sum(1,1:31).*sin089(1,1:31));
        I_sub_total  = I_sub_total/I_sub_total(1);

        EQE_factor_air = pi*sum(I_air_total.*sin089)/90;
        EQE_factor_sub = pi*sum(I_sub_total.*sin089)/90;

        spec_lambda = zeros(wavelength_num,1);
        for i = 1:wavelength_num
            spec_lambda(i,1) = emission_spectrum(i,1)/(wavelength(i)*10^(-9));
        end

        % 축소된 파장 구간에 대한 CE (참고용)
        CE = 683*6.626*10^(-34)*(3*10^8)/(1.6*10^-(19))*sum(V_sel.*spec_lambda.*eta_eff.*P_air(:,1)./Purcell_factor);

        I_FWHM = I_air(:,1);
        I_FWHM = I_FWHM/max(I_FWHM);  % normalized 정면 spectrum
        FWHM = sum(I_FWHM>=0.5);      % 정면 spectrum의 반치폭 (파장 구간 축소 시 의미 제한적)

        k_all = length(d2)*(k1-1)+k2;

        data_matrix(k_all,:) = [d1(k1),d2(k2),FWHM,EQE_air,EQE_sub_confined,EQE_wg,EQE_spp,EQE_abs,EQE_sub,EQE_factor_air,CE,I_sub_sum_30];

    end
end

LEE_out_TE  = EQE_air_matrix_TE./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_out_TMh = EQE_air_matrix_TMh./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_out_TMv = EQE_air_matrix_TMv./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_sub_TE  = EQE_sub_matrix_TE./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_sub_TMh = EQE_sub_matrix_TMh./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_sub_TMv = EQE_sub_matrix_TMv./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_out = EQE_air_matrix./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_sub = EQE_sub_matrix./(lambdaemissioneta_sumUtot)./sumUtot;

U_test = reshape(U_tot,[wavelength_num*u_num,1]);
w  = repmat(wavelength, [u_num 1]);
uu = zeros(wavelength_num*u_num,1);

for i = 1:u_num
    uu(wavelength_num*(i-1)+1:wavelength_num*i,1) = u(i);
end

aa = I_air_total.*sin089;
output = -sum(aa(41:60))/sum(aa(1:90))*EQE_air;

%% Angular-range EQE
% Index convention: sin089 = sind(0:89), index i -> (i-1) deg
% Range [theta1, theta2): indices (theta1+1):(theta2)
aa_total  = sum(aa(1:90));
EQE_0_20  = sum(aa(1:20))  / aa_total * EQE_air;   %  0~19 deg
EQE_20_40 = sum(aa(21:40)) / aa_total * EQE_air;   % 20~39 deg
EQE_40_60 = sum(aa(41:60)) / aa_total * EQE_air;   % 40~59 deg
EQE_60_80 = sum(aa(61:80)) / aa_total * EQE_air;   % 60~79 deg

end


function stop = realtime_pso_monitor(optimValues, state)

% Iteration별 최고 기록을 저장할 persistent 변수
persistent iteration_log nvars_local

stop = false; % 최적화 중지 여부 (기본값)

% --- 1. 'init' 상태 (최적화 시작 시 1회 호출) ---
if strcmp(state, 'init')

    iteration_log = [];

    nvars_local = length(optimValues.bestx);

    if nvars_local > 0
        iteration_log = zeros(0, nvars_local + 1);
    end

    fprintf('--- 실시간 PSO 모니터링 및 중간 저장 시작 ---\n');
    fprintf(' Iter  |  Current Max Power (Totalpower) | 저장 상태\n');
    fprintf('----------------------------------------------------------\n');

    % --- 2. 'iter' 상태 (매 반복 종료 시 호출) ---
elseif strcmp(state, 'iter')

    iter = optimValues.iteration;

    current_best_point = optimValues.bestx;
    current_best_fval  = optimValues.bestfval;

    fprintf(' %5d |      %.6f                 |', iter, -current_best_fval);

    if isempty(nvars_local) || nvars_local == 0
        nvars_local = length(current_best_point);
        iteration_log = zeros(0, nvars_local + 1);
    end

    new_row = [current_best_point, current_best_fval];
    iteration_log(end+1, :) = new_row;

    try
        save('260330_data_in_progress_JinHKim_bare_4060_opt.mat', 'iteration_log');
        fprintf(' 저장됨\n');
    catch ME
        fprintf(' 저장 실패: %s\n', ME.message);
    end

    % --- 3. 'done' 상태 (최적화 완료 시 1회 호출) ---
elseif strcmp(state, 'done')
    fprintf('----------------------------------------------------------\n');
    fprintf('실시간 모니터링 및 중간 저장이 완료되었습니다.\n');
    fprintf('최종 iteration별 최고 기록이 "260330_data_in_progress_JinHKim_bare_4060_opt.mat"에 저장되었습니다.\n');
end
end

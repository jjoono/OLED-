% Made by WC Lee
% clear;
clc;
tic;

load('nk_JH_total.mat')
% load('nk_JHH.mat')
load('Photopic_400_800.mat')
% Copyright �� All Rights Reserved.
load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat')
BSDF=BSDF_MLA(:,:,11);
% <<< CHECK ME >>> the file name reads n_MLA = 1.30 : 0.05 : 2.00, so slice 11 is
% n_MLA = 1.80 while the substrate below is 1.5 (no_bar(end)). Confirm the slice that
% matches your substrate index, or the BSDF describes the wrong interface.

% --- one-off diagnostic on the BSDF (safe to delete): every column
% --- (incidence angle) should sum to 1 if the BSDF conserves energy.
col_sum = sum(BSDF,1);
fprintf('BSDF energy: column sums %.4f .. %.4f\n', min(col_sum), max(col_sum));
% -----------------------------------------------------------------------------

wavelength=550;
wavelength_num=length(wavelength);

emission_spectrum=1;


eta_rad=1;
horizontal_dipole_ratio=2/3;
bottom_air_refractive_index=ones(wavelength_num,1);


%% JO ���� %%
no_bar=[ones(wavelength_num,1) material.l_Ag_McPeak(151) material.l_B3_o_JO(151) material.l_TCTA_B3_o_JO(151) material.l_TAPC_o_JO(151) 1.86+0.003i 1.5];
ne_bar=[ones(wavelength_num,1) material.l_Ag_McPeak(151) material.l_B3_e_JO(151) material.l_TCTA_B3_e_JO(151) material.l_TAPC_e_JO(151) 1.86+0.003i 1.5];

d1= 10:10:500;  %0:5:200; %Anode TAPC
Nd1=length(d1);
d2= 10:10:500;  %0:5:200; %HTL TPBi
Nd2=length(d2);
% data_matrix columns:
%   1 d1  2 d2  3 EQE_air  4 EQE_sub_confined  5 EQE_wg  6 EQE_spp  7 EQE_abs
%   8 EQE_MLA (substrate light extracted through the MLA, R_step round trips)
%   9 EQE_sub (= EQE_air + EQE_sub_confined), so that eta_ext = EQE_MLA/EQE_sub
data_matrix = zeros(Nd1*Nd2,9);
disp_matrix = zeros(Nd1*Nd2,7);
inplane_matrix = zeros(Nd1*Nd2,3001);
for k1=1:length(d1)
    for k2=1:length(d2)
        thickness= [100 d1(k1) 25 d2(k2) 50];

        EML_position=3;
        % <<< CHECK ME >>> layer 4 (TCTA:B3, 25 nm) looks like the emitting layer and
        % z0 = 12.5 = 25/2 is its centre. With EML_position = 3 the emitter instead sits
        % 12.5 nm inside the B3PyMPM layer that touches the Ag. Set EML_position = 4 if
        % the 25 nm layer is the EML. Left at 3 so this file reproduces your current run.

        %         z0=d_slice/2; % emitting position in EML
        z0=12.5;
        u_data_num=1000;
        max_u=3;

        %%

        sin089=sind(0:89);
        cos089=cosd(0:89);
        layer_num=size(no_bar,2);
        u=[(0:u_data_num-1)/u_data_num (u_data_num+1:u_data_num*max_u)/u_data_num];
        u_num=length(u);

        TMF_bottom=TMF_birefringence_whole(no_bar(:,EML_position:layer_num),ne_bar(:,EML_position:layer_num),[thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0],u,wavelength);
        TMF_top=TMF_birefringence_whole(no_bar(:,EML_position:-1:1),ne_bar(:,EML_position:-1:1),[z0 thickness(EML_position-2:-1:1) 0],u,wavelength);

        K_p_v=3/4*real(ne_bar(:,EML_position)./no_bar(:,EML_position)*(u.^2./sqrt(1-u.^2)).*(1+TMF_bottom.r_p).*(1+TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
        K_p_h=real(3./(6*(no_bar(:,EML_position)./ne_bar(:,EML_position)).^2+2)*sqrt(1-u.^2).*(1-TMF_bottom.r_p).*(1-TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
        K_s_h=real(3./(2*(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2+6)*(1./sqrt(1-u.^2)).*(1+TMF_bottom.r_s).*(1+TMF_top.r_s)./(1-TMF_bottom.r_s.*TMF_top.r_s));

        K_p_v2=K_p_v; K_p_h2=K_p_h; K_s_h2=K_s_h;
        K_p_v3=K_p_v; K_p_h3=K_p_h; K_s_h3=K_s_h;

        u_sub_max_p=zeros(wavelength_num,1); u_sub_max_s=zeros(wavelength_num,1);
        u_air_max_p=zeros(wavelength_num,1); u_air_max_s=zeros(wavelength_num,1);

        for i=1:wavelength_num
            if ne_bar(i,layer_num)>ne_bar(i,EML_position), u_sub_max_p(i)=ceil(u_data_num*ne_bar(i,layer_num)/ne_bar(i,EML_position))-1; else, u_sub_max_p(i)=ceil(u_data_num*ne_bar(i,layer_num)/ne_bar(i,EML_position)); end
            exp_phase=ones(1,u_sub_max_p(i));
            if u_sub_max_p(i)>u_data_num, exp_phase(u_data_num+1:u_sub_max_p(i))=exp((-4*pi*no_bar(i,EML_position)*sqrt(u(u_data_num+1:u_sub_max_p(i)).^2-1)*(thickness(EML_position-1)-z0))/wavelength(i)); end
            K_p_v2(i,1:u_sub_max_p(i))=3/8*ne_bar(i,EML_position)*no_bar(i,layer_num)/no_bar(i,EML_position)^2*sqrt(1-(ne_bar(i,EML_position)*u(1:u_sub_max_p(i))/ne_bar(i,layer_num)).^2).*exp_phase.*u(1:u_sub_max_p(i)).^2.*abs((1+TMF_top.r_p(i,1:u_sub_max_p(i))).*TMF_bottom.t_p(i,1:u_sub_max_p(i))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2./abs(1-u(1:u_sub_max_p(i)).^2);
            K_p_h2(i,1:u_sub_max_p(i))=3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2*(1-(ne_bar(i,EML_position)*u(1:u_sub_max_p(i))/ne_bar(i,layer_num)).^2)).*exp_phase.*abs((1-TMF_top.r_p(i,1:u_sub_max_p(i))).*TMF_bottom.t_p(i,1:u_sub_max_p(i))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2/(12*(no_bar(i,EML_position)/ne_bar(i,EML_position))^2+4);

            if no_bar(i,layer_num)>no_bar(i,EML_position), u_sub_max_s(i)=ceil(u_data_num*no_bar(i,layer_num)/no_bar(i,EML_position))-1; else, u_sub_max_s(i)=ceil(u_data_num*no_bar(i,layer_num)/no_bar(i,EML_position)); end
            exp_phase=ones(1,u_sub_max_s(i));
            if u_sub_max_s(i)>u_data_num, exp_phase(u_data_num+1:u_sub_max_s(i))=exp((-4*pi*no_bar(i,EML_position)*sqrt(u(u_data_num+1:u_sub_max_s(i)).^2-1)*(thickness(EML_position-1)-z0))/wavelength(i)); end
            K_s_h2(i,1:u_sub_max_s(i))=3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2-u(1:u_sub_max_s(i)).^2).*exp_phase.*abs((1+TMF_top.r_s(i,1:u_sub_max_s(i))).*TMF_bottom.t_s(i,1:u_sub_max_s(i))./(1-TMF_bottom.r_s(i,1:u_sub_max_s(i)).*TMF_top.r_s(i,1:u_sub_max_s(i)))).^2./((4*(ne_bar(i,EML_position)/no_bar(i,EML_position))^2+12)*abs(1-u(1:u_sub_max_s(i)).^2));

            K_p_v3(i,:)=K_p_v2(i,:); K_p_h3(i,:)=K_p_h2(i,:); K_s_h3(i,:)=K_s_h2(i,:);

            if bottom_air_refractive_index(i)>ne_bar(i,EML_position), u_air_max_p(i)=min(u_sub_max_p(i),ceil(bottom_air_refractive_index(i)*u_data_num/ne_bar(i,EML_position))-1); else, u_air_max_p(i)=min(u_sub_max_p(i),ceil(bottom_air_refractive_index(i)*u_data_num/ne_bar(i,EML_position))); end
            TMF_OLED_bottom_p=TMF_birefringence_whole_p(no_bar(i,layer_num:-1:1),ne_bar(i,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(i,EML_position)*u(1:u_air_max_p(i)),wavelength(i));
            R_p_bottom=abs(TMF_OLED_bottom_p.r_p).^2;
            cos_theta_sub=sqrt(1-(ne_bar(i,EML_position)*u(1:u_air_max_p(i))/ne_bar(i,layer_num)).^2);
            cos_theta_air=sqrt(1-(ne_bar(i,EML_position)*u(1:u_air_max_p(i))/bottom_air_refractive_index(i)).^2);
            r_p=(bottom_air_refractive_index(i)*cos_theta_sub-no_bar(i,layer_num)*cos_theta_air)./(bottom_air_refractive_index(i)*cos_theta_sub+no_bar(i,layer_num)*cos_theta_air);
            R_sub_air_bottom_p=abs(r_p).^2; T_sub_air_bottom_p=1-R_sub_air_bottom_p;

            if bottom_air_refractive_index(i)>no_bar(i,EML_position), u_air_max_s(i)=min(u_sub_max_s(i),ceil(bottom_air_refractive_index(i)*u_data_num/no_bar(i,EML_position))-1); else, u_air_max_s(i)=min(u_sub_max_s(i),ceil(bottom_air_refractive_index(i)*u_data_num/no_bar(i,EML_position))); end
            TMF_OLED_bottom_s=TMF_birefringence_whole_s(no_bar(i,layer_num:-1:1),ne_bar(i,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(i,EML_position)*u(1:u_air_max_s(i)),wavelength(i));
            R_s_bottom=abs(TMF_OLED_bottom_s.r_s).^2;
            cos_theta_sub=sqrt(1-(no_bar(i,EML_position)*u(1:u_air_max_s(i))/no_bar(i,layer_num)).^2);
            cos_theta_air=sqrt(1-(no_bar(i,EML_position)*u(1:u_air_max_s(i))/bottom_air_refractive_index(i)).^2);
            r_s=(no_bar(i,layer_num)*cos_theta_sub-bottom_air_refractive_index(i)*cos_theta_air)./(no_bar(i,layer_num)*cos_theta_sub+bottom_air_refractive_index(i)*cos_theta_air);
            R_sub_air_bottom_s=abs(r_s).^2; T_sub_air_bottom_s=1-R_sub_air_bottom_s;

            K_p_v3(i,1:u_air_max_p(i))=K_p_v2(i,1:u_air_max_p(i)).*T_sub_air_bottom_p./(1-R_p_bottom.*R_sub_air_bottom_p);
            K_p_h3(i,1:u_air_max_p(i))=K_p_h2(i,1:u_air_max_p(i)).*T_sub_air_bottom_p./(1-R_p_bottom.*R_sub_air_bottom_p);
            K_s_h3(i,1:u_air_max_s(i))=K_s_h2(i,1:u_air_max_s(i)).*T_sub_air_bottom_s./(1-R_s_bottom.*R_sub_air_bottom_s);
        end

        const=(1-horizontal_dipole_ratio)*ne_bar(:,EML_position)./(wavelength.^4);
        const2=horizontal_dipole_ratio*no_bar(:,EML_position).*(3+(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2)./(4.*wavelength.^4);
        const3=const*u;
        const4=const2*u;

        U_tot=2*(const3.*K_p_v+const4.*(K_p_h+K_s_h));
        % polarisation-resolved total dissipation (U_tot = U_tot_p + U_tot_s).
        % p and s use DIFFERENT normalisations of u -- u_p = k_x/(k0 n_e) and
        % u_s = k_x/(k0 n_o) -- so with a birefringent EML they have different
        % substrate cut-offs and the channels must be binned separately.
        U_tot_p=2*(const3.*K_p_v+const4.*K_p_h);
        U_tot_s=2*const4.*K_s_h;
        U_bottom_transmit_p=2*(const3.*K_p_v2+const4.*K_p_h2);
        U_bottom_transmit_s=2*const4.*K_s_h2;
        U_bottom_transmit_thick_p=2*(const3.*K_p_v3+const4.*K_p_h3);
        U_bottom_transmit_thick_s=2*const4.*K_s_h3;

        sumUtot=sum(U_tot,2);

        u_free_range = 1:u_data_num;
        K_p_v_free=3/4*real(ne_bar(:,EML_position)./no_bar(:,EML_position)*(u(u_free_range).^2./sqrt(1-u(u_free_range).^2)));
        K_p_h_free=real(3./(6*(no_bar(:,EML_position)./ne_bar(:,EML_position)).^2+2)*sqrt(1-u(u_free_range).^2));
        K_s_h_free=real(3./(2*(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2+6)*(1./sqrt(1-u(u_free_range).^2)));
        const3_free = const * u(u_free_range);
        const4_free = const2 * u(u_free_range);
        U_tot_free = 2*(const3_free.*K_p_v_free + const4_free.*(K_p_h_free + K_s_h_free));
        P_free = sum(U_tot_free, 2);

        Purcell_factor = sumUtot ./ P_free;
        eta_eff = eta_rad*Purcell_factor./(1-eta_rad+eta_rad*Purcell_factor);

        EQE_air_matrix=zeros(wavelength_num,1);
        EQE_sub_confined_matrix=zeros(wavelength_num,1);
        EQE_wg_matrix=zeros(wavelength_num,1);
        EQE_spp_matrix=zeros(wavelength_num,1);
        EQE_abs_matrix=zeros(wavelength_num,1);

        for i=1:wavelength_num
            % the three ranges partition 1:u_num separately for each polarisation
            ip_air=1:u_air_max_p(i);                       is_air=1:u_air_max_s(i);
            ip_sub=1:u_sub_max_p(i);                       is_sub=1:u_sub_max_s(i);
            ip_wg =u_sub_max_p(i)+1:u_data_num;            is_wg =u_sub_max_s(i)+1:u_data_num;
            ip_spp=max(u_sub_max_p(i),u_data_num)+1:u_num; is_spp=max(u_sub_max_s(i),u_data_num)+1:u_num;

            P_air    = sum(U_bottom_transmit_thick_p(i,ip_air))+sum(U_bottom_transmit_thick_s(i,is_air));
            P_sub    = sum(U_bottom_transmit_p(i,ip_sub))      +sum(U_bottom_transmit_s(i,is_sub));
            P_totsub = sum(U_tot_p(i,ip_sub))                  +sum(U_tot_s(i,is_sub));

            EQE_air_matrix(i)          = P_air;
            EQE_sub_confined_matrix(i) = P_sub - P_air;
            EQE_abs_matrix(i)          = P_totsub - P_sub;     % absorbed, both sides
            EQE_wg_matrix(i)           = sum(U_tot_p(i,ip_wg)) +sum(U_tot_s(i,is_wg));
            EQE_spp_matrix(i)          = sum(U_tot_p(i,ip_spp))+sum(U_tot_s(i,is_spp));
        end

        photon_weight = emission_spectrum.*eta_eff./sumUtot;
        % per-wavelength power delivered into the substrate mode (= air + sub-confined);
        % this is P0, the source term of the MLA recycling series below
        EQE_sub_matrix = photon_weight .* (EQE_air_matrix + EQE_sub_confined_matrix);
        EQE_air = sum(photon_weight .* EQE_air_matrix);
        EQE_sub_confined = sum(photon_weight .* EQE_sub_confined_matrix);

        EQE_wg  = sum(photon_weight .* EQE_wg_matrix);
        EQE_spp = sum(photon_weight .* EQE_spp_matrix);
        EQE_abs = sum(photon_weight .* EQE_abs_matrix);

        % With the ranges above the five channels add up to eta_eff by construction,
        % so nothing is rescaled; the residual is printed as a check instead.
        EQE_closure  = EQE_air + EQE_sub_confined + EQE_wg + EQE_spp + EQE_abs;
        EQE_residual = sum(emission_spectrum .* eta_eff) - EQE_closure;

        %% ======================= MLA out-coupling =======================
        % Photon recycling on the extraction surface, following
        % planar_Sweep22_MLA_JH_full_lambda_modified.m:
        %   the substrate-delivered power P0 leaves with angular distribution
        %   Psub_norm, the BSDF transmits part of it and sends the rest back,
        %   the OLED stack reflects the returned light with R_LED, and the
        %   series is summed over R_step terms.

        % --- kernels without the factor 2 of U_bottom_* (reference convention)
        K_bottom_transmit_p       = const3.*K_p_v2 + const4.*K_p_h2;
        K_bottom_transmit_s       = const4.*K_s_h2;
        K_bottom_transmit_thick_p = const3.*K_p_v3 + const4.*K_p_h3;
        K_bottom_transmit_thick_s = const4.*K_s_h3;
        const_free = pi*(const + const2);        % pi * (free-space dissipation per unit u)

        % --- angular distribution in the substrate and in air, on a 0-89 deg grid
        P_sub_p = zeros(wavelength_num,90); P_sub_s = zeros(wavelength_num,90);
        P_air_p = zeros(wavelength_num,90); P_air_s = zeros(wavelength_num,90);
        for i=1:wavelength_num
            rp = ne_bar(i,layer_num)/ne_bar(i,EML_position);
            up = u(1:u_sub_max_p(i));
            P_sub_p(i,:) = rp*interp1(ne_bar(i,EML_position)*up, ...
                sqrt(max(0,rp^2-up.^2)).*K_bottom_transmit_p(i,1:u_sub_max_p(i)), ...
                ne_bar(i,layer_num)*sin089,'spline',0)/const_free(i);

            rs = no_bar(i,layer_num)/no_bar(i,EML_position);
            us = u(1:u_sub_max_s(i));
            P_sub_s(i,:) = rs*interp1(no_bar(i,EML_position)*us, ...
                sqrt(max(0,rs^2-us.^2)).*K_bottom_transmit_s(i,1:u_sub_max_s(i)), ...
                no_bar(i,layer_num)*sin089,'spline',0)/const_free(i);

            ap = bottom_air_refractive_index(i)/ne_bar(i,EML_position);
            uap = u(1:u_air_max_p(i));
            P_air_p(i,:) = ap*interp1(ne_bar(i,EML_position)*uap, ...
                sqrt(max(0,ap^2-uap.^2)).*K_bottom_transmit_thick_p(i,1:u_air_max_p(i)), ...
                bottom_air_refractive_index(i)*sin089,'spline',0)/const_free(i);

            as_ = bottom_air_refractive_index(i)/no_bar(i,EML_position);
            uas = u(1:u_air_max_s(i));
            P_air_s(i,:) = as_*interp1(no_bar(i,EML_position)*uas, ...
                sqrt(max(0,as_^2-uas.^2)).*K_bottom_transmit_thick_s(i,1:u_air_max_s(i)), ...
                bottom_air_refractive_index(i)*sin089,'spline',0)/const_free(i);
        end
        P_sub = P_sub_p + P_sub_s;
        P_air = P_air_p + P_air_s;
        I_sub = P_sub.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);
        I_air = P_air.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);
        I_air_total = sum(I_air,1); I_air_total = I_air_total/I_air_total(1);

        Psub_norm = I_sub.*sin089;                       % solid-angle weight
        Psub_norm = Psub_norm./repmat(sum(Psub_norm,2),1,90);
        Psub_norm(isnan(Psub_norm)) = 0;

        % --- reflectance of the OLED stack seen from inside the substrate
        TMF_OLED_bot_p = TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1), ...
            [0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
        TMF_OLED_bot_s = TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1), ...
            [0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
        ROLED = (abs(TMF_OLED_bot_p.r_p).^2 + abs(TMF_OLED_bot_s.r_s).^2)/2;   % unpolarised

        % --- recycling series, following
        %     planar_Sweep22_MLA_JH_full_lambda_modified.m
        %   R_step is the NUMBER OF TERMS j = 1 .. R_step: j = 1 is the direct
        %   pass through the BSDF, j = k is k-1 round trips, each one a BSDF
        %   reflection followed by the OLED mirror.  The reference writes the
        %   k-th term as  Psub_norm * M^(k-1) * BSDF_T_total'  with
        %   M = BSDF_R' .* R_LED; here the vector is carried forward instead of
        %   raising M to a power, which is the same arithmetic and much faster.
        R_step       = 20;
        BSDF_R       = BSDF(180:-1:91,:);
        BSDF_T_total = sum(BSDF(1:90,:));

        % BSDF as an intensity (per solid angle) rather than a power term,
        % used for the angular profile after the MLA
        BSDF_CE = zeros(180,90);
        for a = 1:90
            for b = 1:180
                BSDF_CE(b,a) = BSDF(b,a)*sind(a-0.5)/sind(b-0.5);
            end
        end
        BSDF_T_CE = BSDF_CE(1:90,:);

        ROLED_3    = repmat(ROLED,1,1,90);
        R_matrix_1 = repmat(reshape(BSDF_R',1,90,90),wavelength_num,1,1).*ROLED_3;
        P0         = EQE_sub_matrix;

        Power      = zeros(1,R_step);
        temp_power = zeros(wavelength_num,R_step);
        P_MLA      = zeros(wavelength_num,90,R_step);

        for i = 1:wavelength_num
            M = reshape(R_matrix_1(i,:,:),90,90);
            v = Psub_norm(i,:);                 % normalised, carries the EQE
            w = P_sub(i,:);                     % unnormalised, carries the angular shape
            for j = 1:R_step
                temp_power(i,j) = v*BSDF_T_total';
                P_MLA(i,:,j)    = w*BSDF_T_CE';
                v = v*M;
                w = w*M;
            end
        end
        for j = 1:R_step
            Power(j) = temp_power(:,j)'*P0;
        end

        EQE_sub    = sum(EQE_sub_matrix);
        EQE_MLA    = sum(Power);                % = TotalPower of the reference
        MLA_tail   = Power(end)/EQE_MLA;        % convergence: last term / total
        P_MLA_total = sum(P_MLA,3);             % angular emission profile after the MLA

        %% ===================== end MLA out-coupling =====================

        % fprintf('\nEQE_air = %d \n', EQE_air);
        % fprintf('EQE_subconfined = %d \n', EQE_sub_confined);
        % fprintf('EQE_wg = %d \n', EQE_wg);
        % fprintf('EQE_spp = %d \n', EQE_spp);
        % fprintf('EQE_abs = %d \n', EQE_abs);
        % % fprintf('EQE_sub = %d \n', EQE_sub);
        % % fprintf('FWHM = %d \n', FWHM);
        % disp('-------------------------------')

        k_all=length(d2)*(k1-1)+k2;
        disp_matrix(k_all,:)=[d1(k1),d2(k2),EQE_air,EQE_sub_confined,EQE_wg,EQE_spp,EQE_abs];
        % inplane_matrix(k_all,:)=[d1(k1),d2(k2),U_tot];
        data_matrix(k_all,:)=[d1(k1),d2(k2),EQE_air,EQE_sub_confined,EQE_wg,EQE_spp,EQE_abs,EQE_MLA,EQE_sub];
        fprintf(['d1=%4d d2=%4d | air=%.4f sub=%.4f wg=%.4f spp=%.4f abs=%.4f ' ...
                 '| sum=%.6f res=%+.1e | EQE_MLA=%.4f eta_ext=%.4f tail=%.1e\n'], ...
            d1(k1),d2(k2),EQE_air,EQE_sub_confined,EQE_wg,EQE_spp,EQE_abs, ...
            EQE_closure,EQE_residual,EQE_MLA,EQE_MLA/EQE_sub,MLA_tail);

    end
end
% --- removed: the block below referenced TE/TM-resolved matrices that this
% --- script does not compute (EQE_air_matrix_TE, lambdaemissioneta_sumUtot, ...)
% LEE_out_TE=EQE_air_matrix_TE./(lambdaemissioneta_sumUtot)./sumUtot;
% LEE_out_TMh=EQE_air_matrix_TMh./(lambdaemissioneta_sumUtot)./sumUtot;
% LEE_out_TMv=EQE_air_matrix_TMv./(lambdaemissioneta_sumUtot)./sumUtot;
% LEE_sub_TE=EQE_sub_matrix_TE./(lambdaemissioneta_sumUtot)./sumUtot;
% LEE_sub_TMh=EQE_sub_matrix_TMh./(lambdaemissioneta_sumUtot)./sumUtot;
% LEE_sub_TMv=EQE_sub_matrix_TMv./(lambdaemissioneta_sumUtot)./sumUtot;
% LEE_out=EQE_air_matrix./(lambdaemissioneta_sumUtot)./sumUtot;
% LEE_sub=EQE_sub_matrix./(lambdaemissioneta_sumUtot)./sumUtot;
% 
% U_test=reshape(U_tot,[wavelength_num*u_num,1]);
% w=repmat(wavelength, [u_num 1]);
% uu=zeros(wavelength_num*u_num,1);
% 
% for i=1:u_num
%     uu(wavelength_num*(i-1)+1:wavelength_num*i,1)=u(i);
% end

aa=I_air_total.*sin089;
output = -sum(aa(40:60))/sum(aa(1:90))*EQE_air;

%% Angular-range EQE
% Index convention: sin089 = sind(0:89), index i -> (i-1) deg
% Range [theta1, theta2): indices (theta1+1):(theta2)
aa_total = sum(aa(1:90));
EQE_0_20  = sum(aa(1:20))  / aa_total * EQE_air;   %  0~19 deg
EQE_20_40 = sum(aa(21:40)) / aa_total * EQE_air;   % 20~39 deg
EQE_40_60 = sum(aa(41:60)) / aa_total * EQE_air;   % 40~59 deg
EQE_60_80 = sum(aa(61:80)) / aa_total * EQE_air;   % 60~79 deg

% fprintf('\n--- Angular-range EQE ---\n');
% fprintf('EQE (0-20 deg)  = %.4f\n', EQE_0_20);
% fprintf('EQE (20-40 deg) = %.4f\n', EQE_20_40);
% fprintf('EQE (40-60 deg) = %.4f\n', EQE_40_60);
% fprintf('EQE (60-80 deg) = %.4f\n', EQE_60_80);
% fprintf('Sum (0-80 deg)  = %.4f  |  EQE_air (0-89 deg) = %.4f\n', ...
%     EQE_0_20+EQE_20_40+EQE_40_60+EQE_60_80, EQE_air);
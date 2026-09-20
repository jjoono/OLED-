% Made by WC Lee
% v2 (2026-09-18): power-dissipation bookkeeping revised so that, for eta_rad = 1,
%   EQE_air + EQE_sub_confined + EQE_wg + EQE_spp + EQE_abs = 1  (exactly, on the u grid).
%
%   Changes with respect to Planar_sweep22_preprint.m (marked "v2" below):
%   (1) Purcell factor from an explicit free-space dissipation integral
%       (as in planar_Sweep22_MJKim_250922_gemini5.m), instead of (const+const2)*u_data_num.
%   (2) Absorption = total dissipated power minus power transmitted into the substrate,
%       evaluated per polarization over the substrate-accessible u range. This counts
%       absorption on BOTH sides of the emitter (top reflector and bottom stack).
%       The old code counted only the bottom-side absorption (U_total - U_transmit),
%       which is why the five channels did not add up to 1.
%   (3) WG and SPP are taken from the total dissipation U_tot, per polarization, with
%       polarization-specific cut-offs; the SPP range starts at max(u_sub_max, u_data_num)+1
%       so that it never overlaps the substrate range when n_sub > n_EML.
%   (4) The three u ranges (substrate-accessible / waveguided / evanescent) partition the
%       whole u grid for each polarization, so sum(channels) = sumUtot by construction.
%       No a-posteriori rescaling is applied; the residual is printed as a check.
%   (5) Diagnostics (not part of the five channels): abs split into top/bottom side,
%       sub_confined split into escape-cone recycling absorption and TIR-trapped power,
%       non-radiative remainder (1 - sum) for eta_rad < 1.
%
% clear;
clc;
tic;

load('nk_JH_total.mat')
% load('nk_JHH.mat')
% load('Photopic_400_800.mat')
% Copyright (c) All Rights Reserved.

% dr4: "parasitic absorption" design rule, simplified stack, single wavelength.
%   Sweeps one parameter (MODE) and reports, per point:
%     eta_sub  = power delivered to the substrate (air + substrate-confined)
%     Aprime   = round-trip loss 1 - <R_LED>, R_LED = reflectance of the OLED stack
%                seen from the substrate, flux-weighted over substrate angles, (p+s)/2
%     eta_ext  = p / (p + (1-p) Aprime)            [photon-recycling law, eq. (2)]
%     EQE      = eta_sub * eta_ext
%   Run e.g.: octave --eval "MODE='kito'; SWEEP=0:0.01:0.08; OUTCSV='dr4_kito.csv'; dr4"
if ~exist('MODE','var'),  MODE='kito'; end
if ~exist('SWEEP','var'), SWEEP=0.02; end
if ~exist('OUTCSV','var'),OUTCSV='dr4.csv'; end
if ~exist('NAG','var'),   NAG=0.044; end    % Re(n) of the Ag reflector
if ~exist('KAG','var'),   KAG=3.5;   end    % Im(n) of the Ag reflector
if ~exist('KITO','var'),  KITO=0.02; end
if ~exist('NTCO','var'),  NTCO=1.8;  end
if ~exist('NORG','var'),  NORG=1.8;  end
if ~exist('NE_ETL','var'),NE_ETL=1.6; end   % extraordinary index of the ETL (n_o stays NORG)
if ~exist('KORG','var'),  KORG=0;    end    % extinction of the transport layers
if ~exist('NSUB','var'),  NSUB=1.8;  end
if ~exist('DETL','var'),  DETL=200;  end
if ~exist('DHTL','var'),  DHTL=200;  end
if ~exist('BOT','var'),   BOT='ito'; end
if ~exist('TOPMAT','var'),TOPMAT='ag'; end  % reflector: 'ag' (McPeak) or 'al' (JO)    % bottom electrode: 'ito' (50 nm TCO) OR 'ag' (10 nm thin Ag)
if ~exist('NAGB','var'),  NAGB=0.044; end   % Re(n) of the thin bottom Ag
if ~exist('DAGB','var'),  DAGB=10;   end
if ~exist('DEML','var'),  DEML=20;   end
if ~exist('DTCO','var'),  DTCO=50;   end
if ~exist('DAG','var'),   DAG=100;   end
if ~exist('PESC','var'),  PESC=[0.4 0.6 0.8]; end   % MLA single-pass escape probability
wavelength=550;
wavelength_num=length(wavelength);

emission_spectrum=1;


eta_rad=1;
horizontal_dipole_ratio=2/3;
bottom_air_refractive_index=ones(wavelength_num,1);


%% JO structure %%
% The bottom electrode is EITHER a 50-nm TCO OR a 10-nm thin Ag, never both.
if strcmp(TOPMAT,'al'), nAgTop=material.l_Al_JO(151); else, nAgTop=material.l_Ag_McPeak(151); end
if strcmp(BOT,'ag'), nBot=NAGB+KAG*1i; dBot=DAGB; else, nBot=NTCO+KITO*1i; dBot=DTCO; end
no_bar=[1 nAgTop NORG+KORG*1i NORG NORG+KORG*1i nBot NSUB];
ne_bar=no_bar;
ne_bar(3)=NE_ETL+KORG*1i;   % uniaxial ETL: n_o = NORG, n_e = NE_ETL; everything else isotropic

d1= SWEEP(:)'; Nd1=length(d1);
d2= 0; Nd2=1;
data_matrix = zeros(Nd1*Nd2,12);
disp_matrix = zeros(Nd1*Nd2,6);
diag_matrix = zeros(Nd1*Nd2,9); ne_rows=zeros(Nd1*Nd2,1); extra_matrix=zeros(Nd1*Nd2,8+2*length(PESC));   % v2: [d1, abs_top, abs_bottom, sub_TIR, sub_recycle_abs, sum5, nonrad, residual, Purcell]
for k1=1:length(d1)
    for k2=1:length(d2)
        thickness=[DAG DETL DEML DHTL dBot];   % Ag(top), ETL, EML, HTL, bottom electrode
        switch MODE
          case 'kito', no_bar(6)=NTCO+d1(k1)*1i;      % extinction of the TCO electrode
          case 'nagb', no_bar(6)=d1(k1)+KAG*1i;       % Re(n) of the thin Ag electrode
          case 'dbot', thickness(5)=d1(k1);           % thickness of the bottom electrode
          case 'nag',  no_bar(2)=d1(k1)+KAG*1i;
          case 'detl', thickness(2)=d1(k1);
          case 'dctl', thickness(2)=d1(k1); thickness(4)=d1(k1);
          case 'netl', NE_ETL=d1(k1);                 % extraordinary index of the ETL
          case 'nsub', NSUB=d1(k1); no_bar(end)=NSUB;
          case 'korg', no_bar(3)=NORG+d1(k1)*1i; no_bar(5)=NORG+d1(k1)*1i;
          otherwise,   error('unknown MODE %s', MODE);
        end
        ne_bar=no_bar; ne_bar(3)=NE_ETL+KORG*1i;

        EML_position=4;

        %         z0=d_slice/2; % emitting position in EML
        z0=DEML/2;   % dipole at the centre of the EML
        if exist('UNUM','var'), u_data_num=UNUM; else, u_data_num=997; end
        if exist('MAXU','var'), max_u=MAXU; else, max_u=3; end

        %%

        sin089=sind(0:89);

        layer_num=size(no_bar,2);

        u=[(0:u_data_num-1)/u_data_num (u_data_num+1:u_data_num*max_u)/u_data_num];

        u_num=length(u);

        TMF_bottom=TMF_birefringence_whole(no_bar(:,EML_position:layer_num),ne_bar(:,EML_position:layer_num),[thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0],u,wavelength);
        TMF_top=TMF_birefringence_whole(no_bar(:,EML_position:-1:1),ne_bar(:,EML_position:-1:1),[z0 thickness(EML_position-2:-1:1) 0],u,wavelength);

        K_p_v=3/4*real(ne_bar(:,EML_position)./no_bar(:,EML_position)*(u.^2./sqrt(1-u.^2)).*(1+TMF_bottom.r_p).*(1+TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
        K_p_h=real(3./(6*(no_bar(:,EML_position)./ne_bar(:,EML_position)).^2+2)*sqrt(1-u.^2).*(1-TMF_bottom.r_p).*(1-TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
        K_s_h=real(3./(2*(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2+6)*(1./sqrt(1-u.^2)).*(1+TMF_bottom.r_s).*(1+TMF_top.r_s)./(1-TMF_bottom.r_s.*TMF_top.r_s));

        K_p_v2=K_p_v;
        K_p_h2=K_p_h;
        K_s_h2=K_s_h;

        K_p_v2_total=K_p_v;
        K_p_h2_total=K_p_h;
        K_s_h2_total=K_s_h;

        K_p_v3=K_p_v;
        K_p_h3=K_p_h;
        K_s_h3=K_s_h;

        u_sub_max_p=zeros(wavelength_num,1);
        u_sub_max_s=zeros(wavelength_num,1);

        u_air_max_p=zeros(wavelength_num,1);
        u_air_max_s=zeros(wavelength_num,1);

        for i=1:wavelength_num

            if ne_bar(i,layer_num)>ne_bar(i,EML_position)

                u_sub_max_p(i)=ceil(u_data_num*ne_bar(i,layer_num)/ne_bar(i,EML_position))-1;

            else

                u_sub_max_p(i)=ceil(u_data_num*ne_bar(i,layer_num)/ne_bar(i,EML_position));

            end

            exp_phase=ones(1,u_sub_max_p(i));

            if u_sub_max_p(i)>u_data_num

                exp_phase(u_data_num+1:u_sub_max_p(i))=exp((-4*pi*no_bar(i,EML_position)*sqrt(u(u_data_num+1:u_sub_max_p(i)).^2-1)*(thickness(EML_position-1)-z0))/wavelength(i));

            end

            K_p_v2(i,1:u_sub_max_p(i))=3/8*ne_bar(i,EML_position)*no_bar(i,layer_num)/no_bar(i,EML_position)^2*sqrt(1-(ne_bar(i,EML_position)*u(1:u_sub_max_p(i))/ne_bar(i,layer_num)).^2).*exp_phase.*u(1:u_sub_max_p(i)).^2.*abs((1+TMF_top.r_p(i,1:u_sub_max_p(i))).*TMF_bottom.t_p(i,1:u_sub_max_p(i))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2./abs(1-u(1:u_sub_max_p(i)).^2);
            K_p_h2(i,1:u_sub_max_p(i))=3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2*(1-(ne_bar(i,EML_position)*u(1:u_sub_max_p(i))/ne_bar(i,layer_num)).^2)).*exp_phase.*abs((1-TMF_top.r_p(i,1:u_sub_max_p(i))).*TMF_bottom.t_p(i,1:u_sub_max_p(i))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2/(12*(no_bar(i,EML_position)/ne_bar(i,EML_position))^2+4);

            if no_bar(i,layer_num)>no_bar(i,EML_position)

                u_sub_max_s(i)=ceil(u_data_num*no_bar(i,layer_num)/no_bar(i,EML_position))-1;

            else

                u_sub_max_s(i)=ceil(u_data_num*no_bar(i,layer_num)/no_bar(i,EML_position));

            end

            exp_phase=ones(1,u_sub_max_s(i));

            if u_sub_max_s(i)>u_data_num

                exp_phase(u_data_num+1:u_sub_max_s(i))=exp((-4*pi*no_bar(i,EML_position)*sqrt(u(u_data_num+1:u_sub_max_s(i)).^2-1)*(thickness(EML_position-1)-z0))/wavelength(i));

            end

            K_s_h2(i,1:u_sub_max_s(i))=3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2-u(1:u_sub_max_s(i)).^2).*exp_phase.*abs((1+TMF_top.r_s(i,1:u_sub_max_s(i))).*TMF_bottom.t_s(i,1:u_sub_max_s(i))./(1-TMF_bottom.r_s(i,1:u_sub_max_s(i)).*TMF_top.r_s(i,1:u_sub_max_s(i)))).^2./((4*(ne_bar(i,EML_position)/no_bar(i,EML_position))^2+12)*abs(1-u(1:u_sub_max_s(i)).^2));

            % net power flux into the bottom half-stack (transmitted + absorbed there); prefactors are half of K_p_v/K_p_h/K_s_h
            K_p_v2_total(i,1:u_sub_max_p(i))=3/8*(u(1:u_sub_max_p(i)).^2).*real((1+TMF_bottom.r_p(i,1:u_sub_max_p(i))).*(1-conj(TMF_bottom.r_p(i,1:u_sub_max_p(i))))./sqrt(1-u(1:u_sub_max_p(i)).^2)).*abs((1+TMF_top.r_p(i,1:u_sub_max_p(i)))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2;
            K_p_h2_total(i,1:u_sub_max_p(i))=3*real((1-TMF_bottom.r_p(i,1:u_sub_max_p(i))).*(1+conj(TMF_bottom.r_p(i,1:u_sub_max_p(i)))).*sqrt(1-u(1:u_sub_max_p(i)).^2)).*abs((1-TMF_top.r_p(i,1:u_sub_max_p(i)))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2/(12*(no_bar(i,EML_position)/ne_bar(i,EML_position))^2+4);
            K_s_h2_total(i,1:u_sub_max_s(i))=3*real((1+TMF_bottom.r_s(i,1:u_sub_max_s(i))).*(1-conj(TMF_bottom.r_s(i,1:u_sub_max_s(i))))./sqrt(1-u(1:u_sub_max_s(i)).^2)).*abs((1+TMF_top.r_s(i,1:u_sub_max_s(i)))./(1-TMF_bottom.r_s(i,1:u_sub_max_s(i)).*TMF_top.r_s(i,1:u_sub_max_s(i)))).^2/(4*(ne_bar(i,EML_position)/no_bar(i,EML_position))^2+12);

            K_p_v3(i,:)=K_p_v2(i,:);
            K_p_h3(i,:)=K_p_h2(i,:);
            K_s_h3(i,:)=K_s_h2(i,:);

            if bottom_air_refractive_index(i)>ne_bar(i,EML_position)

                u_air_max_p(i)=min(u_sub_max_p(i),ceil(bottom_air_refractive_index(i)*u_data_num/ne_bar(i,EML_position))-1);

            else

                u_air_max_p(i)=min(u_sub_max_p(i),ceil(bottom_air_refractive_index(i)*u_data_num/ne_bar(i,EML_position)));

            end

            TMF_OLED_bottom_p=TMF_birefringence_whole_p(no_bar(i,layer_num:-1:1),ne_bar(i,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(i,EML_position)*u(1:u_air_max_p(i)),wavelength(i));

            R_p_bottom=abs(TMF_OLED_bottom_p.r_p).^2;

            cos_theta_sub=sqrt(1-(ne_bar(i,EML_position)*u(1:u_air_max_p(i))/ne_bar(i,layer_num)).^2);
            cos_theta_air=sqrt(1-(ne_bar(i,EML_position)*u(1:u_air_max_p(i))/bottom_air_refractive_index(i)).^2);

            r_p=(bottom_air_refractive_index(i)*cos_theta_sub-no_bar(i,layer_num)*cos_theta_air)./(bottom_air_refractive_index(i)*cos_theta_sub+no_bar(i,layer_num)*cos_theta_air);

            R_sub_air_bottom_p=abs(r_p).^2;
            T_sub_air_bottom_p=1-R_sub_air_bottom_p;

            if bottom_air_refractive_index(i)>no_bar(i,EML_position)

                u_air_max_s(i)=min(u_sub_max_s(i),ceil(bottom_air_refractive_index(i)*u_data_num/no_bar(i,EML_position))-1);

            else

                u_air_max_s(i)=min(u_sub_max_s(i),ceil(bottom_air_refractive_index(i)*u_data_num/no_bar(i,EML_position)));

            end

            TMF_OLED_bottom_s=TMF_birefringence_whole_s(no_bar(i,layer_num:-1:1),ne_bar(i,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(i,EML_position)*u(1:u_air_max_s(i)),wavelength(i));

            R_s_bottom=abs(TMF_OLED_bottom_s.r_s).^2;

            cos_theta_sub=sqrt(1-(no_bar(i,EML_position)*u(1:u_air_max_s(i))/no_bar(i,layer_num)).^2);
            cos_theta_air=sqrt(1-(no_bar(i,EML_position)*u(1:u_air_max_s(i))/bottom_air_refractive_index(i)).^2);

            r_s=(no_bar(i,layer_num)*cos_theta_sub-bottom_air_refractive_index(i)*cos_theta_air)./(no_bar(i,layer_num)*cos_theta_sub+bottom_air_refractive_index(i)*cos_theta_air);

            R_sub_air_bottom_s=abs(r_s).^2;
            T_sub_air_bottom_s=1-R_sub_air_bottom_s;

            K_p_v3(i,1:u_air_max_p(i))=K_p_v2(i,1:u_air_max_p(i)).*T_sub_air_bottom_p./(1-R_p_bottom.*R_sub_air_bottom_p);
            K_p_h3(i,1:u_air_max_p(i))=K_p_h2(i,1:u_air_max_p(i)).*T_sub_air_bottom_p./(1-R_p_bottom.*R_sub_air_bottom_p);
            K_s_h3(i,1:u_air_max_s(i))=K_s_h2(i,1:u_air_max_s(i)).*T_sub_air_bottom_s./(1-R_s_bottom.*R_sub_air_bottom_s);

        end

        const=(1-horizontal_dipole_ratio)*ne_bar(:,EML_position)./(wavelength.^4);
        const2=horizontal_dipole_ratio*no_bar(:,EML_position).*(3+(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2)./(4.*wavelength.^4);
        const3=const*u;
        const4=const2*u;

        U_tot=2*(const3.*K_p_v+const4.*(K_p_h+K_s_h));

        % v2: polarization-resolved total dissipation (U_tot = U_tot_p + U_tot_s)
        U_tot_p=2*(const3.*K_p_v+const4.*K_p_h);
        U_tot_s=2*const4.*K_s_h;

        U_bottom_transmit_p=2*(const3.*K_p_v2+const4.*K_p_h2);
        U_bottom_transmit_ph=2*(const4.*K_p_h2);
        U_bottom_transmit_pv=2*(const3.*K_p_v2);
        U_bottom_transmit_s=2*const4.*K_s_h2;
        U_bottom_transmit=U_bottom_transmit_p+U_bottom_transmit_s;

        U_bottom_transmit_total_p=2*(const3.*K_p_v2_total+const4.*K_p_h2_total);
        U_bottom_transmit_total_s=2*const4.*K_s_h2_total;
        U_bottom_transmit_total=U_bottom_transmit_total_p+U_bottom_transmit_total_s;

        U_bottom_transmit_thick_p=2*(const3.*K_p_v3+const4.*K_p_h3);
        U_bottom_transmit_thick_ph=2*(const4.*K_p_h3);
        U_bottom_transmit_thick_pv=2*(const3.*K_p_v3);
        U_bottom_transmit_thick_s=2*const4.*K_s_h3;
        U_bottom_transmit_thick=U_bottom_transmit_thick_p+U_bottom_transmit_thick_s;

        % v2: free-space dissipation on the same u grid (r = 0, u < 1) for the Purcell factor
        u_free=u(1:u_data_num);
        K_p_v_free=3/4*real(ne_bar(:,EML_position)./no_bar(:,EML_position)*(u_free.^2./sqrt(1-u_free.^2)));
        K_p_h_free=real(3./(6*(no_bar(:,EML_position)./ne_bar(:,EML_position)).^2+2)*sqrt(1-u_free.^2));
        K_s_h_free=real(3./(2*(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2+6)*(1./sqrt(1-u_free.^2)));
        U_free=2*((const*u_free).*K_p_v_free+(const2*u_free).*(K_p_h_free+K_s_h_free));
        P_free=sum(U_free,2);

        const3=repmat(const,1,u_num);
        const4=repmat(const2,1,u_num);

        K_bottom_transmit_p=const3.*K_p_v2+const4.*K_p_h2;
        K_bottom_transmit_s=const4.*K_s_h2;
        K_bottom_transmit=K_bottom_transmit_p+K_bottom_transmit_s;

        K_bottom_transmit_thick_p=const3.*K_p_v3+const4.*K_p_h3;
        K_bottom_transmit_thick_s=const4.*K_s_h3;
        K_bottom_transmit_thick=K_bottom_transmit_thick_p+K_bottom_transmit_thick_s;

        Power_ratio_air_matrix=zeros(wavelength_num,1);
        Power_ratio_air2_matrix=zeros(wavelength_num,1);
        Power_ratio_sub_matrix=zeros(wavelength_num,1);
        Power_ratio_abs_matrix=zeros(wavelength_num,1);
        Power_ratio_wg_matrix=zeros(wavelength_num,1);
        Power_ratio_spp_matrix=zeros(wavelength_num,1);

        EQE_air_matrix=zeros(wavelength_num,1);
        EQE_air2_matrix=zeros(wavelength_num,1);
        EQE_sub_matrix=zeros(wavelength_num,1);
        %%
        EQE_air_matrix_TE=zeros(wavelength_num,1);
        EQE_air_matrix_TMh=zeros(wavelength_num,1);
        EQE_air_matrix_TMv=zeros(wavelength_num,1);
        EQE_sub_matrix_TE=zeros(wavelength_num,1);
        EQE_sub_matrix_TMh=zeros(wavelength_num,1);
        EQE_sub_matrix_TMv=zeros(wavelength_num,1);
        %%
        EQE_abs_matrix=zeros(wavelength_num,1);
        EQE_wg_matrix=zeros(wavelength_num,1);
        EQE_spp_matrix=zeros(wavelength_num,1);
        % v2 diagnostics
        EQE_abs_top_matrix=zeros(wavelength_num,1);
        EQE_abs_bottom_matrix=zeros(wavelength_num,1);
        EQE_closure_matrix=zeros(wavelength_num,1);

        sumUtot=sum(U_tot,2);

        % v2: Purcell factor = total dissipation / free-space dissipation (explicit integral)
        Purcell_factor=sumUtot./P_free;
        eta_eff=eta_rad*Purcell_factor./(1-eta_rad+eta_rad*Purcell_factor);

        emission_spectrum=emission_spectrum/sum(emission_spectrum);

        emissioneta_sumUtot=emission_spectrum.*eta_eff./sumUtot;
        lambdaemissioneta_sumUtot=wavelength.*emissioneta_sumUtot/sum(wavelength.*emission_spectrum);

        P_air_p=zeros(wavelength_num,90);
        P_air_s=zeros(wavelength_num,90);

        P_sub_p=zeros(wavelength_num,90);
        P_sub_s=zeros(wavelength_num,90);

        const3=pi*(const+const2);
        const=bottom_air_refractive_index./ne_bar(:,EML_position);
        const2=bottom_air_refractive_index./no_bar(:,EML_position);

        for i=1:wavelength_num

            % v2: u-ranges per polarization. For each polarization the three ranges
            % (substrate-accessible | waveguided | evanescent) partition 1:u_num.
            ip_air=1:u_air_max_p(i);                       is_air=1:u_air_max_s(i);
            ip_sub=1:u_sub_max_p(i);                       is_sub=1:u_sub_max_s(i);
            ip_wg =u_sub_max_p(i)+1:u_data_num;            is_wg =u_sub_max_s(i)+1:u_data_num;            % empty if n_sub >= n_EML
            ip_spp=max(u_sub_max_p(i),u_data_num)+1:u_num; is_spp=max(u_sub_max_s(i),u_data_num)+1:u_num;

            % v2: channel powers (in units of the dissipation spectrum; normalized by sumUtot below)
            P_air   =sum(U_bottom_transmit_thick_p(i,ip_air))+sum(U_bottom_transmit_thick_s(i,is_air));   % outcoupled to air (thick substrate, multiple reflections)
            P_air2  =sum(U_bottom_transmit_p(i,ip_air))+sum(U_bottom_transmit_s(i,is_air));               % delivered to substrate inside the air escape cone
            P_sub   =sum(U_bottom_transmit_p(i,ip_sub))+sum(U_bottom_transmit_s(i,is_sub));               % delivered to substrate (all angles)
            P_down  =sum(U_bottom_transmit_total_p(i,ip_sub))+sum(U_bottom_transmit_total_s(i,is_sub));   % net flux into the bottom half-stack
            P_totsub=sum(U_tot_p(i,ip_sub))+sum(U_tot_s(i,is_sub));                                        % total dissipation in the substrate-accessible range
            P_abs   =P_totsub-P_sub;      % absorbed anywhere (top reflector + bottom stack) -> THE FIX
            P_abs_b =P_down-P_sub;        % ... of which in the bottom stack (TCO/CTL)  [= old EQE_abs]
            P_abs_t =P_totsub-P_down;     % ... of which on the top/reflector side
            P_wg    =sum(U_tot_p(i,ip_wg))+sum(U_tot_s(i,is_wg));
            P_spp   =sum(U_tot_p(i,ip_spp))+sum(U_tot_s(i,is_spp));

            w_pow=emissioneta_sumUtot(i);          % power weighting
            w_ph =lambdaemissioneta_sumUtot(i);    % photon-number weighting

            Power_ratio_air_matrix(i) =w_pow*P_air;
            Power_ratio_air2_matrix(i)=w_pow*P_air2;
            Power_ratio_sub_matrix(i) =w_pow*P_sub;
            Power_ratio_abs_matrix(i) =w_pow*P_abs;
            Power_ratio_wg_matrix(i)  =w_pow*P_wg;
            Power_ratio_spp_matrix(i) =w_pow*P_spp;

            EQE_air_matrix(i) =w_ph*P_air;
            EQE_air2_matrix(i)=w_ph*P_air2;
            EQE_sub_matrix(i) =w_ph*P_sub;
            EQE_abs_matrix(i) =w_ph*P_abs;
            EQE_wg_matrix(i)  =w_ph*P_wg;
            EQE_spp_matrix(i) =w_ph*P_spp;
            EQE_abs_top_matrix(i)   =w_ph*P_abs_t;
            EQE_abs_bottom_matrix(i)=w_ph*P_abs_b;
            EQE_closure_matrix(i)   =w_ph*sumUtot(i);   % = w_ph*(P_totsub+P_wg+P_spp): what the five channels must add up to
            %%
            EQE_air_matrix_TE(i)=w_ph*(sum(U_bottom_transmit_thick_s(i,is_air)));
            EQE_air_matrix_TMh(i)=w_ph*(sum(U_bottom_transmit_thick_ph(i,ip_air)));
            EQE_air_matrix_TMv(i)=w_ph*(sum(U_bottom_transmit_thick_pv(i,ip_air)));
            %%
            EQE_sub_matrix_TE(i)=w_ph*(sum(U_bottom_transmit_s(i,is_sub)));
            EQE_sub_matrix_TMh(i)=w_ph*(sum(U_bottom_transmit_ph(i,ip_sub)));
            EQE_sub_matrix_TMv(i)=w_ph*(sum(U_bottom_transmit_pv(i,ip_sub)));
            %%

            P_air_p(i,:)=const(i)*spline(ne_bar(i,EML_position)*u(1:u_air_max_p(i)),sqrt(const(i)^2-u(1:u_air_max_p(i)).^2).*K_bottom_transmit_thick_p(i,1:u_air_max_p(i)),bottom_air_refractive_index(i)*sin089)/const3(i);
            P_air_s(i,:)=const2(i)*spline(no_bar(i,EML_position)*u(1:u_air_max_s(i)),sqrt(const2(i)^2-u(1:u_air_max_s(i)).^2).*K_bottom_transmit_thick_s(i,1:u_air_max_s(i)),bottom_air_refractive_index(i)*sin089)/const3(i);

            if bottom_air_refractive_index(i)>ne_bar(i,layer_num)

                P_air_p(i,ceil(asind(ne_bar(i,layer_num)/bottom_air_refractive_index(i)))+1:90)=0;

            end

            if bottom_air_refractive_index(i)>no_bar(i,layer_num)

                P_air_s(i,ceil(asind(no_bar(i,layer_num)/bottom_air_refractive_index(i)))+1:90)=0;

            end

            P_sub_p(i,:)=(ne_bar(i,layer_num)/ne_bar(i,EML_position))*spline(ne_bar(i,EML_position)*u(1:u_sub_max_p(i)),sqrt((ne_bar(i,layer_num)/ne_bar(i,EML_position))^2-u(1:u_sub_max_p(i)).^2).*K_bottom_transmit_p(i,1:u_sub_max_p(i)),ne_bar(i,layer_num)*sin089)/const3(i);
            P_sub_s(i,:)=(no_bar(i,layer_num)/no_bar(i,EML_position))*spline(no_bar(i,EML_position)*u(1:u_sub_max_s(i)),sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2-u(1:u_sub_max_s(i)).^2).*K_bottom_transmit_s(i,1:u_sub_max_s(i)),no_bar(i,layer_num)*sin089)/const3(i);

        end

        Power_ratio_air=sum(Power_ratio_air_matrix);
        Power_ratio_air2=sum(Power_ratio_air2_matrix);
        Power_ratio_sub=sum(Power_ratio_sub_matrix);
        Power_ratio_sub_confined=Power_ratio_sub-Power_ratio_air;
        Power_ratio_abs=sum(Power_ratio_abs_matrix);
        Power_ratio_wg=sum(Power_ratio_wg_matrix);
        Power_ratio_spp=sum(Power_ratio_spp_matrix);

        EQE_air=sum(EQE_air_matrix);
        EQE_air2=sum(EQE_air2_matrix);
        EQE_sub=sum(EQE_sub_matrix);
        EQE_sub_confined=EQE_sub-EQE_air;
        EQE_abs=sum(EQE_abs_matrix);
        EQE_wg=sum(EQE_wg_matrix);
        EQE_spp=sum(EQE_spp_matrix);

        % v2 diagnostics
        EQE_abs_top=sum(EQE_abs_top_matrix);          % reflector-side absorption of substrate-accessible light
        EQE_abs_bottom=sum(EQE_abs_bottom_matrix);    % bottom-stack absorption (what the old code reported as EQE_abs)
        EQE_sub_TIR=EQE_sub-EQE_air2;                 % delivered to substrate outside the escape cone (TIR-trapped)
        EQE_sub_recycle_abs=EQE_air2-EQE_air;         % escape-cone light absorbed in the stack during substrate multiple reflections
        EQE_sum5=EQE_air+EQE_sub_confined+EQE_wg+EQE_spp+EQE_abs;
        EQE_radiative=sum(EQE_closure_matrix);        % = 1 for eta_rad = 1; photon-weighted effective radiative yield otherwise
        EQE_nonrad=1-EQE_radiative;
        EQE_residual=EQE_sum5-EQE_radiative;          % should be ~1e-16

        P_air=P_air_p+P_air_s;
        P_sub=P_sub_p+P_sub_s;
        %         Purcell_factor=ones(401,1);
        I_air=P_air.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);
        I_air_p=P_air_p.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);

        I_air_total=sum(I_air);
        I_air_sum=sum(I_air);
        I_air_total=I_air_total/I_air_total(1);

        I_sub=P_sub.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);
        I_sub_p=P_sub_p.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);
        I_sub_total=sum(I_sub);
        I_sub_sum=sum(I_sub);
        I_sub_total=I_sub_total/I_sub_total(1);

        EQE_factor_air=pi*sum(I_air_total.*sin089)/90;
        EQE_factor_sub=pi*sum(I_sub_total.*sin089)/90;
        for i=1:wavelength_num
            spec_lambda(i,1)=emission_spectrum(i,1)/((i+wavelength(1)-1)*10^(-9));
        end

        I_FWHM=I_air(:,1);
        I_FWHM=I_FWHM/max(I_FWHM); %normalized spectrum
        FWHM=sum(I_FWHM>=0.5);

        fprintf('\nd1 = %g nm\n', d1(k1));
        fprintf('EQE_air = %.6f \n', EQE_air);
        fprintf('EQE_subconfined = %.6f   (TIR-trapped %.6f + escape-cone recycling absorption %.6f)\n', EQE_sub_confined, EQE_sub_TIR, EQE_sub_recycle_abs);
        fprintf('EQE_wg = %.6f \n', EQE_wg);
        fprintf('EQE_spp = %.6f \n', EQE_spp);
        fprintf('EQE_abs = %.6f   (top/reflector side %.6f + bottom stack %.6f)\n', EQE_abs, EQE_abs_top, EQE_abs_bottom);
        fprintf('EQE_sub = %.6f \n', EQE_sub);
        fprintf('Purcell factor = %.4f,  eta_eff = %.4f\n', Purcell_factor(1), eta_eff(1));
        fprintf('Sum of five channels = %.6f   (expected %.6f; residual %.2e; non-radiative %.6f)\n', EQE_sum5, EQE_radiative, EQE_residual, EQE_nonrad);
        fprintf('FWHM = %d \n', FWHM);
        disp('-------------------------------')

        % ---- round-trip loss A' : reflectance of the OLED stack seen from the substrate ----
        th_e=linspace(0,pi/2,1801); th_c=(th_e(1:end-1)+th_e(2:end))/2;
        kap=NSUB*sin(th_c);
        Rp_=abs(TMF_birefringence_whole_p(no_bar(layer_num:-1:1),ne_bar(layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],kap,wavelength).r_p).^2;
        Rs_=abs(TMF_birefringence_whole_s(no_bar(layer_num:-1:1),ne_bar(layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],kap,wavelength).r_s).^2;
        wgt=cos(th_c).*sin(th_c);
        R_LED=sum(0.5*(Rp_+Rs_).*wgt)/sum(wgt);   % flux-weighted isotropic (randomising MLA)
        Aprime=1-R_LED;
        eta_sub=EQE_air+EQE_sub_confined;
        eta_ext=PESC./(PESC+(1-PESC)*Aprime);
        fprintf('  eta_sub = %.4f | A'' = %.4f (R_LED = %.4f) | eta_ext = %s | EQE = %s\n', ...
            eta_sub, Aprime, R_LED, mat2str(round(eta_ext*1e4)/1e4), mat2str(round(eta_sub*eta_ext*1e4)/1e4));

        k_all=length(d2)*(k1-1)+k2;
        disp_matrix(k_all,:)=[d1(k1),EQE_air,EQE_sub_confined,EQE_wg,EQE_spp,EQE_abs]; ne_rows(k_all)=d2(k2);
        extra_matrix(k_all,:)=[d1(k1),eta_sub,Aprime,EQE_wg,EQE_spp,EQE_abs,EQE_abs_top,EQE_abs_bottom,eta_ext,eta_sub*eta_ext];
        diag_matrix(k_all,:)=[d1(k1),EQE_abs_top,EQE_abs_bottom,EQE_sub_TIR,EQE_sub_recycle_abs,EQE_sum5,EQE_nonrad,EQE_residual,Purcell_factor(1)];

%         data_matrix(k_all,:)=[d1(k1),d2(k2),FWHM,EQE_air,EQE_sub_confined,EQE_wg,EQE_spp,EQE_abs,EQE_sub,EQE_factor_air,CE,I_sub_sum_30];

    end
end
LEE_out_TE=EQE_air_matrix_TE./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_out_TMh=EQE_air_matrix_TMh./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_out_TMv=EQE_air_matrix_TMv./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_sub_TE=EQE_sub_matrix_TE./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_sub_TMh=EQE_sub_matrix_TMh./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_sub_TMv=EQE_sub_matrix_TMv./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_out=EQE_air_matrix./(lambdaemissioneta_sumUtot)./sumUtot;
LEE_sub=EQE_sub_matrix./(lambdaemissioneta_sumUtot)./sumUtot;

U_test=reshape(U_tot,[wavelength_num*u_num,1]);
w=repmat(wavelength, [u_num 1]);
uu=zeros(wavelength_num*u_num,1);

for i=1:u_num
    uu(wavelength_num*(i-1)+1:wavelength_num*i,1)=u(i);
end
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

fprintf('\n--- Angular-range EQE ---\n');
fprintf('EQE (0-20 deg)  = %.4f\n', EQE_0_20);
fprintf('EQE (20-40 deg) = %.4f\n', EQE_20_40);
fprintf('EQE (40-60 deg) = %.4f\n', EQE_40_60);
fprintf('EQE (60-80 deg) = %.4f\n', EQE_60_80);
fprintf('Sum (0-80 deg)  = %.4f  |  EQE_air (0-89 deg) = %.4f\n', ...
    EQE_0_20+EQE_20_40+EQE_40_60+EQE_60_80, EQE_air);
csvwrite(OUTCSV,extra_matrix); toc;   % columns: param, eta_sub, Aprime, wg, spp, abs, eta_ext(per p), EQE(per p)

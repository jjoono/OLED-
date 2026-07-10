% Made by WC Lee
% clear;
clc;
tic;

load('nk_JH33.mat')
% load('nk_JHH.mat')
load('Photopic_400_800.mat')
V_301=V_401(1:301,1);
% Copyright �� All Rights Reserved.

wavelength=(400:800).';
wavelength_num=length(wavelength);

% emission_spectrum=I_TCTA_B3PyMPM;
% emission_spectrum=spectrum.FAPbI3_ccho;
% emission_spectrum=spectrum.l_I_Irmphmq2tmd_measure_JH;%eta_rad=0.96, hdr=0.82
% emission_spectrum=spectrum.I_Irppy2acac;
% emission_spectrum=spectrum.InP_HCCho;
% emission_spectrum=spectrum.l_I_IrMDQ2acac; %eta_rad=0.82, hdr=0.76
emission_spectrum=spectrum.l_I_Irdmppyph2tmd; %% JOSong, eta_rad=0.98, hdr=0.865
% emission_spectrum=spectrum.MDQ2acac_andreas;
% emission_spectrum=spectrum.I_Irdmppyphtmd;
% emission_spectrum=spectrum.PNEL_temp;
% emission_spectrum=emission_spectrum(51:71,:);
% emission_spectrum=spectrum.I_Irppy2acac;
% emission_spectrum=emission_spectrum(1:301);
% emission_spectrum=spectrum.Perov_TWLEE;
% emission_spectrum=spectrum.Irppy2acac_ETRI;
emission_spectrum=emission_spectrum(1:401,1);
emission_spectrum=emission_spectrum/sum(emission_spectrum);

eta_rad=1;
horizontal_dipole_ratio=0.865;
bottom_air_refractive_index=ones(wavelength_num,1);

% no_bar=[ones(wavelength_num,1) material.Al 1.8*ones(301,3) material.ITO 1.5*ones(301,1)];
% ne_bar=no_bar;

% no_bar=[ones(wavelength_num,1) material.Al material.B3PyMPM_o material.TCTA_B3PyMPM_o material.TCTA_o material.TAPC_o material.ITO material.Glass_fit];
% ne_bar=[ones(wavelength_num,1) material.Al material.B3PyMPM_e material.TCTA_B3PyMPM_e material.TCTA_e material.TAPC_e material.ITO material.Glass_fit];

% no_bar=[ones(wavelength_num,1) material.l_Al   material.l_MoO3_JH  material.l_TAPC material.l_TCTA_o material.l_TCTA_B3_o material.l_B3_o material.l_Al material.l_Ag_12nm_test material.l_MoO3_JH  material.l_LiF   material.l_Al2O3 material.l_ParyC material.l_Al2O3 material.l_ParyC material.l_Glass];
% ne_bar=[ones(wavelength_num,1) material.l_Al   material.l_MoO3_JH  material.l_TAPC material.l_TCTA_e material.l_TCTA_B3_e_nok material.l_B3_e material.l_Al material.l_Ag_12nm_test material.l_MoO3_JH  material.l_LiF   material.l_Al2O3 material.l_ParyC material.l_Al2O3 material.l_ParyC material.l_Glass];

% no_bar=[ones(wavelength_num,1) material.l_Al   material.l_MoO3_JH  material.l_TAPC material.l_TCTA_o material.l_TCTA_B3_o material.l_B3_o material.l_Al material.l_Ag_12nm_test material.l_MoO3_JH material.l_LiF   material.l_Glass];
% ne_bar=[ones(wavelength_num,1) material.l_Al   material.l_MoO3_JH  material.l_TAPC material.l_TCTA_e material.l_TCTA_B3_e_nok material.l_B3_e material.l_Al material.l_Ag_12nm_test material.l_MoO3_JH material.l_LiF  material.l_Glass];

% no_bar=[ones(wavelength_num,1) material.l_Al material.l_B3_o material.l_TCTA_B3_o material.l_TCTA_o material.l_TAPC material.l_IZO material.l_Glass];
% ne_bar=[ones(wavelength_num,1) material.l_Al material.l_B3_e material.l_TCTA_B3_o material.l_TCTA_e material.l_TAPC material.l_IZO material.l_Glass];

% no_bar=[ones(wavelength_num,1) material.l_Al material.l_TCTA_o material.l_NPB_nok    material.l_B3_o  material.l_NPB_nok material.l_Glass];
% ne_bar=[ones(wavelength_num,1) material.l_Al material.l_TCTA_o material.l_NPB_nok    material.l_B3_o  material.l_NPB_nok material.l_Glass];

% no_bar=[material.l_Air material.l_Al material.l_TCTA_o material.l_NPB material.l_NPB_nok material.l_NPB material.l_NPB_nok  material.l_NPB_nok material.l_Glass];
% ne_bar=[material.l_Air material.l_Al material.l_TCTA_o material.l_NPB material.l_NPB_nok material.l_NPB material.l_NPB_nok  material.l_NPB_nok material.l_Glass];

% no_bar=[material.l_Air material.l_Al material.l_TCTA_o material.l_185 material.l_NPB_nok  material.l_NPB_nok material.l_Glass];
% ne_bar=[material.l_Air material.l_Al material.l_TCTA_o material.l_185 material.l_NPB_nok  material.l_NPB_nok material.l_Glass];

% no_bar=[material.Air material.Ag_NIR material.MoO3_arb material.TFB_NIR real(material.FAPbI3_ccho) material.ZnO_PEIE_NIR real(material.IZO_NIR) 1.5*ones(wavelength_num,1)];
% ne_bar=[material.Air material.Ag_NIR material.MoO3_arb material.TFB_NIR real(material.FAPbI3_ccho) material.ZnO_PEIE_NIR real(material.IZO_NIR) 1.5*ones(wavelength_num,1)];

% no_bar=[ones(wavelength_num,1) material.l_Al  material.l_HATCN material.l_TCTA_o material.l_HATCN material.l_Ag_12nm_test5 material.l_ZnS 1.65*ones(wavelength_num,1)];
% ne_bar=[ones(wavelength_num,1) material.l_Al  material.l_HATCN material.l_TCTA_e material.l_HATCN material.l_Ag_12nm_test5 material.l_ZnS 1.65*ones(wavelength_num,1)];

% no_bar=[ones(wavelength_num,1) material.Ag material.ITO material.l_HATCN(1:301) material.l_NPB(1:301)  material.l_TCTA_o(1:301) material.CBP material.TPBi material.Ag_12nm_test5 material.l_NPB(1:301) material.l_LiF(1:301) ones(301,1)];
% ne_bar=[ones(wavelength_num,1) material.Ag material.ITO material.l_HATCN(1:301) material.l_NPB(1:301)  material.l_TCTA_e(1:301) material.CBP material.TPBi material.Ag_12nm_test5 material.l_NPB(1:301) material.l_LiF(1:301) ones(301,1)];

%% JO ���� %%
no_bar=[ones(wavelength_num,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp  1.5*ones(wavelength_num,1)];
ne_bar=[ones(wavelength_num,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp  1.5*ones(wavelength_num,1)];

% no_bar=[ones(wavelength_num,1) material.Al_JO 1.65*ones(wavelength_num,1) 1.599*ones(wavelength_num,1) 1.5*ones(wavelength_num,1) 1.99*ones(wavelength_num,1)  1.5*ones(wavelength_num,1)];
% ne_bar=[ones(wavelength_num,1) material.Al_JO 1.65*ones(wavelength_num,1) 1.599*ones(wavelength_num,1) 1.5*ones(wavelength_num,1) 1.99*ones(wavelength_num,1)  1.5*ones(wavelength_num,1)];
% no_bar=no_bar(1:301,:);ne_bar=ne_bar(1:301,:);

% no_bar=[ones(wavelength_num,1) material.l_Al_JO  1.4*ones(wavelength_num,1) material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_HATCN material.l_Al_JO  material.l_Ag_12nm_test5 material.l_ZnS 1.65*ones(wavelength_num,1)];
% ne_bar=[ones(wavelength_num,1) material.l_Al_JO  1.4*ones(wavelength_num,1) material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_HATCN material.l_Al_JO  material.l_Ag_12nm_test5 material.l_ZnS 1.65*ones(wavelength_num,1)];

% no_bar=[ones(wavelength_num,1) material.l_Al_JO  material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_HATCN material.l_Ag_12nm_test5 material.l_ZnS ones(wavelength_num,1)];
% ne_bar=[ones(wavelength_num,1) material.l_Al_JO  material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_HATCN material.l_Ag_12nm_test5 material.l_ZnS ones(wavelength_num,1)];

% no_bar=[ones(wavelength_num,1) material.l_Al_JO material.l_MoO3_JH material.l_HATCN material.l_NPB material.l_HATCN material.l_NPB material.l_TCTA_o_JO material.l_B3_o_JO material.l_Al_JO material.l_Ag_12nm_test5 material.l_NPB ones(wavelength_num,1)];
% ne_bar=[ones(wavelength_num,1) material.l_Al_JO material.l_MoO3_JH material.l_HATCN material.l_NPB material.l_HATCN material.l_NPB material.l_TCTA_e_JO material.l_B3_e_JO material.l_Al_JO material.l_Ag_12nm_test5 material.l_NPB ones(wavelength_num,1)];

% no_bar=[ones(wavelength_num,1) material.Al 2.4*ones(wavelength_num,1) 1.8*ones(wavelength_num,1) 1.8*ones(wavelength_num,1) material.ITO 1.5*ones(wavelength_num,1)];
% ne_bar=[ones(wavelength_num,1) material.Al 2.4*ones(wavelength_num,1) 1.8*ones(wavelength_num,1) 1.8*ones(wavelength_num,1) material.ITO 1.5*ones(wavelength_num,1)];

% no_bar=[ones(wavelength_num,1) material.l_Al 2*ones(wavelength_num,1) material.l_TCTA_B3_e 1.5*ones(wavelength_num,1) 2*ones(wavelength_num,1) 1.5*ones(wavelength_num,1)];
% ne_bar=[ones(wavelength_num,1) material.l_Al 2*ones(wavelength_num,1) material.l_TCTA_B3_e 1.5*ones(wavelength_num,1) 2*ones(wavelength_num,1) 1.5*ones(wavelength_num,1)];
% ow=ones(wavelength_num,1);

% no_bar=[ow material.l_Al_JO 1.7*ow 1.7*ow 1.7*ow 2*ow 1.5*ow];
% ne_bar=no_bar;
% ITO_temp=complex(1.9*ow,0.003*ow);
% no_bar=[ones(wavelength_num,1) material.Al_liter material.TPYMB_3 2.301*ones(wavelength_num,1) 1.501*ones(wavelength_num,1) ITO_temp 1.501*ones(wavelength_num,1)];
% ne_bar=[ones(wavelength_num,1) material.Al_liter material.TPYMB_3 2.301*ones(wavelength_num,1) 1.501*ones(wavelength_num,1) ITO_temp 1.501*ones(wavelength_num,1)];

% no_bar=[ones(wavelength_num,1) material.l_Al material.l_TCTA_o material.l_lossy material.l_lossless material.l_lossy material.l_B3_o  material.l_NPB_nok material.l_Glass];
% ne_bar=[ones(wavelength_num,1) material.l_Al material.l_TCTA_e material.l_lossy material.l_lossless material.l_lossy material.l_B3_e  material.l_NPB_nok material.l_Glass];

% no_bar=[ones(wavelength_num,1) material.l_Al material.l_IZO material.l_ZnOPEIE material.l_bebq2_JH real(material.l_bebq2_JH)  material.l_TCTA_o material.l_NPB material.l_HATCN material.l_NPB material.l_HATCN material.l_Ag_12nm_test5 material.l_ZnS ones(401,1)];
% ne_bar=[ones(wavelength_num,1) material.l_Al material.l_IZO material.l_ZnOPEIE material.l_bebq2_JH real(material.l_bebq2_JH)  material.l_TCTA_e material.l_NPB material.l_HATCN material.l_NPB material.l_HATCN material.l_Ag_12nm_test5 material.l_ZnS ones(401,1)];

% no_bar=[ones(wavelength_num,1), material.Al 1.7*ones(wavelength_num,1) material.Ag_12nm_test5 material.TiO2_sputter ones(wavelength_num,1)];
% ne_bar=[ones(wavelength_num,1), material.Al 1.7*ones(wavelength_num,1) material.Ag_12nm_test5 material.TiO2_sputter ones(wavelength_num,1)];

% no_bar=[ones(wavelength_num,1), material.l_Al material.l_LiF material.l_Al_cauchy material.l_ZnMgO_HCCho material.l_InP_HCCho_nok material.l_TFB material.l_AI4083 material.l_ITO  material.l_Glass];
% ne_bar=[ones(wavelength_num,1), material.l_Al material.l_LiF material.l_Al_cauchy material.l_ZnMgO_HCCho material.l_InP_HCCho_nok material.l_TFB material.l_AI4083 material.l_ITO  material.l_Glass];

% no_bar=[ones(wavelength_num,1) material.l_Al material.l_bebq2_JH real(material.l_bebq2_JH) material.l_TCTA_o material.l_NPB material.l_HATCN material.l_Ag_12nm_test5 material.l_ZnS ones(wavelength_num,1)];
% ne_bar=[ones(wavelength_num,1) material.l_Al material.l_bebq2_JH real(material.l_bebq2_JH) material.l_TCTA_e material.l_NPB material.l_HATCN material.l_Ag_12nm_test5 material.l_ZnS ones(wavelength_num,1)];
% no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp  1.5*ones(401,1)];
% ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp  1.5*ones(401,1)];
% no_bar=no_bar(1:301,:);
% ne_bar=ne_bar(1:301,:);
% no_bar=[1.5001*ones(wavelength_num,1) 1.55*ones(wavelength_num,1) 1.5001*ones(wavelength_num,1)];
% ne_bar=no_bar;

% no_bar=[ones(401,1) material.l_Al_JO material.l_B4_o material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO 1.4*ones(wavelength_num,1) material.l_ITO_SNU_temp 2.5*ones(wavelength_num,1) 1.5*ones(401,1)];
% ne_bar=[ones(401,1) material.l_Al_JO material.l_B4_e material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO 1.4*ones(wavelength_num,1) material.l_ITO_SNU_temp 2.5*ones(wavelength_num,1) 1.5*ones(401,1)];

% no_bar=[ones(401,1) material.l_Al_JO material.l_B4_o material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp 1.78*ones(401,1)];
% ne_bar=[ones(401,1) material.l_Al_JO material.l_B4_e material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp 1.78*ones(401,1)];

%% POLED
% no_bar=[ones(301,1) material.Ag material.Bphen real(material.l_NPB(1:301,1)) material.l_TAPC(1:301,1) material.l_Ag_12nm_test5(1:301,1) 1.51*ones(301,1)];
% k_SubPC=imag(material.SubPC);
% SubPC=complex(1.7*ones(301,1),k_SubPC*10);
% no_bar=[ones(301,1) material.Ag material.Bphen real(material.l_NPB(1:301,1)) material.l_TAPC(1:301,1) material.SubPC_temp material.l_TAPC(1:301,1) material.l_Ag_12nm_test5(1:301,1) 1.51*ones(301,1)];
% ne_bar=no_bar;
% no_bar=[ones(401,1) material.l_Al_JO material.l_bebq2_temp material.l_bebq2_nok_temp material.l_NPB material.l_MoO3_JH material.l_IZO_SJ_RTP ones(401,1)];
% ne_bar=[ones(401,1) material.l_Al_JO material.l_bebq2_temp material.l_bebq2_nok_temp material.l_NPB material.l_MoO3_JH material.l_IZO_SJ_RTP ones(401,1)];

% no_bar=[ones(301,1) material.Ag_Palik material.BPhen_CS material.BAlq real(material.NPB_MDQ2acac) material.NPB material.Spiro_TTB material.Spiro_Cl6SubPc material.Spiro_TTB material.Ag_Palik 1.51*ones(301,1)];
% no_bar=[ones(301,1) material.Ag_Palik material.BPhen_CS material.BAlq real(material.NPB_MDQ2) material.NPB material.Spiro_TTB material.Ag_Palik 1.51*ones(301,1)];
% ne_bar=no_bar;
d1= 70;  %0:5:200; %Anode TAPC
Nd1=length(d1);
d2= 197;  %0:5:200; %HTL TPBi
Nd2=length(d2);
data_matrix = zeros(Nd1*Nd2,12);
for k1=1:length(d1)
    for k2=1:length(d2)
        %         fprintf('MoO3 = %d, LiF = %d  ', d1(k1), d2(k2))
        %           thickness= [100 5 7 50 10 15 30 12 100 50];
        thickness= [100 50 25 10 35 150];
%         thickness= [100 d1(k1) 40 d2(k2) 15];
        % thickness= 20;
        % dETL=51;
        % dHTL=10;
        % dHigh=16;
        % dITO=65;
        % thickness=[100 dETL 25 10 dHTL 50 dITO dHigh];
        %         thickness= [100 2.5 70 25 10 90 7 2 15 190];
        %           thickness= [100 70 25 10 90 7 12 180];
        %         thickness= [100 d1(k1) 10 200 150];
        %         thickness= [100 10 7 d1(k1) 7 d1(k1) 15 50 1 15 d2(k2)];
        %         thickness =  [100  90 25  30  65  3  12  d1(k1)  d2(k2)  30  500  30  500]; %% thickness of each layer from bottom to top
        %         thickness =  [100  8   60  20  62   3  12  d1(k1)  d2(k2)];
        %         thickness= [100 d1 100 d1 d2 ];
        %         thickness= [100 7 d1(k1) d2(k2) 30 160];
        %         thickness =  [100 10 30 d1(k1) 25 10 d2(k2) 7 d2(k2) 7 12 210 ];
        %         thickness = [100 d1(k1) 30 10 d2(k2) 150];
        %         thickness =  [100  190 100  100 75];
        %         thickness= [100 10 30 30 25 10 100 7 100 7 10 210];

        EML_position=4;

        %         z0=d_slice/2; % emitting position in EML
        z0=12.5;
        u_data_num=997;
        max_u=3;

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

        sumUtot=sum(U_tot,2);

        Purcell_factor=sumUtot./((const+const2)*u_data_num);
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

            Power_ratio_air_matrix(i)=emissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_p(i,1:u_air_max_p(i)))+sum(U_bottom_transmit_thick_s(i,1:u_air_max_s(i))));
            Power_ratio_air2_matrix(i)=emissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,1:u_air_max_p(i)))+sum(U_bottom_transmit_s(i,1:u_air_max_s(i))));
            Power_ratio_sub_matrix(i)=emissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_s(i,1:u_sub_max_s(i))));
            Power_ratio_abs_matrix(i)=emissioneta_sumUtot(i)*(sum(U_bottom_transmit_total_p(i,1:u_sub_max_p(i))-U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_total_s(i,1:u_sub_max_s(i))-U_bottom_transmit_s(i,1:u_sub_max_s(i)))); % �� ������� ������� abs�� bottom ���� abs�� ���
            Power_ratio_wg_matrix(i)=emissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,u_sub_max_p(i)+1:u_data_num))+sum(U_bottom_transmit_s(i,u_sub_max_s(i)+1:u_data_num)));
            Power_ratio_spp_matrix(i)=emissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,max(u_sub_max_p(i),u_data_num)+1:end))+sum(U_bottom_transmit_s(i,max(u_sub_max_s(i),u_data_num)+1:end)));
            EQE_air_matrix(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_p(i,1:u_air_max_p(i)))+sum(U_bottom_transmit_thick_s(i,1:u_air_max_s(i))));
            %%
            EQE_air_matrix_TE(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_s(i,1:u_air_max_s(i))));
            EQE_air_matrix_TMh(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_ph(i,1:u_air_max_p(i))));
            EQE_air_matrix_TMv(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_pv(i,1:u_air_max_p(i))));
            %%
            EQE_air2_matrix(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,1:u_air_max_p(i)))+sum(U_bottom_transmit_s(i,1:u_air_max_s(i))));
            EQE_sub_matrix(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_s(i,1:u_sub_max_s(i))));
            %%
            EQE_sub_matrix_TE(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_s(i,1:u_sub_max_s(i))));
            EQE_sub_matrix_TMh(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_ph(i,1:u_sub_max_p(i))));
            EQE_sub_matrix_TMv(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_pv(i,1:u_sub_max_p(i))));
            %%
            EQE_abs_matrix(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_total_p(i,1:u_sub_max_p(i))-U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_total_s(i,1:u_sub_max_s(i))-U_bottom_transmit_s(i,1:u_sub_max_s(i)))); % �� ������� ������� abs�� bottom ���� abs�� ���
            EQE_wg_matrix(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,u_sub_max_p(i)+1:u_data_num))+sum(U_bottom_transmit_s(i,u_sub_max_s(i)+1:u_data_num)));
            EQE_spp_matrix(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,max(u_sub_max_p(i),u_data_num)+1:end))+sum(U_bottom_transmit_s(i,max(u_sub_max_s(i),u_data_num)+1:end)));

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
        I_sub_sum_30=sum(I_sub_sum(1,1:31).*sin089(1,1:31));
        I_sub_total=I_sub_total/I_sub_total(1);

        EQE_factor_air=pi*sum(I_air_total.*sin089)/90;
        EQE_factor_sub=pi*sum(I_sub_total.*sin089)/90;
        for i=1:wavelength_num
            spec_lambda(i,1)=emission_spectrum(i,1)/((i+wavelength(1)-1)*10^(-9));
        end
        % 400:700
%         CE= 683*6.626*10^(-34)*(3*10^8)/(1.6*10^-(19))*sum(V_301.*spec_lambda.*eta_eff.*P_air(:,1)./Purcell_factor);
        %         400:800
        %         CE= 683*6.626*10^(-34)*(3*10^8)/(1.6*10^-(19))*sum(V_401.*spec_lambda.*eta_eff.*P_air(:,1)./Purcell_factor);
        %         CE_sub= 683*6.626*10^(-34)*(3*10^8)/(1.6*10^-(19))*sum(V_401.*spec_lambda.*eta_eff.*P_sub(:,1)./Purcell_factor);


        %   I_air_PD=P_air.*repmat(PD_responsivity350850(101:451).*emission_spectrum.*eta_eff./(Purcell_factor),1,90);

        %   I_air_PD_total=sum(I_air_PD);
        %  I_air_PD_total=I_air_PD_total/I_air_PD_total(1);
        %
        %  I_sub_PD=P_sub.*repmat(PD_responsivity350850(101:451).*emission_spectrum.*eta_eff./(Purcell_factor),1,90);

        %  I_sub_PD_total=sum(I_sub_PD);
        %  I_sub_PD_total=I_sub_PD_total/I_sub_PD_total(1);

        I_FWHM=I_air(:,1);
        I_FWHM=I_FWHM/max(I_FWHM); %normalized ���� spectrum
        FWHM=sum(I_FWHM>=0.5); % ���� spectrum�� ��ġ�� ���ϱ�

        fprintf('\nEQE_air = %d \n', EQE_air);
        fprintf('EQE_subconfined = %d \n', EQE_sub_confined);
        fprintf('EQE_wg = %d \n', EQE_wg);
        fprintf('EQE_spp = %d \n', EQE_spp);
        fprintf('EQE_abs = %d \n', EQE_abs);
        fprintf('EQE_sub = %d \n', EQE_sub);
        fprintf('FWHM = %d \n', FWHM);
        disp('-------------------------------')

        k_all=length(d2)*(k1-1)+k2;

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
% Made by JH Kim
% Copyright ??All Rights Reserved.
% clear;
clc;
tic;
clear;
load('nk_JH_total.mat')
% load('Reflectance_data_JH_full_lambda_JO.mat')
% load('BSDF_250103.mat')
% BSDF_MLA=BSDF_matrix;
load('Photopic_400_800.mat')
% BSDF_flat=flat_sub_fresnel(1.5,1);
% BBSDF=BSDF_flat.BSDF;
% BSDF=BBSDF;
% load('n180sub_Lighttools.mat')
% load('BSDF_IMLA_n150.mat')
load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat')
BSDF=BSDF_MLA(:,:,8);
% load('BSDF_SL_S_1to15by2_G_0.5to1by0.1_nSL_150.mat')
% load('BSDF_SL_S_3.35_G_0.9_nSL_150_JO.mat')
% load('hexagonal_half_cone_MLA_BSDF(n_MLA=1.65,scale=0.8-0.05-1,height=1-0.2-2.0).mat')
wavelength=(400:800)';
wavelength_num=length(wavelength);
sin089=sind(0:89);
cos089=cosd(0:89);
emission_spectrum=spectrum.l_I_Irmphmq2tmd_measure_JH; %eta_rad=0.96, hdr=0.82
% emission_spectrum=spectrum.l_I_IrMDQ2acac; %eta_rad=0.82, hdr=0.76
% emission_spectrum=spectrum.I_Irdmppyphtmd; %% JOSong, eta_rad=0.98, hdr=0.865
% emission_spectrum=spectrum.l_I_Irdmppyph2tmd;
% emission_spectrum=zeros(401,1);
% emission_spectrum(181:183,1)=1;
% emission_spectrum=spectrum.Irppy2acac_ETRI;
emission_spectrum=emission_spectrum(1:401);
emission_spectrum=emission_spectrum/sum(emission_spectrum);

eta_rad=0.96;
horizontal_dipole_ratio=0.82;
bottom_air_refractive_index=ones(wavelength_num,1);

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

%% JO 援ъ“ %%
% no_bar=[material.Air material.Al material.B3PyMPM_o material.TCTA_B3PyMPM_1_1_o material.TCTA_o material.TAPC_o material.ITO  material.n165];
% ne_bar=[material.Air material.Al material.B3PyMPM_e material.TCTA_B3PyMPM_1_1_e material.TCTA_e material.TAPC_e material.ITO  material.n165];
%%

% no_bar=[ones(wavelength_num,1) material.l_Al material.l_TCTA_o material.l_lossy material.l_lossless material.l_lossy material.l_B3_o  material.l_NPB_nok material.l_Glass];
% ne_bar=[ones(wavelength_num,1) material.l_Al material.l_TCTA_e material.l_lossy material.l_lossless material.l_lossy material.l_B3_e  material.l_NPB_nok material.l_Glass];

d1=45;  %0:5:200; %Anode TAPC
Nd1=length(d1);

d2=80;  %0:5:200; %HTL TPBi
Nd2=length(d2);
% d_slice=5;

R_index=1;
% n_MLA_list=1.3:0.05:2;
n_MLA_list=1.75;
n_MLA_index=0;
data_matrix= zeros(Nd1*Nd2*size(n_MLA_index,2),13);
Power_matrix=zeros(Nd1*Nd2*size(n_MLA_index,2),13);
for i=1:wavelength_num
    spec_lambda(i,1)=emission_spectrum(i,1)/((i+wavelength(1)-1)*10^(-9));
end
V_301=V_401(1:301,1);
for n_MLA=n_MLA_list
    %     no_bar=[ones(wavelength_num,1) material.l_Al material.l_IZO material.l_ZnOPEIE material.l_B3_o material.l_TCTA_B3_o     material.l_TCTA_o material.l_NPB material.l_HATCN material.l_NPB material.l_HATCN material.l_Ag_12nm_test material.l_NPB n_MLA*ones(401,1)];
    %     ne_bar=[ones(wavelength_num,1) material.l_Al material.l_IZO material.l_ZnOPEIE material.l_B3_e material.l_TCTA_B3_e_nok material.l_TCTA_e material.l_NPB material.l_HATCN material.l_NPB material.l_HATCN material.l_Ag_12nm_test material.l_NPB n_MLA*ones(401,1)];

    %     no_bar=[ones(wavelength_num,1) material.l_Al material.l_IZO material.l_ZnOPEIE material.l_B4_o material.l_bebq2_nok_temp material.l_TCTA_o material.l_NPB material.l_HATCN material.l_NPB material.l_HATCN material.l_Ag_12nm_test material.l_ZnS n_MLA*ones(401,1)];
    %     ne_bar=[ones(wavelength_num,1) material.l_Al material.l_IZO material.l_ZnOPEIE material.l_B4_e material.l_bebq2_nok_temp material.l_TCTA_e material.l_NPB material.l_HATCN material.l_NPB material.l_HATCN material.l_Ag_12nm_test material.l_ZnS n_MLA*ones(401,1)];

    %     no_bar=[ones(wavelength_num,1) material.l_Al material.l_IZO material.l_ZnOPEIE material.l_bebq2_JH real(material.l_bebq2_JH)  material.l_TCTA_o material.l_NPB material.l_HATCN material.l_NPB material.l_HATCN material.l_Ag_12nm_test5 material.l_ZnS material.l_Al2O3 material.l_ParyC material.l_Al2O3 material.l_ParyC n_MLA*ones(401,1)];
    %     ne_bar=[ones(wavelength_num,1) material.l_Al material.l_IZO material.l_ZnOPEIE material.l_bebq2_JH real(material.l_bebq2_JH)  material.l_TCTA_e material.l_NPB material.l_HATCN material.l_NPB material.l_HATCN material.l_Ag_12nm_test5 material.l_ZnS material.l_Al2O3 material.l_ParyC material.l_Al2O3 material.l_ParyC n_MLA*ones(401,1)];

    %     no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp  1.75*ones(401,1)];
    %     ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp  1.75*ones(401,1)];

    % no_bar=[ones(401,1) material.l_Al_JO material.l_TAPC_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO material.l_B3_o_JO material.l_Ag_12nm_test5 1.9*ones(401,1) 1.5*ones(401,1)];
    % ne_bar=[ones(401,1) material.l_Al_JO material.l_TAPC_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO material.l_B3_e_JO material.l_Ag_12nm_test5 1.9*ones(401,1) 1.5*ones(401,1)];
    %% DBR
    % point=[63.5986485202018	84.0002184779090	107.024852922091	111.467523502461	111.375476566327	15	41.9196563306408	2.22273344914414]; % Ag-Ag temp
    % point=[11.4389226329117	125.619774999932	58.6669246832745	148.539765657879	0.222731292461696	15.0001565171575	86.7093383701562	2.40000000000000]; % ITO-Ag temp
    % point=[12	126	59	149	20	15	87	2.4]; % ITO-Ag temp
    % point=[98.5896045828191	136.898008571782	71.0606790753657	119.712406005383	5	86.0300814908642	2.39797636156170]; % ITO-Ag temp
    point=[10	1.593865561014864e+02	1.236739749536133e+02	62.526201369314421];
    high=real(material.l_TiO2_SJ_RTP);
    low=real(material.l_SiO2_SJ_RTP);
    dETL=point(1);
    dHTL=point(2);
    dHigh=point(3);
    dLow=point(4);
    % dIZO=point(5);
    % dAg=point(5);
    % dCap=point(6);
    % nCap=point(7);
    % no_bar=[ones(401,1) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) material.l_Ag_12nm_test5 material.l_B3_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO 0.8*material.l_TAPC_o_JO material.l_Ag_12nm_test5 1.9*ones(401,1) 1.65*ones(401,1)];
    % ne_bar=[ones(401,1) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) material.l_Ag_12nm_test5 material.l_B3_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO 0.8*material.l_TAPC_e_JO material.l_Ag_12nm_test5 1.9*ones(401,1) 1.65*ones(401,1)];
    % thickness=[dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dAg1 80 25 10 55 dAg2 dCap];
    % EML_position=25;
    % no_bar=[ones(401,1) high low high low high low high low high low high low high low material.l_ITO_SNU_temp material.l_TAPC_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO material.l_B3_o_JO material.l_ITO_SNU_temp 2.5*ones(401,1)];
    % ne_bar=[ones(401,1) high low high low high low high low high low high low high low material.l_ITO_SNU_temp material.l_TAPC_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO material.l_B3_o_JO material.l_ITO_SNU_temp 2.5*ones(401,1)];
    % no_bar=[ones(401,1) material.l_Ag_McPeak material.l_TAPC_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO material.l_B3_o_JO material.l_ITO_SNU_temp 1.77*ones(401,1)];
    % ne_bar=[ones(401,1) material.l_Ag_McPeak material.l_TAPC_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO material.l_B3_o_JO material.l_ITO_SNU_temp 1.77*ones(401,1)];

    % thickness=[1000 dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow 50 dHTL 10 25 dETL 50];
        % thickness=[1000 100 dHTL 10 25 dETL 50];


    % no_bar=[ones(401,1) material.l_Si_liter material.l_Al_JO real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) material.l_TAPC_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO material.l_B3_o_JO material.l_Ag_12nm_test5 nCap*ones(401,1) 1.65*ones(401,1)];
    % ne_bar=[ones(401,1) material.l_Si_liter material.l_Al_JO real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) material.l_TAPC_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO material.l_B3_e_JO material.l_Ag_12nm_test5 nCap*ones(401,1) 1.65*ones(401,1)];
    % no_bar=no_bar(1:5:401,:);
    % ne_bar=ne_bar(1:5:401,:);
    % no_bar_OLED=[ones(401,1) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO material.l_TAPC_o_JO material.l_Ag_12nm_test5 nCap*ones(401,1) 1.65*ones(401,1)];
    % ne_bar_OLED=[ones(401,1) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO material.l_TAPC_e_JO material.l_Ag_12nm_test5 nCap*ones(401,1) 1.65*ones(401,1)];
    % no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO  material.l_TCTA_B3_o_JO material.l_TCTA_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
    % ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO  material.l_TCTA_B3_e_JO material.l_TCTA_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];

    % layer_num=size(no_bar,2);

    % *** dETL과 dHTL이 포함된 thickness_OLED 벡터 ***
    % thickness_OLED=[dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dIZO dETL 25 10 dHTL dAg dCap];
    % thickness=[100 100 dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dHTL 10 25 dETL dAg dCap];

    % no_bar=[ones(401,1) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) material.l_ITO_SNU_temp material.l_B3_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO material.l_TAPC_o_JO material.l_Ag_12nm_test5 nCap*ones(401,1) 1.65*ones(401,1)];
    % ne_bar=[ones(401,1) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) material.l_ITO_SNU_temp material.l_B3_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO material.l_TAPC_e_JO material.l_Ag_12nm_test5 nCap*ones(401,1) 1.65*ones(401,1)];

    % no_bar=[ones(401,1) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) material.l_ITO_SNU_temp material.l_B3_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO material.l_TAPC_o_JO material.l_Ag_12nm_test5 nCap*ones(401,1) 1.65*ones(401,1)];
    % ne_bar=[ones(401,1) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) material.l_ITO_SNU_temp material.l_B3_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO material.l_TAPC_e_JO material.l_Ag_12nm_test5 nCap*ones(401,1) 1.65*ones(401,1)];


    % no_bar=[ones(401,1) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) material.l_Ag_12nm_test5 material.l_B3_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO material.l_TAPC_o_JO material.l_Ag_12nm_test5 nCap*ones(401,1) 1.65*ones(401,1)];
    % ne_bar=[ones(401,1) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) real(material.l_LiF) real(material.l_ZnS) material.l_Ag_12nm_test5 material.l_B3_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO material.l_TAPC_e_JO material.l_Ag_12nm_test5 nCap*ones(401,1) 1.65*ones(401,1)];


    % no_bar_OLED=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];
    % ne_bar_OLED=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp 1.51*ones(401,1)];

    % layer_num=size(no_bar,2);

    % *** dETL과 dHTL이 포함된 thickness_OLED 벡터 ***

    % thickness=[dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dIZO dETL 25 10 dHTL dAg dCap];
    % EML_position=5;
    d1=10;
    n_temp=2*ones(401,1);
    k_temp=0.001*ones(401,1);
    ITO_temp=complex(n_temp,k_temp);
    %% Al
    % no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO material.l_TAPC_o_JO ITO_temp 1.5*ones(401,1)];
    % ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO material.l_TAPC_e_JO ITO_temp 1.5*ones(401,1)];

    % no_bar=[ones(401,1) material.l_Al_JO material.l_B3_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO material.l_TAPC_o_JO material.l_ITO_SNU_temp ones(401,1)];
    % ne_bar=[ones(401,1) material.l_Al_JO material.l_B3_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO material.l_TAPC_e_JO material.l_ITO_SNU_temp ones(401,1)];
    % thickness=[100 65 25 10 130 d1(k1)];
    EML_position=4;
    %     no_bar=no_bar(1:301,:);
    %     ne_bar=ne_bar(1:301,:);
        no_bar=[ones(wavelength_num,1) material.l_Al real(material.l_bebq2_JH) real(material.l_bebq2_JH)  material.l_TAPC_o_JO complex(1.6*ones(401,1),0.1*ones(401,1)) complex(2*ones(401,1),0.02*ones(401,1)) ones(401,1)];
        ne_bar=[ones(wavelength_num,1) material.l_Al real(material.l_bebq2_JH) real(material.l_bebq2_JH)  material.l_TAPC_e_JO complex(1.6*ones(401,1),0.1*ones(401,1)) complex(2*ones(401,1),0.02*ones(401,1)) ones(401,1)];
    %     no_bar=[ones(wavelength_num,1) material.l_Al material.l_IZO material.l_ZnOPEIE material.l_B3_o real(material.l_TCTA_B3_o)      material.l_TCTA_o material.l_NPB material.l_HATCN material.l_NPB material.l_HATCN material.l_Ag_12nm_test material.l_ZnS  n_MLA*ones(401,1)];
    %     ne_bar=[ones(wavelength_num,1) material.l_Al material.l_IZO material.l_ZnOPEIE material.l_B3_e real(material.l_TCTA_B3_e_nok)  material.l_TCTA_e material.l_NPB material.l_HATCN material.l_NPB material.l_HATCN material.l_Ag_12nm_test material.l_ZnS  n_MLA*ones(401,1)];

    n_MLA_index=n_MLA_index+1;
    %% single n MLA index
    %     BSDF=normalized_BSDF_matrix(:,:,26);
    %%
    %     BSDF_M=BSDF_MLA(:,:,10);
    %     BSDF=BSDF_M;
    %     BSDF=normalized_BSDF_matrix(:,:,5,5);
    % BSDF=normalized_BSDF_matrix;
    %     BSDF=BSDF_flat;
    %     BSDF=BSDF_MLA;
    BSDF_R(:,:)=BSDF(180:-1:91,:); % ?몃줈異?AOR, 媛?줈異?AOI
    BSDF_T_total=sum(BSDF(1:90,:));

    for k1=1:length(d1)

        for k2=1:length(d2)

            thickness= [100 50 25 65 50 50];
            % thickness=[dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh dLow dHigh d1(k1) dETL 25 10 dHTL 5 dCap];
            z0=12.5;
            u_data_num=998;
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

                    %%%%% 원래 쓰던거 u_sub_max_p(i)=ceil(u_data_num*ne_bar(i,layer_num)/ne_bar(i,EML_position))-1;
                    u_sub_max_p(i)=ceil(u_data_num*ne_bar(i,layer_num)/ne_bar(i,EML_position))-1;

                else

                    %%%%% 원래 쓰던거 u_sub_max_p(i)=ceil(u_data_num*ne_bar(i,layer_num)/ne_bar(i,EML_position));
                    u_sub_max_p(i)=ceil(u_data_num*ne_bar(i,layer_num)/ne_bar(i,EML_position));

                end

                exp_phase=ones(1,u_sub_max_p(i));

                if u_sub_max_p(i)>u_data_num

                    exp_phase(u_data_num+1:u_sub_max_p(i))=exp((-4*pi*no_bar(i,EML_position)*sqrt(u(u_data_num+1:u_sub_max_p(i)).^2-1)*(thickness(EML_position-1)-z0))/wavelength(i));

                end

                K_p_v2(i,1:u_sub_max_p(i))=3/8*ne_bar(i,EML_position)*no_bar(i,layer_num)/no_bar(i,EML_position)^2*sqrt(1-(ne_bar(i,EML_position)*u(1:u_sub_max_p(i))/ne_bar(i,layer_num)).^2).*exp_phase.*u(1:u_sub_max_p(i)).^2.*abs((1+TMF_top.r_p(i,1:u_sub_max_p(i))).*TMF_bottom.t_p(i,1:u_sub_max_p(i))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2./abs(1-u(1:u_sub_max_p(i)).^2);
                K_p_h2(i,1:u_sub_max_p(i))=3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2*(1-(ne_bar(i,EML_position)*u(1:u_sub_max_p(i))/ne_bar(i,layer_num)).^2)).*exp_phase.*abs((1-TMF_top.r_p(i,1:u_sub_max_p(i))).*TMF_bottom.t_p(i,1:u_sub_max_p(i))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2/(12*(no_bar(i,EML_position)/ne_bar(i,EML_position))^2+4);

                if no_bar(i,layer_num)>no_bar(i,EML_position)

                    %%%%%% 원래 쓰던거 % u_sub_max_s(i)=ceil(u_data_num*no_bar(i,layer_num)/no_bar(i,EML_position))-1;
                    u_sub_max_s(i)=ceil(u_data_num*no_bar(i,layer_num)/no_bar(i,EML_position))-1;

                else

                    %%%%%% 원래 쓰던거 u_sub_max_s(i)=ceil(u_data_num*no_bar(i,layer_num)/no_bar(i,EML_position));
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
                Power_ratio_abs_matrix(i)=emissioneta_sumUtot(i)*(sum(U_bottom_transmit_total_p(i,1:u_sub_max_p(i))-U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_total_s(i,1:u_sub_max_s(i))-U_bottom_transmit_s(i,1:u_sub_max_s(i)))); % ??怨꾩궛?쇰줈 ?살뼱吏?뒗 abs??bottom 諛⑺뼢 abs留?怨꾩궛
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
                EQE_abs_matrix(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_total_p(i,1:u_sub_max_p(i))-U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_total_s(i,1:u_sub_max_s(i))-U_bottom_transmit_s(i,1:u_sub_max_s(i)))); % ??怨꾩궛?쇰줈 ?살뼱吏?뒗 abs??bottom 諛⑺뼢 abs留?怨꾩궛
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
                %% P_sub_extrapolation 방지 gemini
                % P_sub_p(i,:)=(ne_bar(i,layer_num)/ne_bar(i,EML_position))*spline(ne_bar(i,EML_position)*u(1:u_sub_max_p(i)),sqrt((ne_bar(i,layer_num)/ne_bar(i,EML_position))^2-u(1:u_sub_max_p(i)).^2).*K_bottom_transmit_p(i,1:u_sub_max_p(i)),ne_bar(i,layer_num)*sin089)/const3(i);
                % P_sub_s(i,:)=(no_bar(i,layer_num)/no_bar(i,EML_position))*spline(no_bar(i,EML_position)*u(1:u_sub_max_s(i)),sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2-u(1:u_sub_max_s(i)).^2).*K_bottom_transmit_s(i,1:u_sub_max_s(i)),no_bar(i,layer_num)*sin089)/const3(i);

                % --- P_sub_p 계산 (외삽 방지) ---

                % 1. 스플라인에 사용할 원본 데이터 (x, y)를 정의합니다.
                u_data_p = u(1:u_sub_max_p(i));
                x_data_p = ne_bar(i,EML_position) * u_data_p;

                % [안전장치] 외삽이 발생하기 직전, 즉 임계각(TIR) 근처에서
                % sqrt() 내부가 음수가 되어 NaN이 발생하는 것을 먼저 방지합니다.
                sqrt_arg_p = (ne_bar(i,layer_num)/ne_bar(i,EML_position))^2 - u_data_p.^2;
                safe_sqrt_arg_p = max(0, sqrt_arg_p); % 음수가 될 경우 0으로 강제

                y_data_p = sqrt(safe_sqrt_arg_p) .* K_bottom_transmit_p(i, 1:u_sub_max_p(i));

                % 2. 스플라인에 요청할 쿼리 지점 (xq) 정의
                xq_data_p = ne_bar(i,layer_num) * sin089;

                % 3. [핵심] 'spline' 대신 'interp1'을 사용합니다.
                %    interp1(x, y, xq, 'method', extrapval)
                interpolated_values_p = interp1(x_data_p, ...  % 원본 x
                    y_data_p, ...  % 원본 y (NaN 방지됨)
                    xq_data_p, ... % 요청 x
                    'spline', ...  % 'spline' 방식을 동일하게 사용
                    0);            % [!!!] 이것이 핵심입니다.
                % x_data_p 범위를 벗어난 모든
                % xq_data_p(예: 89도 값)에 대해
                % 강제로 '0'을 반환시킵니다.

                % 4. 최종 P_sub_p 계산
                P_sub_p(i,:) = (ne_bar(i,layer_num) / ne_bar(i,EML_position)) * interpolated_values_p / const3(i);


                % --- P_sub_s 계산 (동일하게 수정) ---

                u_data_s = u(1:u_sub_max_s(i));
                x_data_s = no_bar(i,EML_position) * u_data_s;

                % [안전장치]
                sqrt_arg_s = (no_bar(i,layer_num)/no_bar(i,EML_position))^2 - u_data_s.^2;
                safe_sqrt_arg_s = max(0, sqrt_arg_s);

                y_data_s = sqrt(safe_sqrt_arg_s) .* K_bottom_transmit_s(i, 1:u_sub_max_s(i));

                xq_data_s = no_bar(i,layer_num) * sin089;

                % [핵심] 외삽 값을 0으로 지정
                interpolated_values_s = interp1(x_data_s, y_data_s, xq_data_s, 'spline', 0);

                P_sub_s(i,:) = (no_bar(i,layer_num) / no_bar(i,EML_position)) * interpolated_values_s / const3(i);
                %%
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
            I_air=P_air.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);

            I_air_total=sum(I_air);
            I_air_total=I_air_total/I_air_total(1);

            I_sub=P_sub.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);

            I_sub_total=sum(I_sub);
            I_sub_total=I_sub_total/I_sub_total(1);

            Psub_norm=I_sub.*sin089;
            temp=repmat(sum(Psub_norm,2),1,90);
            Psub_norm=Psub_norm./temp; clear temp;
            modifiedMatrix = Psub_norm;
            nanIndices = isnan(modifiedMatrix);
            modifiedMatrix(nanIndices) = 0;
            Psub_norm=modifiedMatrix; clear modifiedMatrix nanIndices;
            EQE_factor_air=pi*sum(I_air_total.*sin089)/90;
            EQE_factor_sub=pi*sum(I_sub_total.*sin089)/90;

            I_FWHM=I_air(:,1);
            I_FWHM=I_FWHM/max(I_FWHM); %normalized ?뺣㈃ spectrum
            FWHM=sum(I_FWHM>=0.5); % ?뺣㈃ spectrum??諛섏튂??援ы븯湲?

            fprintf('\nEQE_air = %d \n', EQE_air);
            fprintf('EQE_subconfined = %d \n', EQE_sub_confined);
            fprintf('EQE_wg = %d \n', EQE_wg);
            fprintf('EQE_spp = %d \n', EQE_spp);
            fprintf('EQE_abs = %d \n', EQE_abs);
            fprintf('EQE_sub = %d \n', EQE_sub);
            fprintf('FWHM = %d \n', FWHM);
            disp('-------------------------------')

            k_all=length(d2)*(k1-1)+k2;

            data_matrix(R_index,1:12)=[n_MLA,d1(k1),d2(k2),FWHM,EQE_factor_air,EQE_factor_sub,EQE_air,EQE_sub_confined,EQE_wg,EQE_spp,EQE_abs,EQE_sub];

            TMF_OLED_bottom_p=TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
            TMF_OLED_bottom_s=TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);

            R_p_bottom=abs(TMF_OLED_bottom_p.r_p).^2;
            T_p_bottom=no_bar(:,1)./no_bar(:,layer_num)*(1./cos089).*sqrt(1-(ne_bar(:,layer_num)./ne_bar(:,1)*sin089).^2).*abs(TMF_OLED_bottom_p.t_p).^2;

            R_s_bottom=abs(TMF_OLED_bottom_s.r_s).^2;
            T_s_bottom=no_bar(:,1)./no_bar(:,layer_num)*(1./cos089).*sqrt(1-(no_bar(:,layer_num)./no_bar(:,1)*sin089).^2).*abs(TMF_OLED_bottom_s.t_s).^2;

            for i=1:wavelength_num
                T_p_bottom(i,ceil(asind(ne_bar(i,1)/ne_bar(i,layer_num)))+1:end)=0;
                T_s_bottom(i,ceil(asind(no_bar(i,1)/no_bar(i,layer_num)))+1:end)=0;
            end

            Transmittance=(T_p_bottom+T_s_bottom)/2;
            Reflectance=(R_p_bottom+R_s_bottom)/2;

            BSDF_TT=BSDF';
            BSDF_T=BSDF(1:1:90,:);

            Psub_norm=I_sub.*sin089;
            temp=repmat(sum(Psub_norm,2),1,90);
            Psub_norm=Psub_norm./temp; clear temp;
            modifiedMatrix = Psub_norm;
            nanIndices = isnan(modifiedMatrix);
            modifiedMatrix(nanIndices) = 0;
            Psub_norm=modifiedMatrix; clear modifiedMatrix nanIndices;
            ROLED=Reflectance;
            P0=EQE_sub_matrix;
            %             ROLED=reshape(data_R(R_index,:,1:90),wavelength_num,90,1);
            R_bot=repmat(ROLED,1,1,90);
            Power1=sum(Psub_norm.*repmat(BSDF_T_total,wavelength_num,1),2)'*P0;
            R_matrix_1=repmat(reshape(BSDF_R,1,90,90),wavelength_num,1,1).*R_bot;

            BSDF_CE = zeros(180, 90);
            for i = 1:90
                for j = 1:180
                    BSDF_CE(j,i) = BSDF(j,i)*abs(cosd(j-0.5))/cosd(i-0.5);
                end
            end

            for i=1:wavelength_num
                spec_lambda(i,1)=emission_spectrum(i,1)/((i+wavelength(1)-1)*10^(-9));
            end

            R_step = 30;

            BSDF_R = BSDF(180:-1:91, :);
            BSDF_T_total = sum(BSDF(1:90, :));
            BSDF_T=BSDF(1:1:90,:);
            BSDF_T_CE = BSDF_CE(1:1:90, :);

            P0 = EQE_sub_matrix;
            R_bot = repmat(ROLED, 1, 1, 90);

            for i=1:wavelength_num
                Power1_full_w(i,:)=Psub_norm(i,:)*BSDF_T';
            end

            t_angle = 30;

            Power1_full=sum(Power1_full_w.*repmat(P0,1,90),1);
            Power1_full_t_angle=sum(Power1_full(1:t_angle));

            Power1 = sum(Psub_norm .* repmat(BSDF_T_total, wavelength_num, 1), 2)' * P0;
            R_matrix_1 = repmat(reshape(BSDF_R, 1, 90, 90), wavelength_num, 1, 1) .* R_bot;

            PMLA1 = P_sub*BSDF_T_CE';
            % CE_MLA1= 683*6.626*10^(-34)*(3*10^8)/(1.6*10^-(19))*sum(V_401.*spec_lambda.*eta_eff.*PMLA1(:,1)./Purcell_factor);

            PMLA1 = P_sub*BSDF_T';
            % CE_MLA1= 683*6.626*10^(-34)*(3*10^8)/(1.6*10^-(19))*sum(V_401.*spec_lambda.*eta_eff.*PMLA1(:,1)./Purcell_factor);

            for i = 1:wavelength_num
                temp_psub = Psub_norm(i, :);
                temp_rmat = reshape(R_matrix_1(i, :, :), 90, 90);
                temp_power2(i) = temp_psub * temp_rmat * BSDF_T_total';
                temp_power2_full(:,i)=temp_psub*temp_rmat*BSDF_T';
                temp_power2_forCE(:,i)=temp_psub*temp_rmat*BSDF_T_CE';
            end

            Power2 = temp_power2 * P0;
            Power2_full=temp_power2_full*P0;
            Power2_full_t_angle=sum(Power2_full(1:t_angle));
            PMLA2=temp_power2_forCE'.*P_sub;
            % CE_MLA2= 683*6.626*10^(-34)*(3*10^8)/(1.6*10^-(19))*sum(V_401.*spec_lambda.*eta_eff.*PMLA2(:,1)./Purcell_factor);

            % CE_MLA(1)=CE_MLA1;
            % CE_MLA(2)=CE_MLA2;

            R_matrix = zeros(wavelength_num, 90, 90, 10);
            R_matrix(:, :, :, 1) = repmat(reshape(BSDF_R, 1, 90, 90), wavelength_num, 1, 1) .* R_bot;
            Power = zeros(1, R_step);
            Power(1) = Power1;
            Power(2) = Power2;

            for i = 1:R_step
                for j = 1:wavelength_num
                    R_temp = sum(R_matrix(j, :, :, i), 3);
                    R_temp = reshape(repmat(R_temp, 1, 1, 90), 90, 90)';
                    R_matrix(j, :, :, i + 1) = R_temp .* BSDF_R;
                    R_matrix(j, :, :, i + 1) = R_matrix(j, :, :, i + 1) .* R_bot(j, :, :);
                    temp_power(i + 2, j) = Psub_norm(j, :) * reshape(R_matrix(j, :, :, i + 1), 90, 90) * BSDF_T_total';
                    temp_power_full(i+2,:,j)=Psub_norm(j,:)*reshape(R_matrix(j,:,:,i+1),90,90)*BSDF_T';
                    temp_power_forCE(i+2,:,j)=Psub_norm(j,:)*reshape(R_matrix(j,:,:,i+1),90,90)*BSDF_T_CE';
                end
                Power(i + 2) = temp_power(i + 2, :) * P0;
                % P_MLA(:,:,i+2)=repmat(temp_power_forCE(i+2,:)',1,90).*P_sub;
                % CE_MLA(i+2)= 683*6.626*10^(-34)*(3*10^8)/(1.6*10^-(19))*sum(V_401.*spec_lambda.*eta_eff.*P_MLA(:,1,i+2)./Purcell_factor);
                Power_full(i+2,:)=reshape(temp_power_full(i+2,:,:),90,wavelength_num)*P0;
                Totalpower_t_angle(i+2)=sum(Power_full(i+2,1:t_angle),2);
            end

            Totalpower = sum(Power);

            Totalpower_t_angle(1)=Power1_full_t_angle;
            Totalpower_t_angle(2)=Power2_full_t_angle;
            Power_t_angle=sum(Totalpower_t_angle);

            %             fprintf('EQE_sub = %d, Max EQE = %d, EQE_factor_sub = %d \n',EQE_sub,Totalpower,EQE_factor_sub);
            fprintf('%d of %d \n',R_index,Nd1*Nd2*size(n_MLA_list,2));
            data_matrix(R_index,13)=Totalpower;
            R_index=R_index+1;
        end
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

toc;
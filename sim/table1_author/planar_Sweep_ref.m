% Made by IOEL (J.O.Song W.C.Lee J.H.Park J.H.Kim ...)
% Copyright ⓒ All Rights Reserved.
clear;clc;
tic;

load('nk_JH_total.mat')

wavelength=(400:700).';
wavelength_num=length(wavelength);

emission_spectrum=spectrum.I_Irppy2acac;
emission_spectrum=emission_spectrum/sum(emission_spectrum);

eta_rad=1;
horizontal_dipole_ratio=2/3;
bottom_air_refractive_index=ones(wavelength_num,1);
W=1:301;   % 530-580 nm on the 400-800 nm grid of the material library
no_bar=[ones(wavelength_num,1) material.l_Ag_McPeak(W) material.l_B3_o_JO(W) material.l_TCTA_B3_o_JO(W) material.l_TAPC_o_JO(W) material.l_ITO(W) 1.8*ones(wavelength_num,1)];
ne_bar=[ones(wavelength_num,1) material.l_Ag_McPeak(W) material.l_B3_e_JO(W) material.l_TCTA_B3_e_JO(W) material.l_TAPC_e_JO(W) material.l_ITO(W) 1.8*ones(wavelength_num,1)];

d1= 40; % ex) 0:5:200 
Nd1=length(d1);

d2= 50; % ex) 0:5:200 
Nd2=length(d2);


data_matrix = zeros(Nd1*Nd2,10);

for k1=1:length(d1)
    
    for k2=1:length(d2)
        
        thickness =  [100 300 25 300 50]; %% thickness of each layer from bottom to top
        EML_position=4;
        
        z0=12.5; % emitting position in EML
        
        u_data_num=1000;
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
        U_bottom_transmit_s=2*const4.*K_s_h2;
        U_bottom_transmit=U_bottom_transmit_p+U_bottom_transmit_s;
        
        U_bottom_transmit_total_p=2*(const3.*K_p_v2_total+const4.*K_p_h2_total);
        U_bottom_transmit_total_s=2*const4.*K_s_h2_total;
        U_bottom_transmit_total=U_bottom_transmit_total_p+U_bottom_transmit_total_s;
        
        U_bottom_transmit_thick_p=2*(const3.*K_p_v3+const4.*K_p_h3);
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
            Power_ratio_abs_matrix(i)=emissioneta_sumUtot(i)*(sum(U_bottom_transmit_total_p(i,1:u_sub_max_p(i))-U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_total_s(i,1:u_sub_max_s(i))-U_bottom_transmit_s(i,1:u_sub_max_s(i)))); % 이 계산으로 얻어지는 abs는 bottom 방향 abs만 계산
            Power_ratio_wg_matrix(i)=emissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,u_sub_max_p(i)+1:u_data_num))+sum(U_bottom_transmit_s(i,u_sub_max_s(i)+1:u_data_num)));
            Power_ratio_spp_matrix(i)=emissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,max(u_sub_max_p(i),u_data_num)+1:end))+sum(U_bottom_transmit_s(i,max(u_sub_max_s(i),u_data_num)+1:end)));
            
            EQE_air_matrix(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_thick_p(i,1:u_air_max_p(i)))+sum(U_bottom_transmit_thick_s(i,1:u_air_max_s(i))));
            EQE_air2_matrix(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,1:u_air_max_p(i)))+sum(U_bottom_transmit_s(i,1:u_air_max_s(i))));
            EQE_sub_matrix(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_s(i,1:u_sub_max_s(i))));
            EQE_abs_matrix(i)=lambdaemissioneta_sumUtot(i)*(sum(U_bottom_transmit_total_p(i,1:u_sub_max_p(i))-U_bottom_transmit_p(i,1:u_sub_max_p(i)))+sum(U_bottom_transmit_total_s(i,1:u_sub_max_s(i))-U_bottom_transmit_s(i,1:u_sub_max_s(i)))); % 이 계산으로 얻어지는 abs는 bottom 방향 abs만 계산
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
        
        I_air=P_air.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);
        
        I_air_total=sum(I_air);
        I_air_total=I_air_total/I_air_total(1);
        
        I_sub=P_sub.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);
        
        I_sub_total=sum(I_sub);
        I_sub_total=I_sub_total/I_sub_total(1);
        
        EQE_factor=pi*sum(I_air_total.*sin089)/90;

        I_FWHM=I_air(:,1);
        I_FWHM=I_FWHM/max(I_FWHM); %normalized 정면 spectrum
        FWHM=sum(I_FWHM>=0.5); % 정면 spectrum의 반치폭 구하기
        
        fprintf('EQE_air = %d \n', EQE_air);
        fprintf('FWHM = %d \n', FWHM);
        disp('-------------------------------')
        
        k_all=length(d2)*(k1-1)+k2;
        
        data_matrix(k_all,:)=[d1(k1),d2(k2),0,FWHM,EQE_air,EQE_sub_confined,EQE_wg,EQE_spp,EQE_abs,EQE_sub];
        
    end
end

toc;
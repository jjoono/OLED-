% Table 1 of the manuscript with the author's own model (same machinery as
% Planar_sweep22_preprint_MLA.m): real material stack, hemispherical MLA BSDF
% (slice 11, n_MLA = 1.8), recycling series with R_step terms.
%
% Runs every Table 1 case twice: full wavelength (400-700 nm, Ir(ppy)2acac
% spectrum, photon-number weighting) and single wavelength (550 nm), and writes
% table1_author_model.csv.  Needs nk_JH_total.mat, the BSDF .mat and the
% TMF_birefringence_whole* functions on the path, exactly like the sweep script.
%
% Stack (air side first):
%   air / Ag_McPeak 100 nm / B3PyMPM d_ETL / TCTA:B3PyMPM 25 nm (EML, z0 = 12.5)
%   / TAPC 180 nm / ITO 50 nm (k scaled to 0.0032 or 0.002 at 550 nm) / substrate n = 1.8
clear; clc;
load('nk_JH_total.mat')
load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat')
BSDF = BSDF_MLA(:,:,11);                      % n_MLA = 1.80
BSDF_R = BSDF(180:-1:91,:); BSDF_T_total = sum(BSDF(1:90,:));
R_step = 20; eta_rad = 1;

cases = [];                                   % d_ETL, k_ITO(550), Theta
for d = [200 300]
  for kI = [0.0032 0.002]
    for th = [2/3 0.8 0.9]
      cases(end+1,:) = [d kI th];
    end
  end
end

fid = fopen('table1_author_model.csv','w');
fprintf(fid,'d_ETL_nm,k_ITO_550,Theta,mode,eta_sub,spp,wg,abs,eta_ext,EQE\n');
for mode = 1:2
  if mode == 1, W = 1:301; label = 'full_400_700'; else, W = 151; label = 'single_550'; end
  wavelength = (400:800)'; wavelength = wavelength(W); wavelength_num = length(wavelength);
  emission_spectrum = spectrum.I_Irppy2acac(W); emission_spectrum = emission_spectrum(:)/sum(emission_spectrum);
  bottom_air_refractive_index = ones(wavelength_num,1);
  for c = 1:size(cases,1)
    d_ETL = cases(c,1); k_ITO = cases(c,2); horizontal_dipole_ratio = cases(c,3);
    ito = material.l_ITO(W); ito = ito(:);
    k550 = imag(material.l_ITO(151));           % library k at 550 nm (0.0032285)
    ito = real(ito) + 1i*imag(ito)*(k_ITO/k550); % scale the k curve to the target value at 550 nm
    no_bar = [ones(wavelength_num,1) material.l_Ag_McPeak(W) material.l_B3_o_JO(W) material.l_TCTA_B3_o_JO(W) material.l_TAPC_o_JO(W) ito 1.8*ones(wavelength_num,1)];
    ne_bar = [ones(wavelength_num,1) material.l_Ag_McPeak(W) material.l_B3_e_JO(W) material.l_TCTA_B3_e_JO(W) material.l_TAPC_e_JO(W) ito 1.8*ones(wavelength_num,1)];
    thickness = [100 d_ETL 25 180 50]; EML_position = 4; z0 = 12.5;
    u_data_num = 1000; max_u = 3;
    sin089 = sind(0:89); layer_num = size(no_bar,2);
    u = [(0:u_data_num-1)/u_data_num (u_data_num+1:u_data_num*max_u)/u_data_num]; u_num = length(u);

    TMF_bottom = TMF_birefringence_whole(no_bar(:,EML_position:layer_num),ne_bar(:,EML_position:layer_num),[thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0],u,wavelength);
    TMF_top    = TMF_birefringence_whole(no_bar(:,EML_position:-1:1),ne_bar(:,EML_position:-1:1),[z0 thickness(EML_position-2:-1:1) 0],u,wavelength);
    K_p_v = 3/4*real(ne_bar(:,EML_position)./no_bar(:,EML_position)*(u.^2./sqrt(1-u.^2)).*(1+TMF_bottom.r_p).*(1+TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
    K_p_h = real(3./(6*(no_bar(:,EML_position)./ne_bar(:,EML_position)).^2+2)*sqrt(1-u.^2).*(1-TMF_bottom.r_p).*(1-TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
    K_s_h = real(3./(2*(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2+6)*(1./sqrt(1-u.^2)).*(1+TMF_bottom.r_s).*(1+TMF_top.r_s)./(1-TMF_bottom.r_s.*TMF_top.r_s));
    K_p_v2 = K_p_v; K_p_h2 = K_p_h; K_s_h2 = K_s_h;
    u_sub_max_p = zeros(wavelength_num,1); u_sub_max_s = zeros(wavelength_num,1);
    for i = 1:wavelength_num
      if ne_bar(i,layer_num)>ne_bar(i,EML_position), u_sub_max_p(i)=ceil(u_data_num*ne_bar(i,layer_num)/ne_bar(i,EML_position))-1; else, u_sub_max_p(i)=ceil(u_data_num*ne_bar(i,layer_num)/ne_bar(i,EML_position)); end
      exp_phase = ones(1,u_sub_max_p(i));
      if u_sub_max_p(i)>u_data_num, exp_phase(u_data_num+1:u_sub_max_p(i)) = exp((-4*pi*no_bar(i,EML_position)*sqrt(u(u_data_num+1:u_sub_max_p(i)).^2-1)*(thickness(EML_position-1)-z0))/wavelength(i)); end
      K_p_v2(i,1:u_sub_max_p(i)) = 3/8*ne_bar(i,EML_position)*no_bar(i,layer_num)/no_bar(i,EML_position)^2*sqrt(1-(ne_bar(i,EML_position)*u(1:u_sub_max_p(i))/ne_bar(i,layer_num)).^2).*exp_phase.*u(1:u_sub_max_p(i)).^2.*abs((1+TMF_top.r_p(i,1:u_sub_max_p(i))).*TMF_bottom.t_p(i,1:u_sub_max_p(i))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2./abs(1-u(1:u_sub_max_p(i)).^2);
      K_p_h2(i,1:u_sub_max_p(i)) = 3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2*(1-(ne_bar(i,EML_position)*u(1:u_sub_max_p(i))/ne_bar(i,layer_num)).^2)).*exp_phase.*abs((1-TMF_top.r_p(i,1:u_sub_max_p(i))).*TMF_bottom.t_p(i,1:u_sub_max_p(i))./(1-TMF_bottom.r_p(i,1:u_sub_max_p(i)).*TMF_top.r_p(i,1:u_sub_max_p(i)))).^2/(12*(no_bar(i,EML_position)/ne_bar(i,EML_position))^2+4);
      if no_bar(i,layer_num)>no_bar(i,EML_position), u_sub_max_s(i)=ceil(u_data_num*no_bar(i,layer_num)/no_bar(i,EML_position))-1; else, u_sub_max_s(i)=ceil(u_data_num*no_bar(i,layer_num)/no_bar(i,EML_position)); end
      exp_phase = ones(1,u_sub_max_s(i));
      if u_sub_max_s(i)>u_data_num, exp_phase(u_data_num+1:u_sub_max_s(i)) = exp((-4*pi*no_bar(i,EML_position)*sqrt(u(u_data_num+1:u_sub_max_s(i)).^2-1)*(thickness(EML_position-1)-z0))/wavelength(i)); end
      K_s_h2(i,1:u_sub_max_s(i)) = 3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2-u(1:u_sub_max_s(i)).^2).*exp_phase.*abs((1+TMF_top.r_s(i,1:u_sub_max_s(i))).*TMF_bottom.t_s(i,1:u_sub_max_s(i))./(1-TMF_bottom.r_s(i,1:u_sub_max_s(i)).*TMF_top.r_s(i,1:u_sub_max_s(i)))).^2./((4*(ne_bar(i,EML_position)/no_bar(i,EML_position))^2+12)*abs(1-u(1:u_sub_max_s(i)).^2));
    end
    const = (1-horizontal_dipole_ratio)*ne_bar(:,EML_position)./(wavelength.^4);
    const2 = horizontal_dipole_ratio*no_bar(:,EML_position).*(3+(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2)./(4.*wavelength.^4);
    const3 = const*u; const4 = const2*u;
    U_tot = 2*(const3.*K_p_v+const4.*(K_p_h+K_s_h));
    U_tot_p = 2*(const3.*K_p_v+const4.*K_p_h); U_tot_s = 2*const4.*K_s_h;
    U_bt_p = 2*(const3.*K_p_v2+const4.*K_p_h2); U_bt_s = 2*const4.*K_s_h2;
    sumUtot = sum(U_tot,2);
    ur = 1:u_data_num;
    K_p_v_free = 3/4*real(ne_bar(:,EML_position)./no_bar(:,EML_position)*(u(ur).^2./sqrt(1-u(ur).^2)));
    K_p_h_free = real(3./(6*(no_bar(:,EML_position)./ne_bar(:,EML_position)).^2+2)*sqrt(1-u(ur).^2));
    K_s_h_free = real(3./(2*(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2+6)*(1./sqrt(1-u(ur).^2)));
    P_free = sum(2*((const*u(ur)).*K_p_v_free + (const2*u(ur)).*(K_p_h_free + K_s_h_free)), 2);
    Purcell_factor = sumUtot./P_free; eta_eff = eta_rad*Purcell_factor./(1-eta_rad+eta_rad*Purcell_factor);
    E_sub = zeros(wavelength_num,1); E_wg = E_sub; E_spp = E_sub; E_abs = E_sub;
    for i = 1:wavelength_num
      ip_sub = 1:u_sub_max_p(i); is_sub = 1:u_sub_max_s(i);
      ip_wg = u_sub_max_p(i)+1:u_data_num; is_wg = u_sub_max_s(i)+1:u_data_num;
      ip_spp = max(u_sub_max_p(i),u_data_num)+1:u_num; is_spp = max(u_sub_max_s(i),u_data_num)+1:u_num;
      P_sub = sum(U_bt_p(i,ip_sub))+sum(U_bt_s(i,is_sub)); P_totsub = sum(U_tot_p(i,ip_sub))+sum(U_tot_s(i,is_sub));
      E_sub(i) = P_sub; E_abs(i) = P_totsub-P_sub;
      E_wg(i) = sum(U_tot_p(i,ip_wg))+sum(U_tot_s(i,is_wg)); E_spp(i) = sum(U_tot_p(i,ip_spp))+sum(U_tot_s(i,is_spp));
    end
    photon_weight = wavelength.*emission_spectrum.*eta_eff./sumUtot./sum(wavelength.*emission_spectrum);
    EQE_sub_matrix = photon_weight.*E_sub;
    eta_sub = sum(EQE_sub_matrix); spp = sum(photon_weight.*E_spp); wg = sum(photon_weight.*E_wg); ab = sum(photon_weight.*E_abs);
    % angular distribution in the substrate and the OLED reflectance, as in the sweep script
    K_bt_p = const3.*K_p_v2+const4.*K_p_h2; K_bt_s = const4.*K_s_h2; const_free = pi*(const+const2);
    P_sub_ang = zeros(wavelength_num,90);
    for i = 1:wavelength_num
      rp = ne_bar(i,layer_num)/ne_bar(i,EML_position); up = u(1:u_sub_max_p(i));
      P_sub_ang(i,:) = rp*interp1(ne_bar(i,EML_position)*up, sqrt(max(0,rp^2-up.^2)).*K_bt_p(i,1:u_sub_max_p(i)), ne_bar(i,layer_num)*sin089,'spline',0)/const_free(i);
      rs = no_bar(i,layer_num)/no_bar(i,EML_position); us = u(1:u_sub_max_s(i));
      P_sub_ang(i,:) = P_sub_ang(i,:) + rs*interp1(no_bar(i,EML_position)*us, sqrt(max(0,rs^2-us.^2)).*K_bt_s(i,1:u_sub_max_s(i)), no_bar(i,layer_num)*sin089,'spline',0)/const_free(i);
    end
    I_sub = P_sub_ang.*repmat(emission_spectrum.*eta_eff./Purcell_factor,1,90);
    % P_sub_ang is already the power per unit polar angle (it goes as sin(theta) for a
    % homogeneous medium, see test_psub_shape.m), so it is NOT multiplied by sin089 again
    % as in the sweep script; the BSDF columns take the power per 1-degree bin.
    Psub_norm = I_sub; Psub_norm = Psub_norm./repmat(sum(Psub_norm,2),1,90); Psub_norm(isnan(Psub_norm)) = 0;
    Tp = TMF_birefringence_whole_p(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],ne_bar(:,layer_num)*sin089,wavelength);
    Ts = TMF_birefringence_whole_s(no_bar(:,layer_num:-1:1),ne_bar(:,layer_num:-1:1),[0 thickness(layer_num-2:-1:1) 0],no_bar(:,layer_num)*sin089,wavelength);
    ROLED = (abs(Tp.r_p).^2+abs(Ts.r_s).^2)/2;
    Power = zeros(wavelength_num,R_step);
    for i = 1:wavelength_num
      M = (BSDF_R').*repmat(ROLED(i,:),90,1);   % M(a,b) = BSDF_R(b,a) * R_LED(b)  (row vector convention v = v*M)
      v = Psub_norm(i,:);
      for j = 1:R_step, Power(i,j) = v*BSDF_T_total'; v = v*M; end
    end
    EQE_MLA = sum(sum(Power,2).*EQE_sub_matrix); eta_ext = EQE_MLA/eta_sub;
    fprintf('%s d_ETL %3d k %.4f Theta %.3f | eta_sub %.4f spp %.4f wg %.4f abs %.4f | eta_ext %.4f EQE %.4f\n', label, d_ETL, k_ITO, horizontal_dipole_ratio, eta_sub, spp, wg, ab, eta_ext, EQE_MLA);
    fprintf(fid,'%d,%.4f,%.3f,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n', d_ETL, k_ITO, horizontal_dipole_ratio, label, eta_sub, spp, wg, ab, eta_ext, EQE_MLA);
  end
end
fclose(fid); disp('table1_author_model.csv written');

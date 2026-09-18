% Power-dissipation spectrum U_tot(u, d_ETL) for the generic Ag stack.
% Writes pds_map.csv: row 1 = u grid, column 1 = d_ETL (nm), rest = U_tot (a.u.).
load('nk_JH_total.mat')
wavelength=550; wavelength_num=1;
horizontal_dipole_ratio=2/3;
no_bar=[1 material.l_Ag_McPeak(151) 1.8 1.8 1.8 1.8+0.02i 1.8];
ne_bar=no_bar;                       % isotropic ETL (n_e = n_o = 1.8)
layer_num=size(no_bar,2); EML_position=4; z0=10;
u_data_num=997; max_u=1.5;
u=[(0:u_data_num-1)/u_data_num (u_data_num+1:u_data_num*max_u)/u_data_num];
u_num=length(u);
d1=10:5:500;                          % d_ETL
M=zeros(length(d1),u_num);
for k=1:length(d1)
  thickness=[100 d1(k) 20 50 50];     % Ag, ETL, EML, HTL, TCO
  TMF_bottom=TMF_birefringence_whole(no_bar(:,EML_position:layer_num),ne_bar(:,EML_position:layer_num),[thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0],u,wavelength);
  TMF_top   =TMF_birefringence_whole(no_bar(:,EML_position:-1:1),ne_bar(:,EML_position:-1:1),[z0 thickness(EML_position-2:-1:1) 0],u,wavelength);
  K_p_v=3/4*real(ne_bar(:,EML_position)./no_bar(:,EML_position)*(u.^2./sqrt(1-u.^2)).*(1+TMF_bottom.r_p).*(1+TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
  K_p_h=real(3./(6*(no_bar(:,EML_position)./ne_bar(:,EML_position)).^2+2)*sqrt(1-u.^2).*(1-TMF_bottom.r_p).*(1-TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
  K_s_h=real(3./(2*(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2+6)*(1./sqrt(1-u.^2)).*(1+TMF_bottom.r_s).*(1+TMF_top.r_s)./(1-TMF_bottom.r_s.*TMF_top.r_s));
  const =(1-horizontal_dipole_ratio)*ne_bar(:,EML_position)./(wavelength.^4);
  const2=horizontal_dipole_ratio*no_bar(:,EML_position).*(3+(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2)./(4.*wavelength.^4);
  M(k,:)=2*((const*u).*K_p_v+(const2*u).*(K_p_h+K_s_h));
end
csvwrite('pds_map.csv',[0 u; d1(:) M]);
printf('done: %d x %d, u max %.3f, range %.3g .. %.3g\n', size(M,1), size(M,2), max(u), min(M(:)), max(M(:)));

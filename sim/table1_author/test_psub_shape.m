% Shape test of the substrate angular distribution used to feed the BSDF series:
% homogeneous medium (every layer n = 1.8, no absorption, 'substrate' n = 1.8).
% An isotropic dipole then radiates uniformly per solid angle, so the power per
% unit polar angle must go as sin(theta).  Which of P_sub_ang and P_sub_ang.*sin
% follows sin(theta) tells whether the extra sin089 factor in Psub_norm is right.
wavelength = 550; wavelength_num = 1; horizontal_dipole_ratio = 2/3;
no_bar = 1.8*ones(1,7); ne_bar = no_bar; thickness = [100 200 25 180 50]; EML_position = 4; z0 = 12.5;
u_data_num = 1000; max_u = 3; sin089 = sind(0:89); layer_num = 7;
u = [(0:u_data_num-1)/u_data_num (u_data_num+1:u_data_num*max_u)/u_data_num];
TMF_bottom = TMF_birefringence_whole(no_bar(:,EML_position:layer_num),ne_bar(:,EML_position:layer_num),[thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0],u,wavelength);
TMF_top    = TMF_birefringence_whole(no_bar(:,EML_position:-1:1),ne_bar(:,EML_position:-1:1),[z0 thickness(EML_position-2:-1:1) 0],u,wavelength);
i = 1; usm = u_data_num;   % n_sub = n_EML: substrate range is u < 1
K_p_v2 = 3/8*ne_bar(i,EML_position)*no_bar(i,layer_num)/no_bar(i,EML_position)^2*sqrt(1-(ne_bar(i,EML_position)*u(1:usm)/ne_bar(i,layer_num)).^2).*u(1:usm).^2.*abs((1+TMF_top.r_p(i,1:usm)).*TMF_bottom.t_p(i,1:usm)./(1-TMF_bottom.r_p(i,1:usm).*TMF_top.r_p(i,1:usm))).^2./abs(1-u(1:usm).^2);
K_p_h2 = 3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2*(1-(ne_bar(i,EML_position)*u(1:usm)/ne_bar(i,layer_num)).^2)).*abs((1-TMF_top.r_p(i,1:usm)).*TMF_bottom.t_p(i,1:usm)./(1-TMF_bottom.r_p(i,1:usm).*TMF_top.r_p(i,1:usm))).^2/(12*(no_bar(i,EML_position)/ne_bar(i,EML_position))^2+4);
K_s_h2 = 3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2-u(1:usm).^2).*abs((1+TMF_top.r_s(i,1:usm)).*TMF_bottom.t_s(i,1:usm)./(1-TMF_bottom.r_s(i,1:usm).*TMF_top.r_s(i,1:usm))).^2./((4*(ne_bar(i,EML_position)/no_bar(i,EML_position))^2+12)*abs(1-u(1:usm).^2));
const = (1-horizontal_dipole_ratio)*ne_bar(i,EML_position)/wavelength^4; const2 = horizontal_dipole_ratio*no_bar(i,EML_position)*(3+(ne_bar(i,EML_position)/no_bar(i,EML_position))^2)/(4*wavelength^4);
K_bt_p = const*u(1:usm).*K_p_v2 + const2*u(1:usm).*K_p_h2; K_bt_s = const2*u(1:usm).*K_s_h2; const_free = pi*(const+const2);
rp = 1; up = u(1:usm);
P_sub_ang = rp*interp1(ne_bar(i,EML_position)*up, sqrt(max(0,rp^2-up.^2)).*K_bt_p, ne_bar(i,layer_num)*sin089,'spline',0)/const_free ...
          + rp*interp1(no_bar(i,EML_position)*up, sqrt(max(0,rp^2-up.^2)).*K_bt_s, no_bar(i,layer_num)*sin089,'spline',0)/const_free;
th = 0:89; ref = sind(th);
a = P_sub_ang/max(P_sub_ang); b = (P_sub_ang.*sin089)/max(P_sub_ang.*sin089); r = ref/max(ref);
fprintf('theta   P_sub_ang(norm)   P_sub_ang*sin(norm)   sin(theta)(norm)\n');
for t = [11 21 31 46 61 76 86], fprintf('%4d      %.3f             %.3f              %.3f\n', th(t), a(t), b(t), r(t)); end
fprintf('rms deviation from sin(theta):  P_sub_ang %.4f   P_sub_ang*sin %.4f\n', sqrt(mean((a-r).^2)), sqrt(mean((b-r).^2)));

wavelength=550; u_data_num=997; max_u=3;
u=[(0:u_data_num-1)/u_data_num (u_data_num+1:u_data_num*max_u)/u_data_num];
n_Al=complex(0.9656319200398573,6.458058124032318);
no_bar=[1 n_Al 1.75 1.9 1.5]; ne_bar=no_bar;
thickness=[200 100 100]; EML_position=3; z0=50; layer_num=5;
TMF_bottom=TMF_birefringence_whole(no_bar(:,EML_position:layer_num),ne_bar(:,EML_position:layer_num),[thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0],u,wavelength);
TMF_top   =TMF_birefringence_whole(no_bar(:,EML_position:-1:1),     ne_bar(:,EML_position:-1:1),     [z0 thickness(EML_position-2:-1:1) 0],u,wavelength);
rad=u<1; ev=u>1;
for nm={'as-is',''}
  K_p_v=3/4*real((u.^2./sqrt(1-u.^2)).*(1+TMF_bottom.r_p).*(1+TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
  K_s_h=3/8*real((1./sqrt(1-u.^2)).*(1+TMF_bottom.r_s).*(1+TMF_top.r_s)./(1-TMF_bottom.r_s.*TMF_top.r_s));
  K_p_h=3/8*real(sqrt(1-u.^2).*(1-TMF_bottom.r_p).*(1-TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
  fprintf('radiative(u<1): F_v=%.3f  F_sh=%.3f  F_ph=%.3f\n', sum(K_p_v(rad))/u_data_num, sum(K_s_h(rad))/u_data_num, sum(K_p_h(rad))/u_data_num);
  fprintf('evanescent(u>1): F_v=%.3f  F_sh=%.3f  F_ph=%.3f\n', sum(K_p_v(ev))/u_data_num, sum(K_s_h(ev))/u_data_num, sum(K_p_h(ev))/u_data_num);
  break
end
% branch test: conjugate r in evanescent region (equivalent to opposite sqrt branch inside stack)
rp_b=TMF_bottom.r_p; rp_t=TMF_top.r_p; rs_b=TMF_bottom.r_s; rs_t=TMF_top.r_s;
rp_b(ev)=conj(rp_b(ev)); rp_t(ev)=conj(rp_t(ev)); rs_b(ev)=conj(rs_b(ev)); rs_t(ev)=conj(rs_t(ev));
K_p_v=3/4*real((u.^2./sqrt(1-u.^2)).*(1+rp_b).*(1+rp_t)./(1-rp_b.*rp_t));
K_s_h=3/8*real((1./sqrt(1-u.^2)).*(1+rs_b).*(1+rs_t)./(1-rs_b.*rs_t));
K_p_h=3/8*real(sqrt(1-u.^2).*(1-rp_b).*(1-rp_t)./(1-rp_b.*rp_t));
fprintf('conj-branch evanescent: F_v=%.3f  F_sh=%.3f  F_ph=%.3f\n', sum(K_p_v(ev))/u_data_num, sum(K_s_h(ev))/u_data_num, sum(K_p_h(ev))/u_data_num);

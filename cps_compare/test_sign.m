% quick diagnostic: which r_p sign convention makes physical results?
wavelength=550; u_data_num=997; max_u=3;
u=[(0:u_data_num-1)/u_data_num (u_data_num+1:u_data_num*max_u)/u_data_num];
n_Al=complex(0.9656319200398573,6.458058124032318);
no_bar=[1 n_Al 1.75 1.9 1.5]; ne_bar=no_bar;
thickness=[200 100 100]; EML_position=3; z0=50; layer_num=5;

TMF_bottom=TMF_birefringence_whole(no_bar(:,EML_position:layer_num),ne_bar(:,EML_position:layer_num),[thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0],u,wavelength);
TMF_top   =TMF_birefringence_whole(no_bar(:,EML_position:-1:1),     ne_bar(:,EML_position:-1:1),     [z0 thickness(EML_position-2:-1:1) 0],u,wavelength);

fprintf('u=0 : r_p_top=%.4f%+.4fi  r_s_top=%.4f%+.4fi\n', real(TMF_top.r_p(1)),imag(TMF_top.r_p(1)),real(TMF_top.r_s(1)),imag(TMF_top.r_s(1)));
fprintf('u=0 : r_p_bot=%.4f%+.4fi  r_s_bot=%.4f%+.4fi\n', real(TMF_bottom.r_p(1)),imag(TMF_bottom.r_p(1)),real(TMF_bottom.r_s(1)),imag(TMF_bottom.r_s(1)));

for sign_flip=[1 -1]
  rp_b=sign_flip*TMF_bottom.r_p; rp_t=sign_flip*TMF_top.r_p;
  K_p_v=3/4*real((u.^2./sqrt(1-u.^2)).*(1+rp_b).*(1+rp_t)./(1-rp_b.*rp_t));
  K_p_h=real(3/8*sqrt(1-u.^2).*(1-rp_b).*(1-rp_t)./(1-rp_b.*rp_t));
  K_s_h=real(3/8*(1./sqrt(1-u.^2)).*(1+TMF_bottom.r_s).*(1+TMF_top.r_s)./(1-TMF_bottom.r_s.*TMF_top.r_s));
  F_v=sum(K_p_v)/u_data_num; F_h=sum(K_p_h+K_s_h)/u_data_num;
  fprintf('sign=%+d : Purcell_v=%.3f  Purcell_h=%.3f (expect O(1), positive)\n', sign_flip, F_v, F_h);
end

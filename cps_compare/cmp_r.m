wavelength=550; u_test=[0.5 0.95 1.02 1.05 1.5];
n_Al=complex(0.9656319200398573,6.458058124032318);
no_bar=[1 n_Al 1.75 1.9 1.5]; ne_bar=no_bar;
thickness=[200 100 100]; EML_position=3; layer_num=5; z0=50;
TMF_bottom=TMF_birefringence_whole(no_bar(:,EML_position:layer_num),ne_bar(:,EML_position:layer_num),[thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0],u_test,wavelength);
TMF_top   =TMF_birefringence_whole(no_bar(:,EML_position:-1:1),     ne_bar(:,EML_position:-1:1),     [z0 thickness(EML_position-2:-1:1) 0],u_test,wavelength);
for k=1:5
 fprintf('u=%5.2f p: a_b=%.5f%+.5fi  a_t=%.5f%+.5fi\n',u_test(k),real(TMF_bottom.r_p(k)),imag(TMF_bottom.r_p(k)),real(TMF_top.r_p(k)),imag(TMF_top.r_p(k)));
 fprintf('u=%5.2f s: a_b=%.5f%+.5fi  a_t=%.5f%+.5fi\n',u_test(k),real(TMF_bottom.r_s(k)),imag(TMF_bottom.r_s(k)),real(TMF_top.r_s(k)),imag(TMF_top.r_s(k)));
end

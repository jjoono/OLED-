load('nk_JH_total.mat');
wavelength=550; wavelength_num=1; horizontal_dipole_ratio=2/3;
u_data_num=997; max_u=3;
u=[(0:u_data_num-1)/u_data_num (u_data_num+1:u_data_num*max_u)/u_data_num];
function probe(no_bar, ne_bar, thickness, EML_position, z0, u, wavelength, label)
  layer_num=size(no_bar,2); u_data_num=997;
  TMF_bottom=TMF_birefringence_whole(no_bar(:,EML_position:layer_num),ne_bar(:,EML_position:layer_num),[thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0],u,wavelength);
  TMF_top=TMF_birefringence_whole(no_bar(:,EML_position:-1:1),ne_bar(:,EML_position:-1:1),[z0 thickness(EML_position-2:-1:1) 0],u,wavelength);
  K_p_v=3/4*real(ne_bar(:,EML_position)./no_bar(:,EML_position)*(u.^2./sqrt(1-u.^2)).*(1+TMF_bottom.r_p).*(1+TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
  K_p_h=real(3./(6*(no_bar(:,EML_position)./ne_bar(:,EML_position)).^2+2)*sqrt(1-u.^2).*(1-TMF_bottom.r_p).*(1-TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
  K_s_h=real(3./(2*(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2+6)*(1./sqrt(1-u.^2)).*(1+TMF_bottom.r_s).*(1+TMF_top.r_s)./(1-TMF_bottom.r_s.*TMF_top.r_s));
  hdr=2/3; const=(1-hdr)*ne_bar(:,EML_position)./(wavelength.^4); const2=hdr*no_bar(:,EML_position).*(3+(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2)./(4.*wavelength.^4);
  U=2*((const*u).*K_p_v+(const2*u).*(K_p_h+K_s_h)); Uv=2*(const*u).*K_p_v; Uh=2*(const2*u).*(K_p_h+K_s_h);
  S=sum(U);
  seg=@(a,b) sum(U(u>=a & u<b))/S;
  printf('%s\n  total=%.4g  u<0.556(air):%.3f  0.556-0.833(sub):%.3f  0.833-1(wg):%.3f  1-1.11:%.3f  1.11-1.5:%.3f  1.5-3:%.3f\n', label, S, seg(0,0.5556), seg(0.5556,0.8334), seg(0.8334,1), seg(1,1.1112), seg(1.1112,1.5), seg(1.5,3));
  idx=find(u>1.0 & u<1.3); [mx,im]=max(abs(U(idx))); printf('  largest |U| for 1<u<1.3 at u=%.3f: U=%.3g  (Im r_p top=%.3g, Im r_p bottom=%.3g)\n', u(idx(im)), U(idx(im)), imag(TMF_top.r_p(idx(im))), imag(TMF_bottom.r_p(idx(im))));
end
Al=material.l_Al_JO(151);
probe([1 Al 1.8 2+0.05i 1.5],[1 Al 1.8 2+0.05i 1.5],[100 50 50],3,25,u,wavelength,'A) as in preprint: Al(+k) / EML / TCO(2+0.05i) / glass');
probe([1 conj(Al) 1.8 2-0.05i 1.5],[1 conj(Al) 1.8 2-0.05i 1.5],[100 50 50],3,25,u,wavelength,'B) conjugated n,k: Al(-k) / TCO(2-0.05i)');
probe([1 Al 1.8 1.8 1.5],[1 Al 1.8 1.8 1.5],[100 50 50],3,25,u,wavelength,'C) Al(+k), lossless index-matched bottom (no TCO)');
probe([1 conj(Al) 1.8 1.8 1.5],[1 conj(Al) 1.8 1.8 1.5],[100 50 50],3,25,u,wavelength,'D) Al(-k), lossless bottom');
probe([1 Al 1.8 2+0.05i 1.5],[1 Al 1.8 2+0.05i 1.5],[100 200 50],3,100,u,wavelength,'E) as A but EML 200 nm, z0=100');

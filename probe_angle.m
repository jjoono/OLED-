% angle-resolved TOP (Ag) absorption at 550 nm, real Ag stack
clear; load('nk_JH_total.mat')
wl=550; iw=151;
no_bar=[1 material.l_Ag_McPeak(iw) material.l_B3_o_JO(iw) material.l_TCTA_o_JO(iw) material.l_TCTA_B3_o_JO(iw) material.l_TAPC_o_JO(iw) material.l_HATCN(iw) material.l_TAPC_o_JO(iw) material.l_HATCN(iw) material.l_TAPC_o_JO(iw) material.l_HATCN(iw) material.l_ITO_SNU_temp(iw) 1.8];
ne_bar=[1 material.l_Ag_McPeak(iw) material.l_B3_e_JO(iw) material.l_TCTA_e_JO(iw) material.l_TCTA_B3_e_JO(iw) material.l_TAPC_e_JO(iw) material.l_HATCN(iw) material.l_TAPC_e_JO(iw) material.l_HATCN(iw) material.l_TAPC_e_JO(iw) material.l_HATCN(iw) material.l_ITO_SNU_temp(iw) 1.8];
thickness=[100 500 25 10 500 3 500 3 500 3 150];
EML_position=4; z0=12.5; u_data_num=997; max_u=3;
layer_num=size(no_bar,2);
u=[(0:u_data_num-1)/u_data_num (u_data_num+1:u_data_num*max_u)/u_data_num];
u_num=length(u); hdr=0.76;
TMF_bottom=TMF_birefringence_whole(no_bar(:,EML_position:layer_num),ne_bar(:,EML_position:layer_num),[thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0],u,wl);
TMF_top=TMF_birefringence_whole(no_bar(:,EML_position:-1:1),ne_bar(:,EML_position:-1:1),[z0 thickness(EML_position-2:-1:1) 0],u,wl);
rbp=TMF_bottom.r_p; rtp=TMF_top.r_p; rbs=TMF_bottom.r_s; rts=TMF_top.r_s;
K_p_v=3/4*real(ne_bar(EML_position)/no_bar(EML_position)*(u.^2./sqrt(1-u.^2)).*(1+rbp).*(1+rtp)./(1-rbp.*rtp));
K_p_h=real(3/(6*(no_bar(EML_position)/ne_bar(EML_position))^2+2)*sqrt(1-u.^2).*(1-rbp).*(1-rtp)./(1-rbp.*rtp));
K_s_h=real(3/(2*(ne_bar(EML_position)/no_bar(EML_position))^2+6)*(1./sqrt(1-u.^2)).*(1+rbs).*(1+rts)./(1-rbs.*rts));
sq=sqrt(1-u.^2);
K_p_v_b=3/8*u.^2.*real((1+rbp).*(1-conj(rbp))./sq).*abs((1+rtp)./(1-rbp.*rtp)).^2;
K_p_h_b=3*real((1-rbp).*(1+conj(rbp)).*sq).*abs((1-rtp)./(1-rbp.*rtp)).^2/(12*(no_bar(EML_position)/ne_bar(EML_position))^2+4);
K_s_h_b=3*real((1+rbs).*(1-conj(rbs))./sq).*abs((1+rts)./(1-rbs.*rts)).^2/(4*(ne_bar(EML_position)/no_bar(EML_position))^2+12);
c1=(1-hdr)*ne_bar(EML_position); c2=hdr*no_bar(EML_position)*(3+(ne_bar(EML_position)/no_bar(EML_position))^2)/4;
U_tot=2*(c1*u.*K_p_v+c2*u.*(K_p_h+K_s_h));
U_top=U_tot-2*(c1*u.*K_p_v_b+c2*u.*(K_p_h_b+K_s_h_b));
U_top_p=2*(c1*u.*(K_p_v-K_p_v_b)+c2*u.*(K_p_h-K_p_h_b));
sU=sum(U_tot);
theta=asind(min(u,1)); % internal angle in EML
fprintf('top(Ag) share at 550nm = %.4f of dissipated power\n', sum(U_top)/sU);
fprintf('  from u<1 (propagating) = %.4f | u>1 (evanescent) = %.4f\n', sum(U_top(u<1))/sU, sum(U_top(u>=1))/sU);
for band=[0 30; 30 50; 50 70; 70 80; 80 90]'
    idx=u<1 & theta>=band(1) & theta<band(2);
    fprintf('  theta %2d-%2d deg: %.4f (p-pol %.4f)\n', band(1), band(2), sum(U_top(idx))/sU, sum(U_top_p(idx))/sU);
end

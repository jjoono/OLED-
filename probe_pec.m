% PEC-mirror check for the n=1.8, ETL/HTL 500 nm, Ag stack:
% where does every unit of emitted power go?
% Channels: nonradiative (eta_eff), bottom propagating into substrate,
% bottom absorbed (between EML and substrate), bottom evanescent-absorbed,
% and TOP-side dissipation (Ag absorption) which the standard EQE channel
% printout does not include.
clear;
load('nk_JH_total.mat')
wavelength=(400:800)';
wavelength_num=length(wavelength);

emission_spectrum=spectrum.l_Irppy2acac;
emission_spectrum=emission_spectrum(1:401);
emission_spectrum=emission_spectrum/sum(emission_spectrum);

eta_rad=0.98;
hdr=0.76;

no_bar=[ones(401,1) complex(0.001*ones(401,1),12*ones(401,1)) material.l_B3_o_JO material.l_TCTA_o_JO material.l_TCTA_B3_o_JO material.l_TAPC_o_JO material.l_HATCN material.l_TAPC_o_JO material.l_HATCN material.l_TAPC_o_JO material.l_HATCN material.l_ITO_SNU_temp 1.8*ones(401,1)];
ne_bar=[ones(401,1) complex(0.001*ones(401,1),12*ones(401,1)) material.l_B3_e_JO material.l_TCTA_e_JO material.l_TCTA_B3_e_JO material.l_TAPC_e_JO material.l_HATCN material.l_TAPC_e_JO material.l_HATCN material.l_TAPC_e_JO material.l_HATCN material.l_ITO_SNU_temp 1.8*ones(401,1)];
thickness=[100 500 25 10 500 3 500 3 500 3 150];
EML_position=4;
z0=12.5;
u_data_num=997;
max_u=3;
layer_num=size(no_bar,2);
u=[(0:u_data_num-1)/u_data_num (u_data_num+1:u_data_num*max_u)/u_data_num];
u_num=length(u);

TMF_bottom=TMF_birefringence_whole(no_bar(:,EML_position:layer_num),ne_bar(:,EML_position:layer_num),[thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0],u,wavelength);
TMF_top=TMF_birefringence_whole(no_bar(:,EML_position:-1:1),ne_bar(:,EML_position:-1:1),[z0 thickness(EML_position-2:-1:1) 0],u,wavelength);

rbp=TMF_bottom.r_p; rtp=TMF_top.r_p;
rbs=TMF_bottom.r_s; rts=TMF_top.r_s;

% total dissipation (same as main script)
K_p_v=3/4*real(ne_bar(:,EML_position)./no_bar(:,EML_position)*(u.^2./sqrt(1-u.^2)).*(1+rbp).*(1+rtp)./(1-rbp.*rtp));
K_p_h=real(3./(6*(no_bar(:,EML_position)./ne_bar(:,EML_position)).^2+2)*sqrt(1-u.^2).*(1-rbp).*(1-rtp)./(1-rbp.*rtp));
K_s_h=real(3./(2*(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2+6)*(1./sqrt(1-u.^2)).*(1+rbs).*(1+rts)./(1-rbs.*rts));

% net power flux into the BOTTOM half space at the dipole plane, ALL u
sq=sqrt(1-repmat(u,wavelength_num,1).^2);
K_p_v_bot=3/8*repmat(u,wavelength_num,1).^2.*real((1+rbp).*(1-conj(rbp))./sq).*abs((1+rtp)./(1-rbp.*rtp)).^2;
K_p_h_bot=3*real((1-rbp).*(1+conj(rbp)).*sq).*abs((1-rtp)./(1-rbp.*rtp)).^2./repmat(12*(no_bar(:,EML_position)./ne_bar(:,EML_position)).^2+4,1,u_num);
K_s_h_bot=3*real((1+rbs).*(1-conj(rbs))./sq).*abs((1+rts)./(1-rbs.*rts)).^2./repmat(4*(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2+12,1,u_num);

const=(1-hdr)*ne_bar(:,EML_position)./(wavelength.^4);
const2=hdr*no_bar(:,EML_position).*(3+(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2)./(4.*wavelength.^4);

U_tot=2*((const*u).*K_p_v+(const2*u).*(K_p_h+K_s_h));
U_bot=2*((const*u).*K_p_v_bot+(const2*u).*(K_p_h_bot+K_s_h_bot));

sumUtot=sum(U_tot,2);
sumUbot=sum(U_bot,2);

Purcell=sumUtot./((const+const2)*u_data_num);
eta_eff=eta_rad*Purcell./(1-eta_rad+eta_rad*Purcell);

w=wavelength.*emission_spectrum/sum(wavelength.*emission_spectrum); % photon weighting as in EQE

eta_eff_bar=sum(w.*eta_eff);
frac_bot=sum(w.*eta_eff.*sumUbot./sumUtot);
frac_top=eta_eff_bar-frac_bot;

fprintf('Purcell F (emission-weighted)      = %.4f\n', sum(w.*Purcell));
fprintf('eta_eff (emission-weighted)        = %.4f  -> nonradiative loss = %.4f\n', eta_eff_bar, 1-eta_eff_bar);
fprintf('bottom-halfspace share (all u)     = %.4f\n', frac_bot);
fprintf('TOP-halfspace share (Ag side)      = %.4f\n', frac_top);

% angle-resolved top absorption at 550 nm for the REAL Ag stack is printed
% by probe_energy_angle.m

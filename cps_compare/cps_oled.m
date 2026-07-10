%% CPS power-dissipation calculation for the FDTD benchmark OLED stack.
% Adapted from planar_Sweep22.m (WC Lee / JOSong group code): identical
% K-integrand formulas and mode bookkeeping, with the stack, wavelength and
% dipole parameters replaced to match oled_fdtd/oled_meep*.py:
%
%   air | Al 200nm | organic n=1.75 100nm (dipole at centre) | ITO n=1.9
%       100nm | glass n=1.5 (semi-infinite)
%   lambda = 550 nm, eta_rad = 1, isotropic no birefringence (ne = no)
%
% Al refractive index is taken from the SAME material model as the FDTD run
% (Meep materials library, Rakic Drude-Lorentz fit) evaluated at 550 nm, so
% the two methods share every physical input.
%
% Outputs power fractions of the total dissipated power:
%   sub  : entering the glass substrate      (== FDTD frac_glass)
%   wg   : guided, non-evanescent in EML     (~= FDTD frac_wg)
%   abs  : parasitic absorption of radiating modes (Al, here)
%   spp  : evanescent / surface plasmon      (abs+spp ~= FDTD frac_metal)
%   air  : escaping glass->air after ray-optics correction (extra info)
% for pure-horizontal (hdr=1), pure-vertical (hdr=0) and isotropic (2/3).

clc;

wavelength = 550;                 % single wavelength [nm]
wavelength_num = 1;

% ---- Al(550 nm) from Meep's Rakic Drude-Lorentz model (set by run script)
n_Al = complex(0.9656319200398573, 6.458058124032318);

% ---- stack: [air  Al  EML  ITO  glass], EML first column convention below
% A small extinction k is added to the confined layers: with strictly
% lossless indices the ITO/organic guided-mode poles sit ON the real-u axis
% and the discrete CPS u-sum diverges. k = 2e-3 moves them off axis; the
% integrated mode power is insensitive to the exact k once the u grid
% resolves the Lorentzian (hence u_data_num below is 10x the default).
kk = 2e-3;
no_bar = [1  n_Al  1.75+1i*kk  1.9+1i*kk  1.5];
ne_bar = no_bar;                  % no birefringence
thickness = [200 100 100];        % Al, EML, ITO  [nm]
EML_position = 3;
z0 = 50;                          % dipole->Al-side interface distance [nm]
layer_num = size(no_bar, 2);

bottom_air_refractive_index = 1;
eta_rad = 1;

u_data_num = 9970;
max_u = 3;
u = [(0:u_data_num-1)/u_data_num (u_data_num+1:u_data_num*max_u)/u_data_num];
u_num = length(u);

% ---- whole-cavity half-space reflection/transmission (as in planar_Sweep22)
TMF_bottom = TMF_birefringence_whole(no_bar(:,EML_position:layer_num), ne_bar(:,EML_position:layer_num), [thickness(EML_position-1)-z0 thickness(EML_position:layer_num-2) 0], u, wavelength);
TMF_top    = TMF_birefringence_whole(no_bar(:,EML_position:-1:1),      ne_bar(:,EML_position:-1:1),      [z0 thickness(EML_position-2:-1:1) 0],                                u, wavelength);

K_p_v = 3/4*real(ne_bar(:,EML_position)./no_bar(:,EML_position)*(u.^2./sqrt(1-u.^2)).*(1+TMF_bottom.r_p).*(1+TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
K_p_h = real(3./(6*(no_bar(:,EML_position)./ne_bar(:,EML_position)).^2+2)*sqrt(1-u.^2).*(1-TMF_bottom.r_p).*(1-TMF_top.r_p)./(1-TMF_bottom.r_p.*TMF_top.r_p));
K_s_h = real(3./(2*(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2+6)*(1./sqrt(1-u.^2)).*(1+TMF_bottom.r_s).*(1+TMF_top.r_s)./(1-TMF_bottom.r_s.*TMF_top.r_s));

K_p_v2 = K_p_v;  K_p_h2 = K_p_h;  K_s_h2 = K_s_h;
K_p_v2_total = K_p_v;  K_p_h2_total = K_p_h;  K_s_h2_total = K_s_h;
K_p_v3 = K_p_v;  K_p_h3 = K_p_h;  K_s_h3 = K_s_h;

i = 1;   % single wavelength

if real(ne_bar(i,layer_num)) > real(ne_bar(i,EML_position))
    u_sub_max_p = ceil(u_data_num*real(ne_bar(i,layer_num))/real(ne_bar(i,EML_position))) - 1;
else
    u_sub_max_p = ceil(u_data_num*real(ne_bar(i,layer_num))/real(ne_bar(i,EML_position)));
end
exp_phase = ones(1, u_sub_max_p);
if u_sub_max_p > u_data_num
    exp_phase(u_data_num+1:u_sub_max_p) = exp((-4*pi*no_bar(i,EML_position)*sqrt(u(u_data_num+1:u_sub_max_p).^2-1)*(thickness(EML_position-1)-z0))/wavelength(i));
end
K_p_v2(i,1:u_sub_max_p) = 3/8*ne_bar(i,EML_position)*no_bar(i,layer_num)/no_bar(i,EML_position)^2*sqrt(1-(ne_bar(i,EML_position)*u(1:u_sub_max_p)/ne_bar(i,layer_num)).^2).*exp_phase.*u(1:u_sub_max_p).^2.*abs((1+TMF_top.r_p(i,1:u_sub_max_p)).*TMF_bottom.t_p(i,1:u_sub_max_p)./(1-TMF_bottom.r_p(i,1:u_sub_max_p).*TMF_top.r_p(i,1:u_sub_max_p))).^2./abs(1-u(1:u_sub_max_p).^2);
K_p_h2(i,1:u_sub_max_p) = 3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2*(1-(ne_bar(i,EML_position)*u(1:u_sub_max_p)/ne_bar(i,layer_num)).^2)).*exp_phase.*abs((1-TMF_top.r_p(i,1:u_sub_max_p)).*TMF_bottom.t_p(i,1:u_sub_max_p)./(1-TMF_bottom.r_p(i,1:u_sub_max_p).*TMF_top.r_p(i,1:u_sub_max_p))).^2/(12*(no_bar(i,EML_position)/ne_bar(i,EML_position))^2+4);

if real(no_bar(i,layer_num)) > real(no_bar(i,EML_position))
    u_sub_max_s = ceil(u_data_num*real(no_bar(i,layer_num))/real(no_bar(i,EML_position))) - 1;
else
    u_sub_max_s = ceil(u_data_num*real(no_bar(i,layer_num))/real(no_bar(i,EML_position)));
end
exp_phase = ones(1, u_sub_max_s);
if u_sub_max_s > u_data_num
    exp_phase(u_data_num+1:u_sub_max_s) = exp((-4*pi*no_bar(i,EML_position)*sqrt(u(u_data_num+1:u_sub_max_s).^2-1)*(thickness(EML_position-1)-z0))/wavelength(i));
end
K_s_h2(i,1:u_sub_max_s) = 3*sqrt((no_bar(i,layer_num)/no_bar(i,EML_position))^2-u(1:u_sub_max_s).^2).*exp_phase.*abs((1+TMF_top.r_s(i,1:u_sub_max_s)).*TMF_bottom.t_s(i,1:u_sub_max_s)./(1-TMF_bottom.r_s(i,1:u_sub_max_s).*TMF_top.r_s(i,1:u_sub_max_s))).^2./((4*(ne_bar(i,EML_position)/no_bar(i,EML_position))^2+12)*abs(1-u(1:u_sub_max_s).^2));

K_p_v2_total(i,1:u_sub_max_p) = 3/8*(u(1:u_sub_max_p).^2).*real((1+TMF_bottom.r_p(i,1:u_sub_max_p)).*(1-conj(TMF_bottom.r_p(i,1:u_sub_max_p)))./sqrt(1-u(1:u_sub_max_p).^2)).*abs((1+TMF_top.r_p(i,1:u_sub_max_p))./(1-TMF_bottom.r_p(i,1:u_sub_max_p).*TMF_top.r_p(i,1:u_sub_max_p))).^2;
K_p_h2_total(i,1:u_sub_max_p) = 3*real((1-TMF_bottom.r_p(i,1:u_sub_max_p)).*(1+conj(TMF_bottom.r_p(i,1:u_sub_max_p))).*sqrt(1-u(1:u_sub_max_p).^2)).*abs((1-TMF_top.r_p(i,1:u_sub_max_p))./(1-TMF_bottom.r_p(i,1:u_sub_max_p).*TMF_top.r_p(i,1:u_sub_max_p))).^2/(12*(no_bar(i,EML_position)/ne_bar(i,EML_position))^2+4);
K_s_h2_total(i,1:u_sub_max_s) = 3*real((1+TMF_bottom.r_s(i,1:u_sub_max_s)).*(1-conj(TMF_bottom.r_s(i,1:u_sub_max_s)))./sqrt(1-u(1:u_sub_max_s).^2)).*abs((1+TMF_top.r_s(i,1:u_sub_max_s))./(1-TMF_bottom.r_s(i,1:u_sub_max_s).*TMF_top.r_s(i,1:u_sub_max_s))).^2/(4*(ne_bar(i,EML_position)/no_bar(i,EML_position))^2+12);

K_p_v3(i,:) = K_p_v2(i,:);
K_p_h3(i,:) = K_p_h2(i,:);
K_s_h3(i,:) = K_s_h2(i,:);

if bottom_air_refractive_index(i) > real(ne_bar(i,EML_position))
    u_air_max_p = min(u_sub_max_p, ceil(bottom_air_refractive_index(i)*u_data_num/real(ne_bar(i,EML_position))) - 1);
else
    u_air_max_p = min(u_sub_max_p, ceil(bottom_air_refractive_index(i)*u_data_num/real(ne_bar(i,EML_position))));
end
TMF_OLED_bottom_p = TMF_birefringence_whole_p(no_bar(i,layer_num:-1:1), ne_bar(i,layer_num:-1:1), [0 thickness(layer_num-2:-1:1) 0], ne_bar(i,EML_position)*u(1:u_air_max_p), wavelength(i));
R_p_bottom = abs(TMF_OLED_bottom_p.r_p).^2;
cos_theta_sub = sqrt(1-(ne_bar(i,EML_position)*u(1:u_air_max_p)/ne_bar(i,layer_num)).^2);
cos_theta_air = sqrt(1-(ne_bar(i,EML_position)*u(1:u_air_max_p)/bottom_air_refractive_index(i)).^2);
r_p = (bottom_air_refractive_index(i)*cos_theta_sub-no_bar(i,layer_num)*cos_theta_air)./(bottom_air_refractive_index(i)*cos_theta_sub+no_bar(i,layer_num)*cos_theta_air);
R_sub_air_bottom_p = abs(r_p).^2;
T_sub_air_bottom_p = 1 - R_sub_air_bottom_p;

if bottom_air_refractive_index(i) > real(no_bar(i,EML_position))
    u_air_max_s = min(u_sub_max_s, ceil(bottom_air_refractive_index(i)*u_data_num/real(no_bar(i,EML_position))) - 1);
else
    u_air_max_s = min(u_sub_max_s, ceil(bottom_air_refractive_index(i)*u_data_num/real(no_bar(i,EML_position))));
end
TMF_OLED_bottom_s = TMF_birefringence_whole_s(no_bar(i,layer_num:-1:1), ne_bar(i,layer_num:-1:1), [0 thickness(layer_num-2:-1:1) 0], no_bar(i,EML_position)*u(1:u_air_max_s), wavelength(i));
R_s_bottom = abs(TMF_OLED_bottom_s.r_s).^2;
cos_theta_sub = sqrt(1-(no_bar(i,EML_position)*u(1:u_air_max_s)/no_bar(i,layer_num)).^2);
cos_theta_air = sqrt(1-(no_bar(i,EML_position)*u(1:u_air_max_s)/bottom_air_refractive_index(i)).^2);
r_s = (no_bar(i,layer_num)*cos_theta_sub-bottom_air_refractive_index(i)*cos_theta_air)./(no_bar(i,layer_num)*cos_theta_sub+bottom_air_refractive_index(i)*cos_theta_air);
R_sub_air_bottom_s = abs(r_s).^2;
T_sub_air_bottom_s = 1 - R_sub_air_bottom_s;

K_p_v3(i,1:u_air_max_p) = K_p_v2(i,1:u_air_max_p).*T_sub_air_bottom_p./(1-R_p_bottom.*R_sub_air_bottom_p);
K_p_h3(i,1:u_air_max_p) = K_p_h2(i,1:u_air_max_p).*T_sub_air_bottom_p./(1-R_p_bottom.*R_sub_air_bottom_p);
K_s_h3(i,1:u_air_max_s) = K_s_h2(i,1:u_air_max_s).*T_sub_air_bottom_s./(1-R_s_bottom.*R_sub_air_bottom_s);

% ---- assemble per-orientation results for hdr = 1, 0, 2/3 ---------------
for hdr = [1 0 2/3]

    horizontal_dipole_ratio = hdr;

    const  = (1-horizontal_dipole_ratio)*ne_bar(:,EML_position)./(wavelength.^4);
    const2 = horizontal_dipole_ratio*no_bar(:,EML_position).*(3+(ne_bar(:,EML_position)./no_bar(:,EML_position)).^2)./(4.*wavelength.^4);
    const3 = const*u;
    const4 = const2*u;

    U_tot = 2*(const3.*K_p_v+const4.*(K_p_h+K_s_h));

    U_bottom_transmit_p = 2*(const3.*K_p_v2+const4.*K_p_h2);
    U_bottom_transmit_s = 2*const4.*K_s_h2;

    U_bottom_transmit_total_p = 2*(const3.*K_p_v2_total+const4.*K_p_h2_total);
    U_bottom_transmit_total_s = 2*const4.*K_s_h2_total;

    U_bottom_transmit_thick_p = 2*(const3.*K_p_v3+const4.*K_p_h3);
    U_bottom_transmit_thick_s = 2*const4.*K_s_h3;

    sumUtot = sum(U_tot,2);
    Purcell_factor = sumUtot./((const+const2)*u_data_num);
    eta_eff = eta_rad*Purcell_factor./(1-eta_rad+eta_rad*Purcell_factor);
    emissioneta_sumUtot = eta_eff./sumUtot;

    P_air = emissioneta_sumUtot*(sum(U_bottom_transmit_thick_p(1,1:u_air_max_p))+sum(U_bottom_transmit_thick_s(1,1:u_air_max_s)));
    P_sub = emissioneta_sumUtot*(sum(U_bottom_transmit_p(1,1:u_sub_max_p))+sum(U_bottom_transmit_s(1,1:u_sub_max_s)));
    P_abs = emissioneta_sumUtot*(sum(U_bottom_transmit_total_p(1,1:u_sub_max_p)-U_bottom_transmit_p(1,1:u_sub_max_p))+sum(U_bottom_transmit_total_s(1,1:u_sub_max_s)-U_bottom_transmit_s(1,1:u_sub_max_s)));
    P_wg  = emissioneta_sumUtot*(sum(U_bottom_transmit_p(1,u_sub_max_p+1:u_data_num))+sum(U_bottom_transmit_s(1,u_sub_max_s+1:u_data_num)));
    P_spp = emissioneta_sumUtot*(sum(U_bottom_transmit_p(1,max(u_sub_max_p,u_data_num)+1:end))+sum(U_bottom_transmit_s(1,max(u_sub_max_s,u_data_num)+1:end)));

    fprintf('hdr=%.4f  Purcell=%.4f | air %.4f  sub %.4f  wg %.4f  spp %.4f  abs %.4f | sum %.4f\n', ...
        hdr, Purcell_factor, P_air, P_sub, P_wg, P_spp, P_abs, P_sub+P_wg+P_spp+P_abs);
end

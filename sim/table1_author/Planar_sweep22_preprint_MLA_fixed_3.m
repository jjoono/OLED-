%% MLA_pass_resolved.m
% OLED on a substrate with an external outcoupling structure (BSDF): substrate-delivered power,
% multi-pass extraction by the matrix series, and pass-by-pass quantities
%   - eta_sub, channel budget (air / substrate-confined / waveguide / evanescent (SPP) / absorption)
%   - eta_ext (matrix series) and eq. (3) with the Lambertian p and A'
%   - per pass j: power arriving at the structure, escape probability p_j, round-trip loss A'_j,
%     angular intensity in the substrate (arriving light) and in air (extracted light)
%
% Rewritten from Planar_sweep22_preprint_MLA.m (WC Lee). Differences from that script:
%   1. P_sub(theta) is power per unit polar angle (sin theta included); it is used as is
%      (no second sin weight).
%   2. Air-side intensity per solid angle: BSDF(b,a) / sin(theta_b)  (no sin(theta_a) factor).
%   3. R_LED is applied at the angle the light is returned to (after the BSDF reflection).
%   4. One source vector per wavelength, weighted by the substrate-delivered power P0(lambda).
%   5. Angular quantities follow planar_Sweep.m: sin089 grid (0:89 deg), I_sub and I_air are
%      intensities per solid angle weighted by emission_spectrum*eta_eff/Purcell (as in planar_Sweep.m).
%   6. Pass-resolved intensities I_sub_k (substrate, light arriving at pass k) and I_air_k (air, light
%      extracted at pass k) in the same units as I_sub and I_air; I_sub_k(:,:,1) = I_sub exactly.
%
% Needs in the path: nk_JH_total.mat, the BSDF file, TMF_birefringence_whole(_p,_s).m
% (the version with the evanescent-branch fix, imag(n cos) >= 0).
% Runs in MATLAB and Octave (the polar plots are skipped in Octave).

clear; clc;

%% ------------------------------------------------------------------ settings
nkfile   = 'nk_JH_total.mat';
bsdffile = 'hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat';
n_MLA_grid = 1.30:0.05:2.00;           % slices of BSDF_MLA

lam      = (400:700)';                 % wavelength (nm), column
spec_name = 'I_Irppy2acac';            % emission spectrum in the library
Theta    = 2/3;                        % horizontal dipole ratio
PLQY     = 1;                          % radiative efficiency eta_rad
n_sub    = 1.8;                        % substrate (= MLA) index
N_pass   = 60;                         % number of passes in the series (j = 1 is the direct pass)

% stack from the reflector (top) to the substrate (bottom); EML is one of the layers
%            name                     n_o field                 n_e field              d (nm)
stack = { 'Ag',                      'l_Ag_McPeak',           'l_Ag_McPeak',          100 ;
          'B3PyMPM (ETL)',           'l_B3_o_JO',             'l_B3_e_JO',            300 ;
          'TCTA:B3PyMPM (EML)',      'l_TCTA_B3_o_JO',        'l_TCTA_B3_e_JO',        25 ;
          'TAPC (HTL)',              'l_TAPC_o_JO',           'l_TAPC_e_JO',          300 ;
          'ITO',                     'l_ITO',                 'l_ITO',                 50 };
iEML = 3;                              % row of the EML in 'stack'
z0   = 12.5;                           % dipole position, measured from the EML/ETL interface (nm)

%% ------------------------------------------------------------------ inputs
L = load(nkfile); mat = L.material; spec = L.spectrum;
Nl = numel(lam); iw = lam - 399;      % library grid starts at 400 nm, 1 nm steps
S  = spec.(spec_name); S = S(iw); S = S(:) / sum(S);

nL = size(stack,1);
no = zeros(Nl, nL+2); ne = zeros(Nl, nL+2);          % [air, stack..., substrate]
no(:,1) = 1; ne(:,1) = 1;
for k = 1:nL
    a = mat.(stack{k,2}); b = mat.(stack{k,3});
    no(:,k+1) = a(iw); ne(:,k+1) = b(iw);
end
no(:,end) = n_sub; ne(:,end) = n_sub;
d  = cell2mat(stack(:,4))';            % layer thicknesses
E  = iEML + 1;                         % EML column in no/ne (column 1 is air)
NC = nL + 2;                           % total columns

B = load(bsdffile); B = B.BSDF_MLA;
[~, islice] = min(abs(n_MLA_grid - n_sub));
BSDF = B(:,:,islice);
fprintf('BSDF slice %d (n_MLA = %.2f); column sums %.4f..%.4f\n', islice, n_MLA_grid(islice), min(sum(BSDF)), max(sum(BSDF)));

sin089 = sind(0:89);                   % angle grid of planar_Sweep.m (0:89 deg); BSDF bin j <-> grid point j
th  = (0:89)';
s089 = sin089; s089(1) = sind(0.5);    % solid-angle weight without 0/0 at the normal
B_T = sum(BSDF(1:90,:), 1);            % angle-resolved escape probability B_T(theta_i), 1 x 90
B_R = BSDF(180:-1:91, :);              % returned power: row = returned angle, column = incident angle
B_T_air = BSDF(1:90,:) ./ sind((0.5:1:89.5)');   % air intensity per solid angle: row = air angle, column = incident angle

%% ------------------------------------------------------------------ dipole model (CPS)
Nu = 1000; umax = 3;
u  = [(0:Nu-1)/Nu, (Nu+1:Nu*umax)/Nu];  % in-plane wavevector / (k0 n_EML)
nu = numel(u);

% below the dipole: EML part, layers towards the substrate, substrate
bot = E:NC;  dbot = [d(iEML) - z0, d(iEML+1:end), 0];
% above the dipole: EML part, layers towards the reflector, air
top = E:-1:1; dtop = [z0, d(iEML-1:-1:1), 0];
TB = TMF_birefringence_whole(no(:,bot), ne(:,bot), dbot, u, lam);
TT = TMF_birefringence_whole(no(:,top), ne(:,top), dtop, u, lam);

noE = no(:,E); neE = ne(:,E);
den_p = 1 - TB.r_p .* TT.r_p;  den_s = 1 - TB.r_s .* TT.r_s;
Kpv = 3/4*real((neE./noE) * (u.^2 ./ sqrt(1-u.^2)) .* (1+TB.r_p).*(1+TT.r_p) ./ den_p);
Kph = real(3./(6*(noE./neE).^2+2) * sqrt(1-u.^2) .* (1-TB.r_p).*(1-TT.r_p) ./ den_p);
Ksh = real(3./(2*(neE./noE).^2+6) * (1./sqrt(1-u.^2)) .* (1+TB.r_s).*(1+TT.r_s) ./ den_s);

cv = (1-Theta) * neE ./ lam.^4;                              % vertical-dipole weight
ch = Theta * noE .* (3 + (neE./noE).^2) ./ (4*lam.^4);        % horizontal-dipole weight
Up = 2*((cv*u).*Kpv + (ch*u).*Kph);  Us = 2*(ch*u).*Ksh;      % total dissipation per u, p and s
Utot = sum(Up + Us, 2);

% power transmitted into the (semi-infinite) substrate, per u
Tp = zeros(Nl, nu); Ts = zeros(Nl, nu); isub_p = zeros(Nl,1); isub_s = zeros(Nl,1);
for i = 1:Nl
    for pol = 'ps'
        if pol == 'p', nref = neE(i); ns = ne(i,end); else, nref = noE(i); ns = no(i,end); end
        im = ceil(Nu*ns/nref) - (ns > nref);                  % last u index that propagates in the substrate
        uu = u(1:im);
        ph = ones(1, im);                                     % near field through the rest of the EML (u > 1)
        if im > Nu
            ph(Nu+1:im) = exp(-4*pi*noE(i)*sqrt(uu(Nu+1:im).^2-1)*(d(iEML)-z0)/lam(i));
        end
        if pol == 'p'
            A1 = abs((1+TT.r_p(i,1:im)).*TB.t_p(i,1:im)./den_p(i,1:im)).^2;
            A2 = abs((1-TT.r_p(i,1:im)).*TB.t_p(i,1:im)./den_p(i,1:im)).^2;
            kv = 3/8*neE(i)*no(i,end)/noE(i)^2 * sqrt(1-(neE(i)*uu/ne(i,end)).^2) .* ph .* uu.^2 .* A1 ./ abs(1-uu.^2);
            kh = 3*sqrt((no(i,end)/noE(i))^2*(1-(neE(i)*uu/ne(i,end)).^2)) .* ph .* A2 / (12*(noE(i)/neE(i))^2+4);
            Tp(i,1:im) = 2*(cv(i)*uu.*kv + ch(i)*uu.*kh); isub_p(i) = im;
        else
            A3 = abs((1+TT.r_s(i,1:im)).*TB.t_s(i,1:im)./den_s(i,1:im)).^2;
            ks = 3*sqrt((no(i,end)/noE(i))^2 - uu.^2) .* ph .* A3 ./ ((4*(neE(i)/noE(i))^2+12)*abs(1-uu.^2));
            Ts(i,1:im) = 2*ch(i)*uu.*ks; isub_s(i) = im;
        end
    end
end

% channel budget per wavelength (fractions of Utot)
ch_sub = zeros(Nl,1); ch_wg = zeros(Nl,1); ch_spp = zeros(Nl,1); ch_abs = zeros(Nl,1);
for i = 1:Nl
    ip = 1:isub_p(i); is = 1:isub_s(i);
    ch_sub(i) = sum(Tp(i,ip)) + sum(Ts(i,is));
    ch_abs(i) = sum(Up(i,ip)) + sum(Us(i,is)) - ch_sub(i);
    ch_wg(i)  = sum(Up(i,isub_p(i)+1:Nu)) + sum(Us(i,isub_s(i)+1:Nu));
    ch_spp(i) = sum(Up(i,max(isub_p(i),Nu)+1:end)) + sum(Us(i,max(isub_s(i),Nu)+1:end));
end
ch_sub = ch_sub./Utot; ch_abs = ch_abs./Utot; ch_wg = ch_wg./Utot; ch_spp = ch_spp./Utot;

% Purcell factor and effective radiative efficiency
Kfree = 3/4*(neE./noE)*(u(1:Nu).^2./sqrt(1-u(1:Nu).^2)) .* (cv*u(1:Nu))*2 ...
      + 2*(ch*u(1:Nu)).*(3./(6*(noE./neE).^2+2)*sqrt(1-u(1:Nu).^2) + 3./(2*(neE./noE).^2+6)*(1./sqrt(1-u(1:Nu).^2)));
F = Utot ./ sum(Kfree, 2);
eta_eff = PLQY*F ./ (1 - PLQY + PLQY*F);

w_ph = lam.*S / sum(lam.*S);           % photon-number weight of the spectrum
P0   = w_ph .* eta_eff .* ch_sub;      % substrate-delivered power per wavelength (EQE units)
eta_sub = sum(P0);
budget  = [sum(w_ph.*eta_eff.*ch_sub), sum(w_ph.*eta_eff.*ch_wg), sum(w_ph.*eta_eff.*ch_spp), sum(w_ph.*eta_eff.*ch_abs)];

%% ------------------------------------------------------------------ I_sub(theta) and R_LED(theta)  (planar_Sweep.m definitions)
% intensity per solid angle in the substrate, per wavelength: P_sub as in planar_Sweep.m, then
% I_sub = P_sub .* (emission_spectrum .* eta_eff ./ Purcell)
uu0 = u; uu0(1) = 1;                                    % u = 0 column is 0 anyway
K_T_p = Tp ./ (2*uu0);  K_T_s = Ts ./ (2*uu0);          % = const.*K_v2 + const2.*K_h2 of planar_Sweep (no factor 2, no u)
c_free = pi*(cv + ch);                                   % pi * (free-space dissipation per unit u)
P_sub = zeros(Nl, 90);
for i = 1:Nl
    rp = ne(i,end)/neE(i); ip = 1:isub_p(i);
    Pp = rp*spline(neE(i)*u(ip), sqrt(max(0, rp^2 - u(ip).^2)).*K_T_p(i,ip), ne(i,end)*sin089)/c_free(i);
    rs = no(i,end)/noE(i); is = 1:isub_s(i);
    Ps = rs*spline(noE(i)*u(is), sqrt(max(0, rs^2 - u(is).^2)).*K_T_s(i,is), no(i,end)*sin089)/c_free(i);
    P_sub(i,:) = Pp + Ps;
end
P_sub(P_sub < 0) = 0;
I_sub = P_sub .* repmat(S .* eta_eff ./ F, 1, 90);        % Nl x 90, same as planar_Sweep.m

% reflectance of the whole OLED seen from inside the substrate (unpolarised), on sin089
allc = NC:-1:1; dall = [0, d(end:-1:1), 0];
Rp = TMF_birefringence_whole_p(no(:,allc), ne(:,allc), dall, ne(:,end)*sin089, lam);
Rs = TMF_birefringence_whole_s(no(:,allc), ne(:,allc), dall, no(:,end)*sin089, lam);
R_LED = (abs(Rp.r_p).^2 + abs(Rs.r_s).^2)/2;           % Nl x 90

%% ------------------------------------------------------------------ multi-pass series (planar_Sweep_passes logic)
% w = I_sub .* sin(theta) is power per unit polar angle; each pass: w <- R_LED(theta_r) .* (B_R * w)
% I_sub_k = w ./ sin(theta) is the intensity of the light arriving at pass k (k = 1 is I_sub),
% I_air_k = B_T_air * w is the intensity in air of the light extracted at pass k.
% The EQE bookkeeping uses the photon-number weight: each wavelength is rescaled so that its
% substrate-delivered power equals P0(lambda); the intensities keep the I_sub units.
I_sub_k = zeros(Nl, 90, N_pass);
I_air_k = zeros(Nl, 90, N_pass);
V_arr  = zeros(90, N_pass);   % power per degree arriving at the structure, pass j (EQE units)
V_ret  = zeros(90, N_pass);   % returned by the structure, before the OLED reflection
Pout   = zeros(1, N_pass);    % extracted power, pass j (EQE units)
for i = 1:Nl
    w = I_sub(i,:) .* s089;                          % I_sub units x sin
    c = P0(i) / max(sum(w), realmin);                % -> EQE units for the bookkeeping
    for j = 1:N_pass
        I_sub_k(i,:,j) = w ./ s089;
        I_air_k(i,:,j) = (B_T_air * w.').';
        V_arr(:,j) = V_arr(:,j) + c*w.';
        Pout(j)    = Pout(j) + c*(B_T*w.');
        r = (B_R * w.').';  V_ret(:,j) = V_ret(:,j) + c*r.';
        w = R_LED(i,:) .* r;                          % R_LED at the returned angle
    end
end
I_sub_k(:,1,2:end) = I_sub_k(:,2,2:end);             % normal point for passes >= 2 taken from 1 deg
I_sub_k_total = squeeze(sum(I_sub_k, 1));            % 90 x N_pass, summed over wavelength
I_air_k_total = squeeze(sum(I_air_k, 1));
I_air = I_air_k_total;                               % kept for the plots below

U    = sum(V_arr, 1);                                  % power arriving at pass j (= remaining eta_sub)
p_j  = Pout ./ U;                                      % escape probability of pass j
A_j  = 1 - sum(V_arr(:,2:end),1) ./ sum(V_ret(:,1:end-1),1);   % round-trip loss after pass j
eta_ext = sum(Pout) / eta_sub;
EQE     = sum(Pout);
fprintf('check pass 1: max|I_sub_k(:,:,1) - I_sub| = %.2e\n', max(max(abs(I_sub_k(:,2:end,1) - I_sub(:,2:end)))));

% eq. (3) with the Lambertian weight
wL = cosd(th).*sind(th); wL = wL/sum(wL);
p_lamb = B_T*wL;
A_lamb = 1 - sum(w_ph.*eta_eff.*ch_sub .* (R_LED*wL)) / eta_sub;   % spectrum-weighted
eta_ext_eq3 = p_lamb / (p_lamb + (1-p_lamb)*A_lamb);
p_eff = sum(Pout)/sum(U);                              % arrival-weighted mean escape probability
A_eff = 1 - sum(sum(V_arr(:,2:end)))/sum(sum(V_ret(:,1:end-1)));   % return-weighted mean round-trip loss

%% ------------------------------------------------------------------ report
fprintf('\neta_sub = %.4f  (waveguide %.4f, SPP/evanescent %.4f, absorption %.4f)\n', budget);
fprintf('eta_ext  series = %.4f   eq.(3) Lambertian = %.4f   (p_Lamb %.4f, A''_Lamb %.4f)\n', eta_ext, eta_ext_eq3, p_lamb, A_lamb);
fprintf('EQE      series = %.4f   eq.(3) = %.4f\n', EQE, eta_sub*eta_ext_eq3);
fprintf('arrival-weighted p = %.4f (%.1f%% of p_Lamb), return-weighted A'' = %.4f -> eq.(3) with both: %.4f\n', p_eff, 100*p_eff/p_lamb, A_eff, p_eff/(p_eff+(1-p_eff)*A_eff));
fprintf('\n pass  arriving  p_j     p_j/p_Lamb  A''_j    frac>60deg\n');
f60 = sum(V_arr(th>60,:),1) ./ U;
for j = [1:6 8 11 16 21 31]
    fprintf('%4d   %.4f   %.4f   %.3f      %.4f  %.3f\n', j, U(j), p_j(j), p_j(j)/p_lamb, A_j(min(j,N_pass-1)), f60(j));
end

save('MLA_pass_resolved_out.mat', 'th', 'lam', 'eta_sub', 'budget', 'eta_ext', 'EQE', 'eta_ext_eq3', 'p_lamb', 'A_lamb', 'p_eff', 'A_eff', ...
     'U', 'p_j', 'A_j', 'Pout', 'V_arr', 'V_ret', 'I_sub', 'I_sub_k_total', 'I_air_k_total', 'R_LED', 'B_T', 'P0');
csvwrite('MLA_pass_table.csv', [(1:N_pass)', U', p_j', (p_j/p_lamb)', [A_j NaN]', Pout', cumsum(Pout)'/eta_sub, f60']);
csvwrite('MLA_I_sub_k_total.csv', [th I_sub_k_total]);
csvwrite('MLA_I_air_k_total.csv', [th I_air_k_total sum(I_air_k_total,2)]);

%% ------------------------------------------------------------------ plots
if ~exist('OCTAVE_VERSION', 'builtin')
    ks  = [1 2 3 6 11 21];
    ang = deg2rad([-flipud(th); th]);
    cols = parula(numel(ks)+1);
    IL = cosd(th);                                      % Lambertian, 1 at the normal

    figure('Name','substrate intensity');
    pax = polaraxes; hold(pax,'on');
    for n = 1:numel(ks)
        r = I_sub_k_total(:,ks(n)) / I_sub_k_total(1,1);   % normalised to I_sub at 0 deg (absolute decay visible)
        polarplot(pax, ang, [flipud(r); r], 'Color', cols(n,:), 'LineWidth', 1.4, 'DisplayName', sprintf('pass %d', ks(n)));
    end
    polarplot(pax, ang, [flipud(IL); IL], 'k--', 'DisplayName', 'Lambertian');
    pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise'; pax.ThetaLim = [-90 90];
    legend('Location','southoutside'); title('Light arriving at the outcoupling structure (in the substrate)');

    figure('Name','air intensity');
    pax = polaraxes; hold(pax,'on');
    for n = 1:numel(ks)
        r = I_air(:,ks(n)) / max(I_air(:,1));
        polarplot(pax, ang, [flipud(r); r], 'Color', cols(n,:), 'LineWidth', 1.4, 'DisplayName', sprintf('pass %d', ks(n)));
    end
    r = sum(I_air,2)/max(sum(I_air,2));
    polarplot(pax, ang, [flipud(r); r], 'k', 'LineWidth', 2, 'DisplayName', 'total (normalised)');
    pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise'; pax.ThetaLim = [-90 90];
    legend('Location','southoutside'); title('Extracted light (in air)');

    figure('Name','escape probability');
    plot(1:30, p_j(1:30), 'o-', 1:30, p_lamb*ones(1,30), 'k--'); grid on
    xlabel('pass j'); ylabel('escape probability'); legend('p_j (series)', 'p, Lambertian');
end

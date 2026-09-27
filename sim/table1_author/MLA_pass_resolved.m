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
%   5. All angular quantities on the BSDF bin centres, 0.5 : 1 : 89.5 deg.
%
% Needs in the path: nk_JH_total.mat, the BSDF file, TMF_birefringence_whole(_p,_s).m
% (the version with the evanescent-branch fix, imag(n cos) >= 0).
% Runs in MATLAB and Octave (the polar plots are skipped in Octave).

if exist('cfg','var') ~= 1, cfg = struct(); end   % optional overrides, see below
clc;

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
bsdf_mode = 'MLA';                     % 'MLA' (BSDF file) or 'flat' (planar substrate/air interface, Fresnel)
n_top    = 1.0;                        % medium above the top layer
tag      = '';                         % suffix for the output files

% stack from the reflector (top) to the substrate (bottom); EML is one of the layers
%            name                     n_o field                 n_e field              d (nm)
stack = { 'Ag',                      'l_Ag_McPeak',           'l_Ag_McPeak',          100 ;
          'B3PyMPM (ETL)',           'l_B3_o_JO',             'l_B3_e_JO',            300 ;
          'TCTA:B3PyMPM (EML)',      'l_TCTA_B3_o_JO',        'l_TCTA_B3_e_JO',        25 ;
          'TAPC (HTL)',              'l_TAPC_o_JO',           'l_TAPC_e_JO',          300 ;
          'ITO',                     'l_ITO',                 'l_ITO',                 50 };
iEML = 3;                              % row of the EML in 'stack'
z0   = 12.5;                           % dipole position, measured from the EML/ETL interface (nm)
% a field of cfg overrides the setting of the same name (e.g. cfg.n_sub = 1.5; cfg.bsdf_mode = 'flat';)
fn = fieldnames(cfg); for q = 1:numel(fn), eval([fn{q} ' = cfg.(fn{q});']); end

%% ------------------------------------------------------------------ inputs
L = load(nkfile); mat = L.material; spec = L.spectrum;
Nl = numel(lam); iw = lam - 399;      % library grid starts at 400 nm, 1 nm steps
S  = spec.(spec_name); S = S(iw); S = S(:) / sum(S);

nL = size(stack,1);
no = zeros(Nl, nL+2); ne = zeros(Nl, nL+2);          % [air, stack..., substrate]
no(:,1) = n_top; ne(:,1) = n_top;
for k = 1:nL
    if ischar(stack{k,2}), a = mat.(stack{k,2}); a = a(iw); else, a = stack{k,2}*ones(Nl,1); end   % library field or constant
    if ischar(stack{k,3}), b = mat.(stack{k,3}); b = b(iw); else, b = stack{k,3}*ones(Nl,1); end
    no(:,k+1) = a; ne(:,k+1) = b;
end
no(:,end) = n_sub; ne(:,end) = n_sub;
d  = cell2mat(stack(:,4))';            % layer thicknesses
E  = iEML + 1;                         % EML column in no/ne (column 1 is air)
NC = nL + 2;                           % total columns

th  = (0.5:1:89.5)';                   % angle bins in the substrate / air (deg)
if strcmp(bsdf_mode, 'MLA')
    B = load(bsdffile); B = B.BSDF_MLA;
    [~, islice] = min(abs(n_MLA_grid - n_sub));
    BSDF = B(:,:,islice);
    fprintf('BSDF slice %d (n_MLA = %.2f); column sums %.4f..%.4f\n', islice, n_MLA_grid(islice), min(sum(BSDF)), max(sum(BSDF)));
else                                   % planar substrate/air interface: Fresnel, specular
    BSDF = zeros(180, 90);
    for a = 1:90
        ci = cosd(th(a)); st = n_sub*sind(th(a));
        if st < 1
            ct = sqrt(1-st^2);
            rs = ((n_sub*ci-ct)/(n_sub*ci+ct))^2; rp = ((ci-n_sub*ct)/(ci+n_sub*ct))^2; Rf = (rs+rp)/2;
            % the 1-deg substrate bin maps to a wider air range near the critical angle: spread it
            lo = n_sub*sind(th(a)-0.5); hi = min(1, n_sub*sind(th(a)+0.5));
            e1 = asind(lo); e2 = asind(hi);
            for bb = floor(e1)+1 : min(90, ceil(e2))
                ov = max(0, min(bb, e2) - max(bb-1, e1));
                BSDF(bb, a) = BSDF(bb, a) + (1 - Rf) * ov/(e2 - e1);
            end
        else
            Rf = 1;
        end
        BSDF(181-a, a) = Rf;           % specular: returned at the same angle (row 181-a <-> angle a)
    end
    fprintf('flat substrate/air interface, n_sub = %.2f\n', n_sub);
end
B_T = sum(BSDF(1:90,:), 1);            % angle-resolved escape probability B_T(theta_i), 1 x 90
B_R = BSDF(180:-1:91, :);              % returned power: row = returned angle, column = incident angle
B_T_air = BSDF(1:90,:) ./ sind(th);    % air intensity per solid angle: row = air angle, column = incident angle

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

%% ------------------------------------------------------------------ P_sub(theta) and R_LED(theta)
% power per unit polar angle in the substrate: dP/dtheta = (dP/du)(du/dtheta)
Psub = zeros(Nl, 90);
for i = 1:Nl
    for pol = 'ps'
        if pol == 'p', nref = neE(i); ns = ne(i,end); T = Tp; im = isub_p(i);
        else,          nref = noE(i); ns = no(i,end); T = Ts; im = isub_s(i); end
        ut  = ns*sind(th)'/nref;                        % u of each substrate angle
        dud = ns*cosd(th)'/nref * pi/180;               % du/dtheta per degree
        Psub(i,:) = Psub(i,:) + interp1(u(1:im), T(i,1:im)*Nu, ut, 'linear', 0) .* dud;   % T*Nu = density per unit u
    end
end
Psub = Psub ./ sum(Psub, 2);                            % angular shape per wavelength (sum = 1)
Psub(~isfinite(Psub)) = 0;

% reflectance of the whole OLED seen from inside the substrate (unpolarised)
allc = NC:-1:1; dall = [0, d(end:-1:1), 0];
Rp = TMF_birefringence_whole_p(no(:,allc), ne(:,allc), dall, ne(:,end)*sind(th)', lam);
Rs = TMF_birefringence_whole_s(no(:,allc), ne(:,allc), dall, no(:,end)*sind(th)', lam);
R_LED = (abs(Rp.r_p).^2 + abs(Rs.r_s).^2)/2;           % Nl x 90

%% ------------------------------------------------------------------ multi-pass series
V_arr  = zeros(90, N_pass);   % power per degree arriving at the structure, pass j (EQE units)
V_ret  = zeros(90, N_pass);   % returned by the structure, before the OLED reflection
I_air  = zeros(90, N_pass);   % extracted light, intensity per solid angle in air (EQE units / sr, arb.)
Pout   = zeros(1, N_pass);    % extracted power, pass j
for i = 1:Nl
    v = P0(i) * Psub(i,:)';                            % column, EQE units
    for j = 1:N_pass
        V_arr(:,j) = V_arr(:,j) + v;
        Pout(j)    = Pout(j) + B_T*v;
        I_air(:,j) = I_air(:,j) + B_T_air*v;
        r = B_R*v;  V_ret(:,j) = V_ret(:,j) + r;
        v = R_LED(i,:)' .* r;                          % R_LED at the returned angle
    end
end
U    = sum(V_arr, 1);                                  % power arriving at pass j (= remaining eta_sub)
p_j  = Pout ./ U;                                      % escape probability of pass j
A_j  = 1 - sum(V_arr(:,2:end),1) ./ sum(V_ret(:,1:end-1),1);   % round-trip loss after pass j
eta_ext = sum(Pout) / eta_sub;
EQE     = sum(Pout);

% intensity per solid angle in the substrate of the light arriving at pass j
I_sub = V_arr ./ sind(th);                             % EQE units / sr (per degree bin)

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

save(['MLA_pass_resolved_out' tag '.mat'], 'th', 'lam', 'eta_sub', 'budget', 'eta_ext', 'EQE', 'eta_ext_eq3', 'p_lamb', 'A_lamb', 'p_eff', 'A_eff', ...
     'U', 'p_j', 'A_j', 'Pout', 'V_arr', 'V_ret', 'I_sub', 'I_air', 'Psub', 'R_LED', 'B_T', 'P0');
csvwrite(['MLA_pass_table' tag '.csv'], [(1:N_pass)', U', p_j', (p_j/p_lamb)', [A_j NaN]', Pout', cumsum(Pout)'/eta_sub, f60']);
csvwrite(['MLA_I_sub' tag '.csv'], [th I_sub]);
csvwrite(['MLA_I_air' tag '.csv'], [th I_air sum(I_air,2)]);

%% ------------------------------------------------------------------ plots
if ~exist('OCTAVE_VERSION', 'builtin')
    ks  = [1 2 3 6 11 21];
    ang = deg2rad([-flipud(th); th]);
    cols = parula(numel(ks)+1);
    IL = cosd(th); IL = IL/sum(IL.*sind(th));          % Lambertian, same normalisation

    figure('Name','substrate intensity');
    pax = polaraxes; hold(pax,'on');
    for n = 1:numel(ks)
        r = I_sub(:,ks(n)) / U(ks(n));                 % shape of pass j (area-normalised)
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

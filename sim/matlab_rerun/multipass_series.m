function [EQE, eta_ext, n_used, tail] = multipass_series(Psub_norm, R_OLED, BSDF, w_sub, tol, n_max)
% MULTIPASS_SERIES  Light escaping through an outcoupling structure after
% repeated round trips between the structure and the OLED (trans-scale model).
%
%   Psub_norm : [n_lambda x 90]  angular distribution of the power entering the
%               substrate (power per 1-degree polar bin, each row sums to 1)
%   R_OLED    : [n_lambda x 90]  reflectance of the OLED stack seen from the substrate
%   BSDF      : [180 x 90]       BSDF of the outcoupling structure, columns = angle of
%               incidence (0..89 deg); rows 1:90 = transmitted into air, rows 91:180
%               = reflected back into the substrate (row 180 = 0 deg)
%   w_sub     : [n_lambda x 1]   photons delivered to the substrate per exciton at each
%               wavelength (sum(w_sub) = eta_sub)
%   tol       : stop when the last term is < tol of the running sum (default 1e-6)
%   n_max     : maximum number of round trips (default 2000)
%
%   EQE       = sum_lambda w_sub * eta_ext(lambda)
%   eta_ext   = EQE / sum(w_sub)
%   n_used    = number of round trips actually summed (largest over wavelengths)
%   tail      = largest (last term / sum) over wavelengths
if nargin < 5 || isempty(tol),   tol = 1e-6; end
if nargin < 6 || isempty(n_max), n_max = 2000; end
BSDF_R = BSDF(180:-1:91,:);          % rows: angle of reflection, columns: angle of incidence
BSDF_T = sum(BSDF(1:90,:), 1);       % total transmission for each angle of incidence
nl = size(Psub_norm,1);
ext = zeros(nl,1); n_used = 0; tail = 0;
for i = 1:nl
    M = BSDF_R' .* repmat(R_OLED(i,:), 90, 1);   % one round trip (row-vector convention v = v*M)
    v = Psub_norm(i,:); t = 0;
    for j = 1:n_max
        term = v*BSDF_T';
        t = t + term;
        v = v*M;
        if t > 0 && term < tol*t, break; end
    end
    ext(i) = t; n_used = max(n_used, j);
    if t > 0, tail = max(tail, term/t); end
end
EQE = sum(w_sub(:).*ext);
eta_ext = EQE/sum(w_sub);
end

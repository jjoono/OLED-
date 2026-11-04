% Extended Data Fig. 7d: EQE against ETL thickness for four ETL extinction coefficients.
% Stack (air side first): Ag (100 nm) / B3PyMPM (d_ETL, k_ETL added to n_o and n_e)
%   / TCTA:B3PyMPM (25 nm, emitter at the centre) / TAPC (200 nm) / ITO (50 nm)
%   / substrate n = 1.8 + close-packed hemispherical MLA (LightTools BSDF, n_MLA = 1.80)
% Ir(ppy)2acac spectrum, 400-700 nm in 1-nm steps, Theta = 2/3, PLQY = 1,
% series summed to convergence (last term < 1e-6 of the sum).
% Needs nk_JH_total.mat and the BSDF .mat on the path (both in ../table1_author).
clear; addpath('../table1_author');
load('nk_JH_total.mat'); load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat');
BSDF = BSDF_MLA(:,:,11);                         % n_MLA = 1.80
W = 1:301; wl = (400:700)'; nl = numel(wl); one = ones(nl,1);
spec = spectrum.I_Irppy2acac(W);
k_list = [0 3e-4 1e-3 3e-3]; d_list = 20:10:400;
if ~isempty(getenv('K_SEL')), k_list = str2num(getenv('K_SEL')); end   % optional: run one k per process
out = 'ED7d_etl_absorption.csv'; if ~isempty(getenv('OUT')), out = getenv('OUT'); end
fid = fopen(out,'w');
fprintf(fid,'k_ETL,d_ETL_nm,eta_sub,eta_ext,EQE,n_terms\n');
for k = k_list
  for d = d_list
    no_bar = [one material.l_Ag_McPeak(W) material.l_B3_o_JO(W)+1i*k material.l_TCTA_B3_o_JO(W) material.l_TAPC_o_JO(W) material.l_ITO(W) 1.8*one];
    ne_bar = [one material.l_Ag_McPeak(W) material.l_B3_e_JO(W)+1i*k material.l_TCTA_B3_e_JO(W) material.l_TAPC_e_JO(W) material.l_ITO(W) 1.8*one];
    o = planar_inputs(no_bar, ne_bar, [100 d 25 200 50], 4, 12.5, wl, spec, 2/3, 1);
    [EQE, ext, n] = multipass_series(o.Psub_norm, o.R_OLED, BSDF, o.w_sub);
    fprintf(fid,'%g,%d,%.6f,%.6f,%.6f,%d\n', k, d, o.eta_sub, ext, EQE, n);
    fprintf('k %.0e d %3d | eta_sub %.4f eta_ext %.4f EQE %.4f (%d terms)\n', k, d, o.eta_sub, ext, EQE, n);
  end
end
fclose(fid);

% Sanity check: planar_inputs + multipass_series (20 terms) must reproduce
% ../table1_author/table1_author_model.csv (d_ETL 300, k_ITO 0.0032, Theta 2/3, 400-700 nm).
% Required in the same folder (or on the MATLAB path):
%   planar_inputs.m, multipass_series.m, TMF_birefringence_whole.m, TMF_birefringence_whole_p.m,
%   TMF_birefringence_whole_s.m, nk_JH_total.mat, hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat
here = fileparts(mfilename('fullpath')); if ~isempty(here), addpath(here); end
if exist(fullfile(here,'..','table1_author'),'dir'), addpath(fullfile(here,'..','table1_author')); end
req = {'planar_inputs','multipass_series','TMF_birefringence_whole','TMF_birefringence_whole_p','TMF_birefringence_whole_s','nk_JH_total.mat','hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat'};
miss = req(cellfun(@(f) ~exist(f,'file'), req)); if ~isempty(miss), error('missing: %s', strjoin(miss, ', ')); end

load('nk_JH_total.mat'); load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat');
W = 1:301; wl = (400:700)'; nl = numel(wl);
ito = material.l_ITO(W); ito = ito(:);
no_bar = [ones(nl,1) material.l_Ag_McPeak(W) material.l_B3_o_JO(W) material.l_TCTA_B3_o_JO(W) material.l_TAPC_o_JO(W) ito 1.8*ones(nl,1)];
ne_bar = [ones(nl,1) material.l_Ag_McPeak(W) material.l_B3_e_JO(W) material.l_TCTA_B3_e_JO(W) material.l_TAPC_e_JO(W) ito 1.8*ones(nl,1)];
o = planar_inputs(no_bar, ne_bar, [100 300 25 180 50], 4, 12.5, wl, spectrum.I_Irppy2acac(W), 2/3, 1);
[E20, x20] = multipass_series(o.Psub_norm, o.R_OLED, BSDF_MLA(:,:,11), o.w_sub, 0, 20);
[Ec, xc, n] = multipass_series(o.Psub_norm, o.R_OLED, BSDF_MLA(:,:,11), o.w_sub);
fprintf('eta_sub %.4f | 20 terms: eta_ext %.4f EQE %.4f | converged (%d terms): EQE %.4f\n', o.eta_sub, x20, E20, n, Ec);

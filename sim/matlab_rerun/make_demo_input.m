% Builds release_minimal/demo_input.mat: the inputs of the multi-pass series for one
% device of the paper (Ag / B3PyMPM 300 nm / TCTA:B3PyMPM 25 nm / TAPC 180 nm / ITO 50 nm
% / n = 1.8 substrate + hemispherical MLA, Ir(ppy)2acac, Theta = 2/3, PLQY = 1).
% Only derived quantities are saved (no refractive-index library, no dipole code).
% Required in the same folder (or on the MATLAB path):
%   planar_inputs.m, multipass_series.m, TMF_birefringence_whole.m, TMF_birefringence_whole_p.m,
%   TMF_birefringence_whole_s.m, nk_JH_total.mat, hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat
here = fileparts(mfilename('fullpath')); if ~isempty(here), addpath(here); end
if exist(fullfile(here,'..','table1_author'),'dir'), addpath(fullfile(here,'..','table1_author')); end
req = {'planar_inputs','multipass_series','TMF_birefringence_whole','TMF_birefringence_whole_p','TMF_birefringence_whole_s','nk_JH_total.mat','hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat'};
miss = req(cellfun(@(f) ~exist(f,'file'), req)); if ~isempty(miss), error('missing: %s', strjoin(miss, ', ')); end

load('nk_JH_total.mat'); load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat');
W = 1:301; wavelength = (400:700)'; nl = numel(wavelength); one = ones(nl,1);
no_bar = [one material.l_Ag_McPeak(W) material.l_B3_o_JO(W) material.l_TCTA_B3_o_JO(W) material.l_TAPC_o_JO(W) material.l_ITO(W) 1.8*one];
ne_bar = [one material.l_Ag_McPeak(W) material.l_B3_e_JO(W) material.l_TCTA_B3_e_JO(W) material.l_TAPC_e_JO(W) material.l_ITO(W) 1.8*one];
o = planar_inputs(no_bar, ne_bar, [100 300 25 180 50], 4, 12.5, wavelength, spectrum.I_Irppy2acac(W), 2/3, 1);
Psub_norm = o.Psub_norm; R_OLED = o.R_OLED; w_sub = o.w_sub; BSDF = BSDF_MLA(:,:,11);
save('-v7','release_minimal/demo_input.mat','wavelength','Psub_norm','R_OLED','w_sub','BSDF');

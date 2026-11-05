% Fig. 4d (right) near its maximum: Ag 100 / B3PyMPM d_ETL / TCTA:B3PyMPM 25 / TAPC d_HTL / ITO 50 / n 1.8 + MLA (slice 11)
% Ir(ppy)2acac 400-700 nm (1 nm), Theta 2/3, PLQY 1; EQE after 20 round trips and converged.
clear;
here = fileparts(mfilename('fullpath')); if ~isempty(here), addpath(here); end
if exist(fullfile(here,'..','table1_author'),'dir'), addpath(fullfile(here,'..','table1_author')); end
load('nk_JH_total.mat'); load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat'); BSDF = BSDF_MLA(:,:,11);
W = 1:301; wl = (400:700)'; one = ones(numel(wl),1); spec = spectrum.I_Irppy2acac(W);
E = str2num(getenv('ETL')); H = [100 150 200 250 300];
for d = E
  for h = H
    no = [one material.l_Ag_McPeak(W) material.l_B3_o_JO(W) material.l_TCTA_B3_o_JO(W) material.l_TAPC_o_JO(W) material.l_ITO(W) 1.8*one];
    ne = [one material.l_Ag_McPeak(W) material.l_B3_e_JO(W) material.l_TCTA_B3_e_JO(W) material.l_TAPC_e_JO(W) material.l_ITO(W) 1.8*one];
    o = planar_inputs(no, ne, [100 d 25 h 50], 4, 12.5, wl, spec, 2/3, 1);
    E20 = multipass_series(o.Psub_norm, o.R_OLED, BSDF, o.w_sub, 0, 20);
    [Ec, x, n] = multipass_series(o.Psub_norm, o.R_OLED, BSDF, o.w_sub);
    fprintf('ETL %d HTL %d | eta_sub %.4f | EQE20 %.4f | EQEconv %.4f (%d)\n', d, h, o.eta_sub, E20, Ec, n);
  end
end

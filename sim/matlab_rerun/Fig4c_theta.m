% Fig. 4c: EQE versus horizontal dipole ratio, 550 nm, PLQY 1, n = 1.8 organics, emitter at the EML centre.
%   low-loss    : Ag 100 / ETL 200 / EML 20 / HTL 200 / ITO 50 / n_sub 1.8 + hemispherical MLA (n 1.8)
%   conventional: Al 100 / ETL  80 / EML 20 / HTL 230 / ITO 50 / glass 1.5 + hemispherical MLA (n 1.5)
% BSDFs: original LightTools file (slice 11 = n 1.80, slice 5 = n 1.50) and the ED8 LightTools set (bsdf/AR1p00.mat, n 1.80).
clear;
here = fileparts(mfilename('fullpath')); if ~isempty(here), addpath(here); end
if exist(fullfile(here,'..','table1_author'),'dir'), addpath(fullfile(here,'..','table1_author')); end
load('nk_JH_total.mat'); L = load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat');
B18_old = L.BSDF_MLA(:,:,11); B15_old = L.BSDF_MLA(:,:,5);
S = load(fullfile(here,'bsdf','AR1p00.mat')); B18_new = S.BSDF;
j = 151; ito = material.l_ITO(j);
nb_ll = [1 material.l_Ag_McPeak(j) 1.8 1.8 1.8 ito 1.8];
nb_cv = [1 material.l_Al_JO(j)     1.8 1.8 1.8 ito 1.5];
TH = 0.5:0.01:1.0; out = zeros(numel(TH),5);
for k = 1:numel(TH)
  o1 = planar_inputs(nb_ll, nb_ll, [100 200 20 200 50], 4, 10, 550, 1, TH(k), 1);
  o2 = planar_inputs(nb_cv, nb_cv, [100  80 20 230 50], 4, 10, 550, 1, TH(k), 1);
  out(k,:) = [TH(k), multipass_series(o1.Psub_norm,o1.R_OLED,B18_old,o1.w_sub), multipass_series(o1.Psub_norm,o1.R_OLED,B18_new,o1.w_sub), ...
              multipass_series(o2.Psub_norm,o2.R_OLED,B15_old,o2.w_sub), o1.eta_sub];
end
fid = fopen('Fig4c_theta.csv','w'); fprintf(fid,'Theta,EQE_lowloss_oldBSDF,EQE_lowloss_ED8BSDF,EQE_conventional_oldBSDF_n1.5,eta_sub_lowloss\n');
fprintf(fid,'%.2f,%.5f,%.5f,%.5f,%.5f\n',out'); fclose(fid);
i = find(abs(TH-2/3)==min(abs(TH-2/3)));
fprintf('Theta 0.50: LL old %.4f new %.4f | conv %.4f\n', out(1,2), out(1,3), out(1,4));
fprintf('Theta %.2f: LL old %.4f new %.4f | conv %.4f\n', TH(i), out(i,2), out(i,3), out(i,4));
fprintf('Theta 1.00: LL old %.4f new %.4f | conv %.4f\n', out(end,2), out(end,3), out(end,4));

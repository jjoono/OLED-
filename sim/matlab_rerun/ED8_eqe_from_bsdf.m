% Extended Data Fig. 8: EQE against microlens aspect ratio and fill factor.
% Generic stack (air side first): Ag or Al (100 nm) / ETL n = 1.8 (200 nm)
%   / EML n = 1.8 (20 nm, emitter at the centre) / HTL n = 1.8 (200 nm)
%   / ITO (50 nm, 1.8636 + 0.0032i) / substrate n = 1.8 + microlens array (n = 1.8)
% 550 nm, Theta = 2/3 (isotropic), PLQY = 1, series summed to convergence.
%
% The BSDF of every lens geometry comes from LightTools (one .mat per row of
% ED8_lighttools_geometries.csv, variable BSDF, 180 x 90, same layout as
% hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat).  Put them in ./bsdf/.
clear;
% Required in the same folder (or on the MATLAB path):
%   planar_inputs.m, multipass_series.m, TMF_birefringence_whole.m, TMF_birefringence_whole_p.m,
%   TMF_birefringence_whole_s.m, nk_JH_total.mat, flat_bsdf.m, pss.m, csv2cell_simple.m, ED8_lighttools_geometries.csv
here = fileparts(mfilename('fullpath')); if ~isempty(here), addpath(here); end
if exist(fullfile(here,'..','table1_author'),'dir'), addpath(fullfile(here,'..','table1_author')); end
req = {'planar_inputs','multipass_series','TMF_birefringence_whole','TMF_birefringence_whole_p','TMF_birefringence_whole_s','nk_JH_total.mat','flat_bsdf','pss','csv2cell_simple'};
miss = req(cellfun(@(f) ~exist(f,'file'), req)); if ~isempty(miss), error('missing: %s', strjoin(miss, ', ')); end

load('nk_JH_total.mat');
j = 151; wl = 550;                                       % 550 nm
refl.Ag = material.l_Ag_McPeak(j); refl.Al = material.l_Al_JO(j); ito = material.l_ITO(j);
G = csv2cell_simple('ED8_lighttools_geometries.csv');   % id, family, shape, aspect_ratio, fill_factor, ...
names = {'Ag','Al'}; P = struct();
for q = 1:2
  m = refl.(names{q});
  nb = [1 m 1.8 1.8 1.8 ito 1.8];
  P.(names{q}) = planar_inputs(nb, nb, [100 200 20 200 50], 4, 10, wl, 1, 2/3, 1);
end
fid = fopen('ED8_eqe.csv','w');
fprintf(fid,'id,family,shape,aspect_ratio,fill_factor,p_single_pass,eta_ext_Al,EQE_Al,eta_ext_Ag,EQE_Ag\n');
% planar reference (no microlenses)
Bf = flat_bsdf(1.8); row = sprintf('planar,reference,flat,0,0,%.6f', pss(Bf));
for q = [2 1]
  o = P.(names{q}); [E, x] = multipass_series(o.Psub_norm, o.R_OLED, Bf, o.w_sub); row = [row sprintf(',%.6f,%.6f', x, E)];
end
fprintf(fid,'%s\n',row); disp(row);
for r = 1:size(G,1)
  f = fullfile('bsdf', [G{r,1} '.mat']);
  if ~exist(f,'file'), fprintf('missing %s\n', f); continue; end
  S = load(f); B = S.BSDF;
  row = sprintf('%s,%s,%s,%s,%s,%.6f', G{r,1}, G{r,2}, G{r,3}, G{r,4}, G{r,5}, pss(B));
  for q = [2 1]
    o = P.(names{q}); [E, x] = multipass_series(o.Psub_norm, o.R_OLED, B, o.w_sub); row = [row sprintf(',%.6f,%.6f', x, E)];
  end
  fprintf(fid,'%s\n',row); disp(row);
end
fclose(fid);


% Extended Data Table 1 with the MATLAB model, series summed to convergence.
% Stack: Ag 100 / B3PyMPM d_ETL / TCTA:B3PyMPM 25 (emitter at centre) / TAPC 180 / ITO 50 (k scaled to k_ITO at 550 nm)
%        / substrate n = 1.8 + hemispherical MLA (LightTools BSDF slice n_MLA = 1.80).
clear;
here = fileparts(mfilename('fullpath')); if ~isempty(here), addpath(here); end
if exist(fullfile(here,'..','table1_author'),'dir'), addpath(fullfile(here,'..','table1_author')); end
load('nk_JH_total.mat'); load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat'); BSDF = BSDF_MLA(:,:,11);
fid = fopen('EDTable1_converged.csv','w'); fprintf(fid,'d_ETL,k_ITO,Theta,mode,eta_sub,spp,A_prime,eta_sta,EQE,n_terms\n');
w = cosd(0:89).*sind(0:89);
for mode = 1:2
  if mode == 1, W = 1:301; lab = 'full'; else, W = 151; lab = '550'; end
  wl = (400:800)'; wl = wl(W); nl = numel(wl); one = ones(nl,1); spec = spectrum.I_Irppy2acac(W);
  for d = [200 300]
    for kI = [0.0032 0.002]
      for th = [2/3 0.8 0.9]
        ito = material.l_ITO(W); ito = real(ito) + 1i*imag(ito)*(kI/imag(material.l_ITO(151)));
        no = [one material.l_Ag_McPeak(W) material.l_B3_o_JO(W) material.l_TCTA_B3_o_JO(W) material.l_TAPC_o_JO(W) ito 1.8*one];
        ne = [one material.l_Ag_McPeak(W) material.l_B3_e_JO(W) material.l_TCTA_B3_e_JO(W) material.l_TAPC_e_JO(W) ito 1.8*one];
        o = planar_inputs(no, ne, [100 d 25 180 50], 4, 12.5, wl, spec, th, 1);
        [E, x, n] = multipass_series(o.Psub_norm, o.R_OLED, BSDF, o.w_sub);
        if mode == 2, A = 1 - sum(o.R_OLED.*w)/sum(w); else, A = NaN; end
        fprintf(fid,'%d,%.4f,%.3f,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%d\n', d, kI, th, lab, o.eta_sub, o.spp, A, x, E, n);
        fprintf('%s d %d k %.4f Th %.3f | eta_sub %.4f spp %.4f A %.4f | eta_sta %.4f EQE %.4f (%d)\n', lab, d, kI, th, o.eta_sub, o.spp, A, x, E, n);
      end
    end
  end
end
fclose(fid);

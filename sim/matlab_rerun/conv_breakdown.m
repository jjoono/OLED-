clear; here=pwd; addpath(here); addpath(fullfile(here,'..','table1_author'));
load('nk_JH_total.mat'); L=load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat');
j=151; ito=material.l_ITO(j); w=cosd(0:89).*sind(0:89);
cases={'conv Al n1.5 80/230',[1 material.l_Al_JO(j) 1.8 1.8 1.8 ito 1.5],[100 80 20 230 50],L.BSDF_MLA(:,:,5),1.5;
       'conv Ag n1.5 80/230',[1 material.l_Ag_McPeak(j) 1.8 1.8 1.8 ito 1.5],[100 80 20 230 50],L.BSDF_MLA(:,:,5),1.5;
       'lowloss Ag n1.8 200/200',[1 material.l_Ag_McPeak(j) 1.8 1.8 1.8 ito 1.8],[100 200 20 200 50],L.BSDF_MLA(:,:,11),1.8};
for c=1:3
  for th=[2/3 1]
    nb=cases{c,2}; o=planar_inputs(nb,nb,cases{c,3},4,10,550,1,th,1);
    A=1-sum(o.R_OLED.*w)/sum(w); B=cases{c,4}; pT=sum(sum(B(1:90,:),1).*w)/sum(w);
    Eflat=multipass_series(o.Psub_norm,o.R_OLED,flat_bsdf(cases{c,5}),o.w_sub);
    Emla=multipass_series(o.Psub_norm,o.R_OLED,B,o.w_sub);
    fprintf('%-24s Th %.2f | eta_sub %.3f  A'' %.3f  p %.3f | EQE planar %.3f  EQE MLA %.3f  eta_sta %.3f | closed form %.3f\n',cases{c,1},th,o.eta_sub,A,pT,Eflat,Emla,Emla/o.eta_sub,o.eta_sub*pT/(pT+(1-pT)*A));
  end
end

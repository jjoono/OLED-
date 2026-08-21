% Contour plot of EQE vs (S, g) for the external scattering layer.
% Run after planar_Sweep_SL_JH_full_lambda.m (loads its saved result).
% The g = 1 column (pure forward scattering = no scattering at all) is a
% degenerate point equivalent to a bare substrate and is excluded.

load('SL_sweep_result.mat')

keep=g_list<0.999;
g_plot=g_list(keep);
E_plot=EQE_SL_map(:,keep);

figure;
contourf(g_plot,S_list,E_plot,linspace(EQE_air,EQE_sub,25),'LineStyle','none');
hold on;
[C,h]=contour(g_plot,S_list,E_plot,0.40:0.02:0.62,'LineColor',[0.12 0.23 0.37]);
clabel(C,h,'FontSize',8,'Color',[0.12 0.23 0.37]);
colormap(flipud(bone)); % single-hue, light->dark
caxis([EQE_air EQE_sub]); % absolute scale: planar EQE .. total substrate light
cb=colorbar;
ylabel(cb,'EQE (absolute)');
xlabel('g (Henyey-Greenstein asymmetry)');
ylabel('S (scattering optical thickness \mu_s d)');
title(sprintf('EQE with scattering layer (n_{SL}=n_{sub}=1.51), thick ETL/HTL fixed\nrange %.3f - %.3f (%.1f%% spread), planar %.3f, EQE_{sub} %.3f', ...
    min(E_plot(:)),max(E_plot(:)),100*(max(E_plot(:))-min(E_plot(:)))/max(E_plot(:)),EQE_air,EQE_sub));

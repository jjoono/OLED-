% run_demo: EQE of one OLED of the paper from the precomputed inputs in demo_input.mat.
load('demo_input.mat')
tic
[EQE20, ext20]       = multipass_series(Psub_norm, R_OLED, BSDF, w_sub, 0, 20);   % 20 round trips
[EQE, ext, n, tail]  = multipass_series(Psub_norm, R_OLED, BSDF, w_sub);          % to convergence
fprintf('eta_sub = %.4f\n', sum(w_sub));
fprintf('20 round trips : eta_ext = %.4f, EQE = %.4f\n', ext20, EQE20);
fprintf('converged (%d) : eta_ext = %.4f, EQE = %.4f (last term %.1e of the sum)\n', n, ext, EQE, tail);
toc

% EQE_MLA from the converged planar inputs plus your BSDF.
%   needs  mla_inputs.mat  (d_ETL, d_HTL, theta_deg, P0, P_sub, R_LED)
load('mla_inputs.mat')
load('hexagonal_half_sphere_MLA_BSDF_nMLA_130_0,05_200.mat')
BSDF = BSDF_MLA(:,:,5);          % <- the slice whose n_MLA matches the substrate (1.50)

R_step       = 20;
BSDF_R       = BSDF(180:-1:91,:);
BSDF_T_total = sum(BSDF(1:90,:));

nE = numel(d_ETL); nH = numel(d_HTL);
EQE_MLA  = zeros(nE,nH);
eta_ext  = zeros(nE,nH);
MLA_tail = zeros(nE,nH);
for a = 1:nE
    for b = 1:nH
        v = squeeze(P_sub(a,b,:))'.*sind(theta_deg);   % power per angular bin
        v = v/sum(v);                                  % = Psub_norm
        M = BSDF_R'.*squeeze(R_LED(a,b,:));            % round trip: BSDF then the mirror
        t = 0;
        for j = 1:R_step
            term = v*BSDF_T_total';
            t = t + term;
            v = v*M;
        end
        eta_ext(a,b)  = t;
        EQE_MLA(a,b)  = t*P0(a,b);
        MLA_tail(a,b) = term/t;                        % last term / total
    end
end
fprintf('EQE_MLA : %.4f .. %.4f    eta_ext : %.4f .. %.4f    worst tail : %.1e\n', ...
    min(EQE_MLA(:)),max(EQE_MLA(:)),min(eta_ext(:)),max(eta_ext(:)),max(MLA_tail(:)));
save('-v7','mla_eqe.mat','d_ETL','d_HTL','EQE_MLA','eta_ext','P0','MLA_tail');

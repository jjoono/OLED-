% Monte Carlo BSDF generator for an external scattering layer (SL)
% index-matched to the substrate (n_SL = n_sub).
%
% Model:
%   - Light is incident from the substrate at AOI bins 1..90 deg (bin centers 0.5..89.5).
%   - Substrate/SL interface: index matched -> no Fresnel reflection, no refraction.
%   - Inside the SL, photons scatter with a Henyey-Greenstein phase function
%     (asymmetry parameter g), no absorption (albedo = 1).
%   - S = scattering optical thickness of the layer (mu_s * d, i.e. the layer
%     thickness in units of the scattering mean free path).
%   - SL/air interface: unpolarized Fresnel reflection / TIR.
%   - Photons crossing back into the substrate are recorded as reflection.
%
% Output BSDF convention (identical to BSDF_MLA):
%   BSDF_SL(180, 90, iS, ig)
%     columns        : AOI bin in substrate (1..90 deg)
%     rows 1..90     : transmission into air, row = ceil(theta_air)
%     rows 91..180   : reflection into substrate, row = 181 - ceil(theta_sub)
%   Every column sums to 1 (lossless layer).
%
% Made for JH Kim's planar_Sweep pipeline.

clear;
tic;

n_SL=1.51;              % scattering layer index = substrate index
S_list=1:2:15;          % scattering optical thickness sweep
g_list=0.5:0.1:1;       % HG asymmetry sweep (g=1 -> pure forward = no scattering)
N_photon=1e5;           % photons per AOI bin
rng(1);                 % reproducibility (Octave: use rand('seed',1) if rng missing)

NS=length(S_list);
Ng=length(g_list);

BSDF_SL=zeros(180,90,NS,Ng);

for iS=1:NS
    S=S_list(iS);
    for ig=1:Ng
        g=g_list(ig);

        for aoi=1:90

            theta_i=(aoi-0.5)*pi/180;

            % direction cosines, +z = toward air; slab spans z in [0,S]
            ux=sin(theta_i)*ones(N_photon,1);
            uy=zeros(N_photon,1);
            uz=cos(theta_i)*ones(N_photon,1);
            z=zeros(N_photon,1);

            w_T=zeros(90,1);    % transmitted into air, per air-angle bin
            w_R=zeros(90,1);    % returned into substrate, per substrate-angle bin

            alive=true(N_photon,1);
            iter=0;

            while any(alive) && iter<20000

                iter=iter+1;
                idx=find(alive);
                n_act=length(idx);

                % free path to next scattering event (units of mean free path)
                l=-log(rand(n_act,1));
                z_new=z(idx)+l.*uz(idx);

                hit_top=(uz(idx)>0)&(z_new>=S);
                hit_bot=(uz(idx)<0)&(z_new<=0);
                scat=~(hit_top|hit_bot);

                % --- top boundary: SL/air Fresnel ---
                it=idx(hit_top);
                if ~isempty(it)
                    ci=uz(it);                          % cos(theta) inside SL
                    sin_t2=n_SL^2*(1-ci.^2);            % sin^2(theta_air)
                    R=ones(length(it),1);               % TIR by default
                    ok=sin_t2<1;
                    ct_air=sqrt(1-sin_t2(ok));
                    rs=(n_SL*ci(ok)-ct_air)./(n_SL*ci(ok)+ct_air);
                    rp=(n_SL*ct_air-ci(ok))./(n_SL*ct_air+ci(ok));
                    R(ok)=(rs.^2+rp.^2)/2;

                    refl=rand(length(it),1)<R;

                    % transmitted -> record air angle
                    tr=it(~refl);
                    if ~isempty(tr)
                        sin_air=n_SL*sqrt(1-uz(tr).^2);
                        ang=asind(min(sin_air,1));
                        bins=min(max(ceil(ang),1),90);
                        w_T=w_T+accumarray(bins,1,[90 1]);
                        alive(tr)=false;
                    end

                    % reflected -> flip direction at the boundary
                    rf=it(refl);
                    z(rf)=S;
                    uz(rf)=-uz(rf);
                end

                % --- bottom boundary: exits into substrate (index matched) ---
                ib=idx(hit_bot);
                if ~isempty(ib)
                    ang=acosd(min(max(-uz(ib),-1),1));
                    bins=min(max(ceil(ang),1),90);
                    w_R=w_R+accumarray(bins,1,[90 1]);
                    alive(ib)=false;
                end

                % --- Henyey-Greenstein scattering event ---
                isc=idx(scat);
                if ~isempty(isc)
                    z(isc)=z_new(scat);
                    n_sc=length(isc);

                    if abs(g)<1e-6
                        ct=2*rand(n_sc,1)-1;
                    elseif g>=1-1e-9
                        ct=ones(n_sc,1);            % pure forward: no deflection
                    else
                        s_hg=(1-g^2)./(1-g+2*g*rand(n_sc,1));
                        ct=(1+g^2-s_hg.^2)/(2*g);
                    end

                    st=sqrt(max(0,1-ct.^2));
                    phi=2*pi*rand(n_sc,1);
                    cp=cos(phi);
                    sp=sin(phi);

                    uxo=ux(isc); uyo=uy(isc); uzo=uz(isc);
                    denom=sqrt(max(1e-12,1-uzo.^2));
                    pole=abs(uzo)>0.99999;

                    nx=st.*(uxo.*uzo.*cp-uyo.*sp)./denom+uxo.*ct;
                    ny=st.*(uyo.*uzo.*cp+uxo.*sp)./denom+uyo.*ct;
                    nz=-st.*cp.*denom+uzo.*ct;

                    nx(pole)=st(pole).*cp(pole);
                    ny(pole)=st(pole).*sp(pole);
                    nz(pole)=sign(uzo(pole)).*ct(pole);

                    nrm=sqrt(nx.^2+ny.^2+nz.^2);
                    ux(isc)=nx./nrm;
                    uy(isc)=ny./nrm;
                    uz(isc)=nz./nrm;
                end

            end

            BSDF_SL(1:90,aoi,iS,ig)=w_T/N_photon;
            BSDF_SL(91:180,aoi,iS,ig)=flipud(w_R)/N_photon; % row 180 = theta_sub ~ 0

        end

        fprintf('S = %g, g = %g done (%g s)\n',S,g,toc);

    end
end

save('BSDF_SL_S_1to15by2_G_0.5to1by0.1_nSL_151.mat','BSDF_SL','S_list','g_list','n_SL','N_photon','-v7');

toc;

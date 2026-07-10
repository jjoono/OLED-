%% Reconstruction of the missing TMF_birefringence_whole.m
%
% The planar_Sweep22.m main script needs, for the two half-spaces seen from
% the dipole plane, the cavity reflection r (INCLUDING propagation through
% the partial first layer) for all u — radiative and evanescent — plus the
% amplitude transmission t (used only in the radiative window u < u_sub).
%
% The uploaded TMF_birefringence_whole_p/_s transfer-matrix routines are
% correct for propagating waves, but with MATLAB's principal sqrt their
% evanescent-layer branch makes |r| of an absorbing stack grow with u
% (unphysical). r is therefore computed here with a Parratt recursion on
% the physical branch Im(kz) >= 0 (evanescent waves decay), which coincides
% with the _p/_s result for u < 1 exactly. t is taken from the _p/_s
% routines, which the main script only ever evaluates in the propagating
% window. (Assumes the isotropic case ne == no, as in cps_oled.m.)

function TMF = TMF_birefringence_whole(no_bar, ne_bar, thickness, u, wavelength)

wavelength_num = length(wavelength);
u_num = length(u);

r_p = zeros(wavelength_num, u_num);
t_p = zeros(wavelength_num, u_num);
r_s = zeros(wavelength_num, u_num);
t_s = zeros(wavelength_num, u_num);

for i = 1:wavelength_num
    k0 = 2*pi/wavelength(i);
    ns = no_bar(i, :);                 % ne == no assumed
    kpar_p = real(ne_bar(i,1)) * u * k0;
    kpar_s = real(no_bar(i,1)) * u * k0;

    r_p(i,:) = parratt_r(ns, thickness, kpar_p, wavelength(i), 'p');
    r_s(i,:) = parratt_r(ns, thickness, kpar_s, wavelength(i), 's');

    % t from the original machinery (main script uses it only for u < u_sub)
    Tp = TMF_birefringence_whole_p(no_bar(i,:), ne_bar(i,:), thickness, kpar_p/k0, wavelength(i));
    Ts = TMF_birefringence_whole_s(no_bar(i,:), ne_bar(i,:), thickness, kpar_s/k0, wavelength(i));
    t_p(i,:) = Tp.t_p;
    t_s(i,:) = Ts.t_s;
end

TMF = struct('r_p', r_p, 't_p', t_p, 'r_s', r_s, 't_s', t_s);
end


function kzv = kz_phys(n, kpar, k0)
    v = sqrt((k0*n).^2 - kpar.^2);
    flip = imag(v) < 0;
    v(flip) = -v(flip);
    kzv = v;
end


function r = parratt_r(ns, thickness, kpar, lambda, pol)
% r seen from layer 1 (dipole medium), INCLUDING propagation through the
% partial first layer thickness(1), i.e. the cavity a-coefficient.
    k0 = 2*pi/lambda;
    N = length(ns);
    kzs = cell(1, N);
    for j = 1:N
        kzs{j} = kz_phys(ns(j), kpar, k0);
    end

    function rf = fresnel(a, b)
        if pol == 's'
            rf = (kzs{a} - kzs{b}) ./ (kzs{a} + kzs{b});
        else
            rf = (ns(b)^2*kzs{a} - ns(a)^2*kzs{b}) ./ (ns(b)^2*kzs{a} + ns(a)^2*kzs{b});
        end
    end

    r = fresnel(N-1, N);
    for j = N-2:-1:1
        ph = exp(2i * kzs{j+1} * thickness(j+1));
        rj = fresnel(j, j+1);
        r = (rj + r.*ph) ./ (1 + rj.*r.*ph);
    end
    % propagation from the dipole plane to the first interface
    r = r .* exp(2i * kzs{1} * thickness(1));
end

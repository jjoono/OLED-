function B = flat_bsdf(n_sub)
% FLAT_BSDF  BSDF (180 x 90, same layout as the LightTools files) of a flat
% substrate/air interface.  Column a+1 is the 1-degree bin centred on a degrees
% (the angular grid on which planar_inputs samples P_sub and R_OLED); the
% unpolarized Fresnel transmittance is averaged over the bin (sin-weighted) so
% that the critical angle does not snap to the grid.
B = zeros(180,90);
for a = 0:89
  t = linspace(max(0,a-0.5), a+0.5, 201); w = sind(t); if a == 0, w = w + eps; end
  s = n_sub*sind(t); T = zeros(size(t)); ok = s < 1;
  ci = cosd(t(ok)); ct = sqrt(1-s(ok).^2);
  rs = (n_sub*ci-ct)./(n_sub*ci+ct); rp = (ct-n_sub*ci)./(ct+n_sub*ci);
  T(ok) = 1 - (rs.^2+rp.^2)/2;
  Tm = sum(T.*w)/sum(w);
  if Tm > 0, B(min(90,floor(asind(min(1,n_sub*sind(a))))+1), a+1) = Tm; end   % only the row sum is used
  B(180-a, a+1) = 1 - Tm;                        % specular reflection, row 180 = 0 deg
end
end

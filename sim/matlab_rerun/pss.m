function p = pss(B)
% single-pass transmission averaged over a Lambertian distribution in the substrate
w = cosd(0:89).*sind(0:89); p = sum(sum(B(1:90,:),1).*w)/sum(w);
end

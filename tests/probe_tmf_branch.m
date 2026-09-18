load('nk_JH_total.mat'); wl=550; Al=material.l_Al_JO(151); n1=1.8;
u=[0.5 1.1 1.5];
% bare EML->Al interface (thickness [0]): TMF r_p vs analytic Fresnel
T=TMF_birefringence_whole([n1 Al],[n1 Al],[0],u,wl);
c1=sqrt(1-u.^2); c2=sqrt(1-(n1*u/Al).^2);
r_an=(Al*c1-n1*c2)./(Al*c1+n1*c2);
printf('u     TMF r_p                analytic r_p (principal sqrt)\n');
for k=1:3, printf('%.1f   %+.4f%+.4fi      %+.4f%+.4fi\n',u(k),real(T.r_p(k)),imag(T.r_p(k)),real(r_an(k)),imag(r_an(k))); end
% with a 25 nm EML gap in front of a 100 nm Al film then air
T2=TMF_birefringence_whole([n1 n1 Al 1],[n1 n1 Al 1],[25 100 0],u,wl);
xi=2*pi/wl*n1*c1;  r_an2=r_an.*exp(2i*xi*25);
printf('\nwith 25 nm gap: u   TMF r_p             analytic r*exp(2i xi z0)\n');
for k=1:3, printf('%.1f   %+.4f%+.4fi      %+.4f%+.4fi\n',u(k),real(T2.r_p(k)),imag(T2.r_p(k)),real(r_an2(k)),imag(r_an2(k))); end
% cos-theta branch used inside TMF for the Al layer
printf('\nTMF cos(Al) at u=1.1: %+.4f%+.4fi   analytic: %+.4f%+.4fi\n', real(T.cos_p(2,2)), imag(T.cos_p(2,2)), real(c2(2)), imag(c2(2)));
printf('TMF cos(EML) at u=1.1: %+.4f%+.4fi\n', real(T.cos_p(2,1)), imag(T.cos_p(2,1)));

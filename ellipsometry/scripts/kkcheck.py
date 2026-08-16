"""Diagnostic B: validate numerical KK in genosc_nk against the analytic
Lorentz oscillator, and test a singularity-subtracted KK for higher accuracy."""
import numpy as np, ellipsometry_fit as ef

E=ef._E_KK
A0,E0,C0=5.0,2.8,0.4   # Lorentz: eps = A/(E0^2 - E^2 - iCE)
den=(E0**2-E**2)**2+(C0*E)**2
eps2=A0*C0*E/den
eps1_ana=A0*(E0**2-E**2)/den          # eps1 - 1 analytic

# current method (matrix with zeroed diag +-1)
eps1_cur=ef._H_KK@eps2

# singularity-subtracted KK
dE=E[1]-E[0]
Ej=E[None,:]; Ei=E[:,None]
with np.errstate(divide='ignore',invalid='ignore'):
    M=Ej/(Ej**2-Ei**2)      # includes +-1 neighbours
    S=1.0/(Ej**2-Ei**2)
np.fill_diagonal(M,0.0); np.fill_diagonal(S,0.0)
Ssum=S.sum(axis=1)
a,b=E[0],E[-1]
Iana=(1.0/(2*E))*(np.log(np.abs((b-E)/(b+E)))-np.log(np.abs((a-E)/(a+E))))
d_eps2=np.gradient(eps2,E)
L=(eps2+E*d_eps2)/(2*E)     # j=i limit of regularized integrand
eps1_sub=(2/np.pi)*(dE*(M@eps2) - dE*E*eps2*Ssum + dE*L + E*eps2*Iana)

m=(E>1.0)&(E<5.0)
print('Lorentz KK validation (eps1-1), band 1-5 eV:')
print('  current  : max|err|=%.2e  mean|err|=%.2e'%(np.max(np.abs(eps1_cur[m]-eps1_ana[m])),np.mean(np.abs(eps1_cur[m]-eps1_ana[m]))))
print('  subtracted: max|err|=%.2e  mean|err|=%.2e'%(np.max(np.abs(eps1_sub[m]-eps1_ana[m])),np.mean(np.abs(eps1_sub[m]-eps1_ana[m]))))
print('  eps1 scale ~%.2f'%np.max(np.abs(eps1_ana)))

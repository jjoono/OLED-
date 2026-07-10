import numpy as np
import sys
sys.path.insert(0,'.')
from cps_ref import parratt, kz, K0, N_EML, N_AL, T_AL, T_ITO, D_EML, Z0
for u in (0.5, 0.95, 1.02, 1.05, 1.5):
    kpar=np.array([K0*N_EML*u])
    for pol in ('p','s'):
        r_b=parratt([N_EML,1.9,1.5],[T_ITO],kpar,pol)[0]
        r_t=parratt([N_EML,N_AL,1.0],[T_AL],kpar,pol)[0]
        a_b=r_b*np.exp(2j*kz(N_EML,kpar)[0]*(D_EML-Z0))
        a_t=r_t*np.exp(2j*kz(N_EML,kpar)[0]*Z0)
        print(f"u={u:5.2f} {pol}: a_b={a_b:.5f}  a_t={a_t:.5f}")

"""Corrected transfer matrix.

ellipsometry_fit._tmm builds the propagation matrix as diag(e^{+i b}, e^{-i b}).
For (E0+,E0-) = M (EN+,EN-) the forward wave must pick up e^{-i b} going back up,
so the two entries are swapped.  Consequence: with absorbing layers the backward
wave is amplified instead of attenuated and |r| can exceed 1 - a 200 nm Ag film
returns R = 1.033 instead of the bulk Fresnel value 0.9685.
"""
import numpy as np

def tmm(wl_nm, n_layers, d_list, theta0_deg, flip=False):
    th = np.deg2rad(theta0_deg); k0 = 2*np.pi/wl_nm
    s0 = n_layers[0]*np.sin(th)
    q = []
    for ni in n_layers:
        c = np.sqrt((ni**2 - s0**2).astype(complex))
        q.append(np.where(c.imag < 0, -c, c))
    N = len(wl_nm)
    def eye():
        M = np.zeros((N,2,2), complex); M[:,0,0]=1; M[:,1,1]=1; return M
    mm = lambda A,B: np.einsum('nij,njk->nik', A, B)
    def I_(r, t):
        M = np.zeros((N,2,2), complex)
        M[:,0,0]=1/t; M[:,0,1]=r/t; M[:,1,0]=r/t; M[:,1,1]=1/t; return M
    def Ip(ni,nj,qi,qj):
        den = nj**2*qi + ni**2*qj
        return I_((nj**2*qi - ni**2*qj)/den, 2*ni*nj*qi/den)
    def Is(qi,qj):
        den = qi+qj
        return I_((qi-qj)/den, 2*qi/den)
    def Pr(qj, d):
        b = k0*qj*d
        M = np.zeros((N,2,2), complex)
        s = 1.0 if flip else -1.0            # flip=True reproduces the old (buggy) sign
        M[:,0,0]=np.exp(s*1j*b); M[:,1,1]=np.exp(-s*1j*b); return M
    Mp, Ms = eye(), eye()
    for i, d in enumerate(d_list):
        ni,nj,qi,qj = n_layers[i],n_layers[i+1],q[i],q[i+1]
        Mp = mm(mm(Mp, Ip(ni,nj,qi,qj)), Pr(qj,d))
        Ms = mm(mm(Ms, Is(qi,qj)),       Pr(qj,d))
    ni,nj,qi,qj = n_layers[-2],n_layers[-1],q[-2],q[-1]
    Mp = mm(Mp, Ip(ni,nj,qi,qj)); Ms = mm(Ms, Is(qi,qj))
    return Mp[:,1,0]/Mp[:,0,0], Ms[:,1,0]/Ms[:,0,0]

if __name__ == '__main__':
    wl=np.array([600.0]); Ag=np.array([0.13+3.90j]); gl=np.array([1.52+0j]); air=np.array([1+0j])
    for ang in [0.0, 65.0]:
        rp,rs = tmm(wl,[air,Ag,gl],[200.0],ang)
        th=np.deg2rad(ang); s0=np.sin(th)
        q0=np.sqrt(1-s0**2+0j); q1=np.sqrt(Ag**2-s0**2+0j)
        bp=abs(((Ag**2*q0-q1)/(Ag**2*q0+q1))[0])**2; bs=abs(((q0-q1)/(q0+q1))[0])**2
        print('%3.0f deg  corrected Rp=%.4f Rs=%.4f   bulk Ag Rp=%.4f Rs=%.4f'
              %(ang,abs(rp[0])**2,abs(rs[0])**2,bp,bs))

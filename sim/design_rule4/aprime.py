"""Round-trip loss A' = 1 - <R_LED> of the simplified stack, seen from the substrate.
Plain TMM at 550 nm; validated against the Octave output of dr4e.m."""
import numpy as np
LAM = 550.0
N_SUB, N_ORG = 1.8, 1.8
N_TCO_DEFAULT = 1.9
N_AL = 0.958 + 6.687j
N_AG = 0.044 + 3.819j

def R_stack(nk, d, n_inc, th):
    """s- and p-reflectance of layers nk (list, finite thickness d in nm) on incident medium n_inc."""
    kx = n_inc*np.sin(th)
    ns = np.array([n_inc] + list(nk), dtype=complex)
    ds = np.array([0.0] + list(d))
    kz = np.sqrt(ns[:, None]**2 - kx[None, :]**2 + 0j)
    out = []
    for pol in ('s', 'p'):
        # interface reflection coefficients
        if pol == 's':
            eta = kz
        else:
            eta = ns[:, None]**2 / kz
        r = np.ones_like(kz[0], dtype=complex)*0
        # recurse from the last interface upward
        r = (eta[-2] - eta[-1])/(eta[-2] + eta[-1])
        for i in range(len(ns)-3, -1, -1):
            rij = (eta[i] - eta[i+1])/(eta[i] + eta[i+1])
            ph = np.exp(2j*2*np.pi/LAM*kz[i+1]*ds[i+1])
            r = (rij + r*ph)/(1 + rij*r*ph)
        out.append(np.abs(r)**2)
    return out

def aprime(d_tco, k_tco, metal, n_tco=N_TCO_DEFAULT, n_sub=N_SUB, d_org=(200., 20., 200.)):
    e = np.linspace(0, np.pi/2, 1801); th = 0.5*(e[1:] + e[:-1])
    nk = [n_tco + 1j*k_tco, N_ORG, N_ORG, N_ORG, metal]
    d = [d_tco, d_org[2], d_org[1], d_org[0], 100.]
    Rs, Rp = R_stack(nk, d, n_sub, th)
    w = np.cos(th)*np.sin(th)
    return 1 - np.sum(0.5*(Rs + Rp)*w)/np.sum(w)

if __name__ == '__main__':
    import sys
    ref = {'Al': N_AL, 'Ag': N_AG}
    print('validation against the Octave run (n_sub = 1.8, k_ITO = 0.02 and 0):')
    oct_ = {('Al',0.02,50):0.206390, ('Al',0.0,50):0.136912, ('Ag',0.02,50):0.094848, ('Ag',0.0,50):0.016619,
            ('Al',0.02,150):0.285614, ('Ag',0.02,150):0.189097}
    for (m, k, d), v in oct_.items():
        print(f'  {m} d={d:3.0f} k={k:.2f}:  python {aprime(d,k,ref[m]):.6f}   octave {v:.6f}')

"""Reconstruct n,k from the parameter set CompleteEASE converged to (ITO #1)
and compare with the chain-extracted n,k."""
import numpy as np, csv

HB = 0.6582119569          # eV*fs
HBS = 6.582119569e-16      # eV*s
EPS0 = 8.8541878128e-12

Eg_ = np.linspace(0.05, 30.0, 40000)   # wide grid for KK

def eps2_TL(E, A, Br, Eo, Eg):
    out = np.zeros_like(E)
    m = E > Eg
    Em = E[m]
    out[m] = (A*Eo*Br*(Em-Eg)**2)/(Em*((Em**2-Eo**2)**2 + Br**2*Em**2))
    return out

def eps2_G(E, A, Br, En):
    s = Br/(2*np.sqrt(np.log(2.0)))
    return A*(np.exp(-((E-En)/s)**2) - np.exp(-((E+En)/s)**2))

def kk(E, e2):
    """principal-value KK on a uniform grid, singularity-subtracted"""
    dE = E[1]-E[0]
    e1 = np.empty_like(E)
    for i in range(len(E)):
        num = E*e2
        den = E**2 - E[i]**2
        f = np.where(np.arange(len(E)) == i, 0.0, num/np.where(den == 0, 1, den))
        # subtract the singular part analytically
        a, b = E[0], E[-1]
        Ei = E[i]
        g = E[i]*e2[i]
        f = f - np.where(np.arange(len(E)) == i, 0.0, g/den)
        Ians = (1.0/(2*Ei))*(np.log(abs((b-Ei)/(b+Ei))) - np.log(abs((a-Ei)/(a+Ei))))
        e1[i] = (2/np.pi)*(np.sum(f)*dE + g*Ians)
    return e1

def model(E, einf, rho, tau, TL, G):
    e2 = eps2_TL(Eg_, *TL) + eps2_G(Eg_, *G)
    e1 = kk(Eg_, e2)
    e1i = np.interp(E, Eg_, e1); e2i = np.interp(E, Eg_, e2)
    Br_d = HB/tau
    A_d = HBS**2/(EPS0*(rho/100.0)*(tau*1e-15))     # eV^2
    e1d = -A_d/(E**2 + Br_d**2)
    e2d = A_d*Br_d/(E*(E**2 + Br_d**2))
    eps = (einf + e1i + e1d) + 1j*(e2i + e2d)
    N = np.sqrt(eps)
    return N.real, N.imag, (e2i+e2d)

wl = np.array([300,350,400,450,500,550,600,650,691,700,750,800,900,1000,1078,1200,1500,1800],float)
E = 1239.841984/wl

CE = dict(einf=2.119, rho=5.3893e-4, tau=5.423,
          TL=(134.7806, 6.171, 4.256, 2.953), G=(-0.093778, 1.4943, 1.794))
MY = dict(einf=2.5698279099768726, rho=1.0644047289948939e-3, tau=6.582119569,
          TL=(165.546636633031, 1.0550595369907967, 4.000000000000001, 3.3367948707895163),
          G=(0.07201885654941534, 0.8249493158431357, 1.911146995799902))

nc, kc, e2c = model(E, **CE)
nm, km, e2m = model(E, **MY)

ch = [x for x in csv.reader(open('ITO_IZO_final_nk.csv')) if x and not x[0].startswith('#')]
D = np.array([[float(v) for v in x[:3]] for x in ch[1:] if x[0].replace('.','').isdigit()])

print('  wl      CE-fit n    CE-fit k    eps2_CE  |  mine n    mine k   | chain n  chain k')
for i, w in enumerate(wl):
    if w <= 1078:
        j = np.argmin(abs(D[:,0]-w)); cn, ck = D[j,1], D[j,2]
        s = '%8.3f %8.4f' % (cn, ck)
    else:
        s = '   --       --   (no data)'
    flag = '  <== k<0' if kc[i] < 0 else ''
    print('%6.0f  %8.3f  %+9.4f  %+9.4f  | %7.3f %8.4f  |%s%s'
          % (w, nc[i], kc[i], e2c[i], nm[i], km[i], s, flag))

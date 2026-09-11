"""Cross-check the plasma-energy density estimate against an explicit
effective-medium fit.

The number I reported came from the Drude term only:  wp^2 = N e^2/(eps0 m*),
so (wp/wp_ref)^2 is the free-electron density relative to the reference film.
That is exact for a simple dilution of carriers, but a real porous/island film
mixes DIELECTRIC FUNCTIONS, which is not linear.  So here the n,k of the most
bulk-like film (2-8, 12 nm on HATCN) is taken as the Ag reference and each
sample's n,k is refitted as reference-Ag + void, once with Bruggeman (a
percolating mixture) and once with Maxwell-Garnett (isolated inclusions).
"""
import numpy as np, json
from scipy.optimize import minimize_scalar
import ag_load as L, ag_fit as AF

R = json.load(open('ag_polish_result.json'))
HB = 6.582119569e-16; EPS0 = 8.8541878128e-12
wl = np.linspace(400, 1080, 200)
REF = '2-8'
N_ref = AF.agN(np.array(R[REF]['p']), wl)
e_ref = N_ref**2

def wp_of(p):
    return np.sqrt(HB**2/(EPS0*(p[3]/100.0)*(p[4]*1e-15)))
WP_REF = wp_of(np.array(R[REF]['p']))

import ellipsometry_fit as ef
def brugg(ea, eb, f):
    """library version - it picks the physical root (Im N >= 0); the naive
    principal branch of the square root is wrong for metals"""
    Na = np.sqrt(ea); Na = np.where(Na.imag<0,-Na,Na)
    Nb = np.sqrt(eb); Nb = np.where(Nb.imag<0,-Nb,Nb)
    return ef.bruggeman_ema(Na, Nb, f)**2

def mg(ei, f):                       # inclusions ei in a void host
    beta = (ei-1.0)/(ei+2.0)
    e = (1 + 2*f*beta)/(1 - f*beta)
    N = np.sqrt(e); N = np.where(N.imag<0,-N,N)
    return N**2

def best_f(e_meas, mix):
    def cost(f):
        e = mix(np.clip(f, 1e-4, 1.0))
        N = np.sqrt(e); N = np.where(N.imag<0,-N,N)
        Nm = np.sqrt(e_meas); Nm = np.where(Nm.imag<0,-Nm,Nm)
        return np.mean((N.real-Nm.real)**2 + (N.imag-Nm.imag)**2)
    r = minimize_scalar(cost, bounds=(0.02, 1.0), method='bounded')
    return r.x, np.sqrt(r.fun)

print('%-5s %-6s %4s | %8s | %8s %8s | %8s %8s'%(
      'id','seed','nom','f (wp^2)','f Brugg','rms n,k','f MaxGar','rms n,k'))
for s in L.ORDER:
    if s not in R: continue
    p = np.array(R[s]['p'])
    f_wp = (wp_of(p)/WP_REF)**2
    e = AF.agN(p, wl)**2
    fb, rb = best_f(e, lambda f: brugg(e_ref, np.ones_like(e_ref), f))
    fm, rm = best_f(e, lambda f: mg(e_ref, f))
    print('%-5s %-6s %4d | %8.2f | %8.2f %8.4f | %8.2f %8.4f'%(
          s, R[s]['seed'], R[s]['ag_nom'], f_wp, fb, rb, fm, rm))

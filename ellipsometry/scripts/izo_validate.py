"""Validate forward machinery: reconstruct CompleteEASE's N2_IZO1 result from
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

its .mod (B-spline nodes + EInf + IR Amp + thicknesses) and compare to data."""
import numpy as np, ellipsometry_fit as ef
from scipy.interpolate import PchipInterpolator, BSpline
d=np.load(str(OUT_DIR / 'izo_data.npz'))
wl=d['n2wl']; psi_m=d['n2psi']; del_m=d['n2del']
angles=np.array([45.0,55.0,65.0,75.0])+0.0224
E1,V1=d['ceE1'],d['ceV1']; einf=float(d['ce_einf1']); irA=float(d['ce_ir1'])
D_IZO=189.2952; D_R=1.6013; D_OX=3.00   # nm (from mod, Angstrom/10)

# local KK on wide grid (Drude tail needs low E)
Eg_=np.linspace(0.05,9.5,2400); dE=Eg_[1]-Eg_[0]
Ej=Eg_[None,:]; Ei=Eg_[:,None]
with np.errstate(divide='ignore',invalid='ignore'):
    M=Ej/(Ej**2-Ei**2); S=1.0/(Ej**2-Ei**2)
np.fill_diagonal(M,0.0); np.fill_diagonal(S,0.0)
Ssum=S.sum(1)
a_,b_=Eg_[0],Eg_[-1]
with np.errstate(divide='ignore'):
    Iana=(1/(2*Eg_))*(np.log(np.abs((b_-Eg_)/(b_+Eg_)))-np.log(np.abs((a_-Eg_)/(a_+Eg_))))
Iana[~np.isfinite(Iana)]=0.0
def kk(e2):
    L=(e2+Eg_*np.gradient(e2,Eg_))/(2*Eg_)
    return (2/np.pi)*(dE*(M@e2)-dE*Eg_*e2*Ssum+dE*L+Eg_*e2*Iana)

def eps_from_nodes(mode):
    if mode=='pchip':
        f=PchipInterpolator(E1,V1,extrapolate=False)
        e2=f(Eg_); e2[~np.isfinite(e2)]=0.0
    else:  # bspline coefficients on given knots
        t=np.concatenate(([E1[0]]*3,E1,[E1[-1]]*3))
        spl=BSpline(t,np.concatenate((V1,[V1[-1]]*2)),3,extrapolate=False)
        e2=spl(Eg_); e2[~np.isfinite(e2)]=0.0
    return np.maximum(e2,0.0)

E_wl=1239.84193/wl
N_si=np.sqrt((np.interp(wl,d['si_wl'],d['si_e1'])+1j*np.interp(wl,d['si_wl'],d['si_e2'])).astype(complex))
N_ox=np.sqrt((np.interp(wl,d['ox_wl'],d['ox_e1'])+1j*np.interp(wl,d['ox_wl'],d['ox_e2'])).astype(complex))
n_air=np.ones(len(wl),dtype=complex)

def model(e2g,einf,irA,irsign,d_izo,d_r):
    e1g=einf+irsign*irA/Eg_**2+kk(e2g)
    e1=np.interp(E_wl,Eg_,e1g); e2=np.interp(E_wl,Eg_,e2g)
    Nf=np.sqrt((e1+1j*e2).astype(complex))
    Nf=np.where(Nf.imag<0,np.conj(Nf),Nf)
    layers=[n_air, ef.bruggeman_ema(Nf,n_air,0.5), Nf, N_ox, N_si]
    dl=[d_r,d_izo,D_OX]
    P=[];L=[]
    for ang in angles:
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        P.append(np.degrees(np.arctan(np.abs(rho)))); L.append(np.degrees(np.angle(rho)))
    return np.array(P).T,np.array(L).T

for mode in ['pchip','bspl']:
    e2g=eps_from_nodes(mode)
    for irsign in [-1,+1]:
        pc,dc=model(e2g,einf,irA,irsign,D_IZO,D_R)
        dd=(dc-del_m+180)%360-180
        print('%5s irsign=%+d: RMSE psi=%.3f del=%.3f'%(mode,irsign,
              np.sqrt(np.mean((pc-psi_m)**2)),np.sqrt(np.mean(dd**2))))
# NCS-space too for best combo report later

"""IZO multi-start fit: analytic Drude + UV-interband spline(e2,KK) + eps_inf."""
import sys, numpy as np, ellipsometry_fit as ef
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

d=np.load(str(OUT_DIR / 'izo_data.npz'))

which=sys.argv[1] if len(sys.argv)>1 else 'N2'
if which=='N2':
    wl=d['n2wl']; psi_m=d['n2psi']; del_m=d['n2del']; dep=d['n2dep']
    angles=np.array([45.0,55.0,65.0,75.0]); ce_d=189.3
else:
    wl=d['o2wl']; psi_m=d['o2psi']; del_m=d['o2del']; dep=d['o2dep']
    angles=np.array([45.0,50.0,55.0,60.0,65.0]); ce_d=179.1

E_wl=1239.84193/wl
N_si=np.sqrt((np.interp(wl,d['si_wl'],d['si_e1'])+1j*np.interp(wl,d['si_wl'],d['si_e2'])).astype(complex))
N_ox=np.sqrt((np.interp(wl,d['ox_wl'],d['ox_e1'])+1j*np.interp(wl,d['ox_wl'],d['ox_e2'])).astype(complex))
n_air=np.ones(len(wl),dtype=complex)
W=1.0/(1.0+(np.abs(dep)/3.0)**2)

Eg=np.linspace(0.4,9.0,1600); dE=Eg[1]-Eg[0]
Ej=Eg[None,:]; Ei=Eg[:,None]
with np.errstate(divide='ignore',invalid='ignore'):
    M=Ej/(Ej**2-Ei**2); S=1.0/(Ej**2-Ei**2)
np.fill_diagonal(M,0); np.fill_diagonal(S,0)
Ssum=S.sum(1); a_,b_=Eg[0],Eg[-1]
with np.errstate(divide='ignore'):
    Iana=(1/(2*Eg))*(np.log(np.abs((b_-Eg)/(b_+Eg)))-np.log(np.abs((a_-Eg)/(a_+Eg))))
Iana[~np.isfinite(Iana)]=0
def kk(e2):
    L=(e2+Eg*np.gradient(e2,Eg))/(2*Eg)
    return (2/np.pi)*(dE*(M@e2)-dE*Eg*e2*Ssum+dE*L+Eg*e2*Iana)

En=np.array([2.6,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.6,4.9,5.2,5.6,6.2,7.0,8.0])
NN=len(En)

def eps_layer(p):
    C0,Ad,gd=p[0],p[1],p[2]
    nodes=np.maximum(p[5:5+NN],0.0)
    e2uv=np.maximum(PchipInterpolator(np.concatenate(([Eg[0],2.2],En)),
                    np.concatenate(([0,0],nodes)))(Eg),0.0)
    e1=C0+kk(e2uv)-Ad/(Eg**2+gd**2)
    e2=e2uv+Ad*gd/(Eg*(Eg**2+gd**2))
    e1w=np.interp(E_wl,Eg,e1); e2w=np.interp(E_wl,Eg,e2)
    Nf=np.sqrt((e1w+1j*np.maximum(e2w,0)).astype(complex))
    return np.where(Nf.imag<0,np.conj(Nf),Nf)

def forward(p):
    Nf=eps_layer(p)
    layers=[n_air, ef.bruggeman_ema(Nf,n_air,0.5), Nf, N_ox, N_si]
    dl=[p[4],p[3],3.0]
    P=[];L=[]
    for ang in angles:
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        P.append(np.degrees(np.arctan(np.abs(rho)))); L.append(np.degrees(np.angle(rho)))
    return np.array(P).T,np.array(L).T

def resid(p):
    P,L=forward(p)
    dd=(L-del_m+180)%360-180
    nodes=p[5:5+NN]
    reg=0.3*np.diff(nodes,2)          # curvature regularization
    return np.concatenate([(W*(P-psi_m)).ravel(),(W*dd/3.0).ravel(),reg])

lo=np.array([1.5, 0.0, 0.05, 140.0, 0.0]+[0]*NN)
hi=np.array([5.0, 30.0, 2.0, 230.0, 6.0]+[25]*NN)
nodes0=[0.3,0.5,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.5,2.6,2.7,2.8,2.8]
best=None
for dsd in [155,165,175,185,195,205]:
    p0=np.array([3.2, 2.0, 0.3, dsd, 1.5]+nodes0)
    r=least_squares(resid,p0,bounds=(lo,hi),xtol=1e-10,ftol=1e-10,max_nfev=250)
    P,L=forward(r.x); dd=(L-del_m+180)%360-180
    rms=np.sqrt(np.mean((P-psi_m)**2+dd**2))
    print('seed d=%d -> d=%.1f cost=%.1f rms=%.2f'%(dsd,r.x[3],r.cost,rms))
    if best is None or r.cost<best.cost: best=r
r=least_squares(resid,best.x,bounds=(lo,hi),xtol=1e-12,ftol=1e-12,max_nfev=500)
p=r.x
P,L=forward(p); dd=(L-del_m+180)%360-180
print('\n%s BEST: d=%.2fnm (CE %.1f) rough=%.2f C0=%.3f Drude A=%.3f g=%.3f'%(
      which,p[3],ce_d,p[4],p[0],p[1],p[2]))
print('RMSE psi=%.3f del=%.3f'%(np.sqrt(np.mean((P-psi_m)**2)),np.sqrt(np.mean(dd**2))))
Nf=eps_layer(p)
for tw in [300,350,400,550,800,1200,1600]:
    i=int(np.argmin(abs(wl-tw)))
    print('  %4dnm n=%.3f k=%.4f'%(wl[i],Nf[i].real,Nf[i].imag))
np.savez(str(OUT_DIR / 'izo_fit_%s.npz')%which,p=p,wl=wl,psi_m=psi_m,del_m=del_m,
         psi_c=P,del_c=L,n=Nf.real,k=Nf.imag,W=W,angles=angles,En=En)
print('saved izo_fit_%s.npz'%which)

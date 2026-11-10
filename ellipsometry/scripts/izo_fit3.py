"""Staged IZO fit: S1 transparent-window Cauchy -> S2 full-range optical
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

(d frozen) -> S3 all free.  Drude(analytic) + UV-interband spline(e2,KK) + C0."""
import sys, numpy as np, ellipsometry_fit as ef
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares
d=np.load(str(OUT_DIR / 'izo_data.npz'))

which=sys.argv[1] if len(sys.argv)>1 else 'N2'
if which=='N2':
    wl=d['n2wl']; psi_m=d['n2psi']; del_m=d['n2del']; dep=d['n2dep']
    angles=np.array([45.0,55.0,65.0,75.0])
else:
    wl=d['o2wl']; psi_m=d['o2psi']; del_m=d['o2del']; dep=d['o2dep']
    angles=np.array([45.0,50.0,55.0,60.0,65.0])

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

En=np.array([2.6,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.6,4.9,5.2,5.6,6.2,7.0,8.0]); NN=len(En)

def eps_layer(C0,Ad,gd,nodes):
    nodes=np.maximum(nodes,0.0)
    e2uv=np.maximum(PchipInterpolator(np.concatenate(([Eg[0],2.2],En)),
                    np.concatenate(([0,0],nodes)))(Eg),0.0)
    e1=C0+kk(e2uv)-Ad/(Eg**2+gd**2)
    e2=e2uv+Ad*gd/(Eg*(Eg**2+gd**2))
    e1w=np.interp(E_wl,Eg,e1); e2w=np.interp(E_wl,Eg,e2)
    Nf=np.sqrt((e1w+1j*np.maximum(e2w,0)).astype(complex))
    return np.where(Nf.imag<0,np.conj(Nf),Nf)

def forward(Nf,dz,dr):
    layers=[n_air, ef.bruggeman_ema(Nf,n_air,0.5), Nf, N_ox, N_si]
    dl=[dr,dz,3.0]
    P=[];L=[]
    for ang in angles:
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        P.append(np.degrees(np.arctan(np.abs(rho)))); L.append(np.degrees(np.angle(rho)))
    return np.array(P).T,np.array(L).T

# ---------- Stage 1: transparent window Cauchy ----------
mask=(wl>=600)&(wl<=1100)
def s1res(x):
    dz,dr,A,B=x
    Nf=((A+B/wl**2)+0j).astype(complex)
    P,L=forward(Nf,dz,dr)
    dd=(L-del_m+180)%360-180
    return np.concatenate([(W*(P-psi_m))[mask].ravel(),(W*dd/3.0)[mask].ravel()])
best=None
for d0 in range(140,231,10):
    r=least_squares(s1res,[d0,1.5,2.0,2e4],bounds=([120,0,1.5,0],[260,6,2.7,2e5]))
    if best is None or r.cost<best.cost: best=r
d1,dr1,A1,B1=best.x
n800=A1+B1/800**2
print('[S1] d=%.2f rough=%.2f n550=%.3f n800=%.3f'%(d1,dr1,A1+B1/550**2,n800))

# ---------- Stage 2: full range, d/dr frozen ----------
nodes0=np.array([0.3,0.5,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.5,2.6,2.7,2.8,2.8])
C0_0=n800**2+0.3
def s2res(x):
    C0,Ad,gd=x[0],x[1],x[2]; nodes=x[3:3+NN]
    Nf=eps_layer(C0,Ad,gd,nodes)
    P,L=forward(Nf,d1,dr1)
    dd=(L-del_m+180)%360-180
    reg=0.3*np.diff(nodes,2)
    return np.concatenate([(W*(P-psi_m)).ravel(),(W*dd/3.0).ravel(),reg])
lo2=np.array([1.5,0.0,0.05]+[0]*NN); hi2=np.array([6.0,30.0,2.0]+[25]*NN)
b2=None
for Ad0,gd0 in [(2.0,0.4),(0.5,0.2),(4.0,0.8),(0.05,0.3)]:
    r=least_squares(s2res,np.concatenate(([C0_0,Ad0,gd0],nodes0)),bounds=(lo2,hi2),max_nfev=200)
    if b2 is None or r.cost<b2.cost: b2=r
x2=b2.x
print('[S2] C0=%.3f Ad=%.3f gd=%.3f  cost=%.0f'%(x2[0],x2[1],x2[2],b2.cost))

# ---------- Stage 3: all free ----------
def s3res(x):
    dz,dr,C0,Ad,gd=x[:5]; nodes=x[5:5+NN]
    Nf=eps_layer(C0,Ad,gd,nodes)
    P,L=forward(Nf,dz,dr)
    dd=(L-del_m+180)%360-180
    reg=0.3*np.diff(nodes,2)
    return np.concatenate([(W*(P-psi_m)).ravel(),(W*dd/3.0).ravel(),reg])
lo3=np.concatenate(([d1-15,0,1.5,0,0.05],[0]*NN))
hi3=np.concatenate(([d1+15,6,6.0,30,2.0],[25]*NN))
x30=np.concatenate(([d1,dr1],x2))
r3=least_squares(s3res,x30,bounds=(lo3,hi3),xtol=1e-12,ftol=1e-12,max_nfev=600)
x=r3.x
dz,dr,C0,Ad,gd=x[:5]; nodes=x[5:5+NN]
Nf=eps_layer(C0,Ad,gd,nodes)
P,L=forward(Nf,dz,dr)
dd=(L-del_m+180)%360-180
print('\n[%s FINAL] d=%.2fnm rough=%.2fnm C0=%.3f Drude A=%.3f g=%.3f'%(which,dz,dr,C0,Ad,gd))
print('RMSE psi=%.3f del=%.3f'%(np.sqrt(np.mean((P-psi_m)**2)),np.sqrt(np.mean(dd**2))))
for tw in [300,350,400,550,800,1200,1600]:
    i=int(np.argmin(abs(wl-tw)))
    print('  %4dnm n=%.3f k=%.4f'%(wl[i],Nf[i].real,Nf[i].imag))
np.savez(str(OUT_DIR / 'izo_fit_%s.npz')%which,x=x,wl=wl,psi_m=psi_m,del_m=del_m,
         psi_c=P,del_c=L,n=Nf.real,k=Nf.imag,W=W,angles=angles,En=En,nodes=nodes,
         d=dz,dr=dr,C0=C0,Ad=Ad,gd=gd)
print('saved izo_fit_%s.npz'%which)

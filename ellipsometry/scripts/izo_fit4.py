"""IZO fit v4: NCS-space residuals (bounded; immune to Psi->90 / Delta-wrap
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

pathologies) + full-range e2 spline with KK + C0 + angle offset. Staged."""
import sys, numpy as np, ellipsometry_fit as ef
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares
d=np.load(str(OUT_DIR / 'izo_data.npz'))
which=sys.argv[1] if len(sys.argv)>1 else 'O2'
if which=='N2':
    wl=d['n2wl']; psi_m=d['n2psi']; del_m=d['n2del']; dep=d['n2dep']
    ang0=np.array([45.0,55.0,65.0,75.0])
else:
    wl=d['o2wl']; psi_m=d['o2psi']; del_m=d['o2del']; dep=d['o2dep']
    ang0=np.array([45.0,50.0,55.0,60.0,65.0])
NA=len(ang0)
E_wl=1239.84193/wl
N_si=np.sqrt((np.interp(wl,d['si_wl'],d['si_e1'])+1j*np.interp(wl,d['si_wl'],d['si_e2'])).astype(complex))
N_ox=np.sqrt((np.interp(wl,d['ox_wl'],d['ox_e1'])+1j*np.interp(wl,d['ox_wl'],d['ox_e2'])).astype(complex))
n_air=np.ones(len(wl),dtype=complex)
W=1.0/(1.0+(np.abs(dep)/3.0)**2)

pr=np.deg2rad(psi_m); dr_=np.deg2rad(del_m)
Nm=np.cos(2*pr); Cm=np.sin(2*pr)*np.cos(dr_); Sm=np.sin(2*pr)*np.sin(dr_)

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

En=np.array([0.55,0.7,0.85,1.0,1.15,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.85,3.0,3.15,
             3.3,3.45,3.6,3.75,3.9,4.05,4.2,4.4,4.6,4.9,5.2,5.6,6.2,7.0,8.2]); NN=len(En)

def eps_layer(C0,nodes):
    nodes=np.maximum(nodes,0.0)
    e2=np.maximum(PchipInterpolator(En,nodes,extrapolate=False)(Eg),0.0)
    e2[~np.isfinite(e2)]=0.0
    e1=C0+kk(e2)
    e1w=np.interp(E_wl,Eg,e1); e2w=np.interp(E_wl,Eg,e2)
    Nf=np.sqrt((e1w+1j*np.maximum(e2w,0)).astype(complex))
    return np.where(Nf.imag<0,np.conj(Nf),Nf)

def ncs_model(Nf,dz,drg,dth):
    layers=[n_air, ef.bruggeman_ema(Nf,n_air,0.5), Nf, N_ox, N_si]
    dl=[drg,dz,3.0]
    out=[]
    for ang in ang0+dth:
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        t=np.abs(rho); psi2=2*np.arctan(t)
        Np_=np.cos(psi2); Cp=np.sin(psi2)*np.cos(np.angle(rho)); Sp=np.sin(psi2)*np.sin(np.angle(rho))
        out.append((Np_,Cp,Sp))
    return out

def resid_full(x,freeze_geo=None):
    if freeze_geo is None:
        dz,drg,dth,C0=x[0],x[1],x[2],x[3]; nodes=x[4:4+NN]
    else:
        dz,drg,dth=freeze_geo; C0=x[0]; nodes=x[1:1+NN]
    Nf=eps_layer(C0,nodes)
    mm=ncs_model(Nf,dz,drg,dth)
    res=[]
    for j,(Np_,Cp,Sp) in enumerate(mm):
        res+= [W[:,j]*(Np_-Nm[:,j]),W[:,j]*(Cp-Cm[:,j]),W[:,j]*(Sp-Sm[:,j])]
    reg=0.10*np.diff(np.maximum(nodes,0),2)
    return np.concatenate(res+[reg])

# ---- S1: window Cauchy in NCS ----
mask=(wl>=600)&(wl<=1050)
def s1(x):
    dz,drg,dth,A,B=x
    Nf=((A+B/wl**2)+0j).astype(complex)
    mm=ncs_model(Nf,dz,drg,dth)
    res=[]
    for j,(Np_,Cp,Sp) in enumerate(mm):
        res+=[ (W[:,j]*(Np_-Nm[:,j]))[mask],(W[:,j]*(Cp-Cm[:,j]))[mask],(W[:,j]*(Sp-Sm[:,j]))[mask] ]
    return np.concatenate(res)
b1=None
for d0 in range(150,221,10):
    r=least_squares(s1,[d0,1.0,0.0,2.0,2e4],bounds=([130,0,-0.5,1.5,0],[250,8,0.5,2.7,2e5]))
    if b1 is None or r.cost<b1.cost: b1=r
dz1,dr1,dth1,A1,B1=b1.x
print('[S1] d=%.2f rough=%.2f dth=%.3f n550=%.3f n800=%.3f cost=%.4f'%(dz1,dr1,dth1,A1+B1/550**2,A1+B1/800**2,b1.cost))

# ---- S2: nodes with geometry frozen ----
n800=A1+B1/800**2
nodes0=np.full(NN,0.05); nodes0[En>3.2]=1.0; nodes0[En>4.0]=2.0
x0=np.concatenate(([n800**2],nodes0))
lo=np.concatenate(([1.5],np.zeros(NN))); hi=np.concatenate(([6.0],np.full(NN,25)))
r2=least_squares(lambda x: resid_full(x,freeze_geo=(dz1,dr1,dth1)),x0,bounds=(lo,hi),max_nfev=300)
print('[S2] C0=%.3f cost=%.4f'%(r2.x[0],r2.cost))

# ---- S3: all free ----
x0=np.concatenate(([dz1,dr1,dth1],r2.x))
lo3=np.concatenate(([dz1-12,0,-0.5,1.5],np.zeros(NN)))
hi3=np.concatenate(([dz1+12,8,0.5,6.0],np.full(NN,25)))
r3=least_squares(resid_full,x0,bounds=(lo3,hi3),xtol=1e-12,ftol=1e-12,max_nfev=800)
x=r3.x
dz,drg,dth,C0=x[0],x[1],x[2],x[3]; nodes=x[4:4+NN]
Nf=eps_layer(C0,nodes)
mm=ncs_model(Nf,dz,drg,dth)
# report in psi/del for familiarity
P=[];L=[]
for (Np_,Cp,Sp) in mm:
    psi=0.5*np.degrees(np.arccos(np.clip(Np_,-1,1)))
    dl_=np.degrees(np.arctan2(Sp,Cp))
    P.append(psi);L.append(dl_)
P=np.array(P).T;L=np.array(L).T
dd=(L-del_m+180)%360-180
ncs_rms=np.sqrt(2*r3.cost/ (3*NA*len(wl)))
print('\n[%s v4 FINAL] d=%.2f rough=%.2f dtheta=%.3f C0=%.3f  NCS-rms=%.5f'%(which,dz,drg,dth,C0,ncs_rms))
print('psi/del RMSE: %.3f / %.3f'%(np.sqrt(np.mean((P-psi_m)**2)),np.sqrt(np.mean(dd**2))))
for tw in [300,350,400,550,800,1200,1600]:
    i=int(np.argmin(abs(wl-tw)))
    print('  %4dnm n=%.3f k=%.4f'%(wl[i],Nf[i].real,Nf[i].imag))
np.savez(str(OUT_DIR / 'izo_fit4_%s.npz')%which,x=x,wl=wl,psi_m=psi_m,del_m=del_m,
         psi_c=P,del_c=L,n=Nf.real,k=Nf.imag,W=W,angles=ang0,En=En,nodes=nodes,d=dz,dr=drg,dth=dth,C0=C0)
print('saved izo_fit4_%s.npz'%which)

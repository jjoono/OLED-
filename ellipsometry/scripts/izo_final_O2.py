"""O2_IZO3 final: (1) parametric KK model (C0+Drude+TL) fitted to chain nk,
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

(2) joint NCS polish on data, (3) full report + save."""
import numpy as np, ellipsometry_fit as ef
from scipy.optimize import least_squares
d=np.load(str(OUT_DIR / 'izo_data.npz'))
ch=np.load(str(OUT_DIR / 'izo_chain_O2.npz'))
wl=d['o2wl']; psi_m=d['o2psi']; del_m=d['o2del']; dep=d['o2dep']
ang0=np.array([45.0,50.0,55.0,60.0,65.0]); NA=5
E_wl=1239.84193/wl
N_si=np.sqrt((np.interp(wl,d['si_wl'],d['si_e1'])+1j*np.interp(wl,d['si_wl'],d['si_e2'])).astype(complex))
N_ox=np.sqrt((np.interp(wl,d['ox_wl'],d['ox_e1'])+1j*np.interp(wl,d['ox_wl'],d['ox_e2'])).astype(complex))
n_air=np.ones(len(wl),dtype=complex)
W=1.0/(1.0+(np.abs(dep)/3.0)**2)
pr=np.deg2rad(psi_m); dr_=np.deg2rad(del_m)
Nm=np.cos(2*pr); Cm=np.sin(2*pr)*np.cos(dr_); Sm=np.sin(2*pr)*np.sin(dr_)

def tl_eps2(E,A,E0,C,Egp):
    den=(E**2-E0**2)**2+C**2*E**2
    return np.where(E>Egp, A*E0*C*(E-Egp)**2/(den*np.maximum(E,1e-9)), 0.0)

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

def eps_param(p,Earr):
    C0,Ad,gd,At,E0t,Ct,Egt,Ag,Eg_c,Br=p
    e2i=tl_eps2(Eg,At,E0t,Ct,Egt)+Ag*(np.exp(-((Eg-Eg_c)/Br)**2)-np.exp(-((Eg+Eg_c)/Br)**2))
    e2i=np.maximum(e2i,0)
    e1g=C0+kk(e2i)-Ad/(Eg**2+gd**2)
    e2g=e2i+Ad*gd/(Eg*(Eg**2+gd**2))
    return np.interp(Earr,Eg,e1g), np.interp(Earr,Eg,np.maximum(e2g,0))

# ---- (1) fit parametric to chain nk ----
wlc=ch['wl']; nc=ch['n']; kc=ch['k']; resc=ch['res']
Ec=1239.84193/wlc
e1c=nc**2-kc**2; e2c=2*nc*kc
wq=1.0/(0.05+resc)          # trust low-res chain points
def r1(p):
    e1,e2=eps_param(p,Ec)
    return np.concatenate([wq*(e1-e1c),wq*(e2-e2c)*2.0])
p0=[4.4,1.73,0.03, 8.0,4.6,1.0,3.6, 0.3,4.2,0.4]
lo=[2.0,0.0,0.004, 0.0,3.6,0.1,2.8, 0.0,3.4,0.05]
hi=[6.5,6.0,0.5, 200.0,7.5,3.0,4.4, 5.0,5.2,1.2]
b=None
for gd0 in [0.01,0.03,0.1]:
    p0[2]=gd0
    r=least_squares(r1,p0,bounds=(lo,hi),max_nfev=2000)
    if b is None or r.cost<b.cost: b=r
p1=b.x
print('[param->chain] C0=%.3f Drude A=%.3f g=%.4f | TL A=%.1f E0=%.2f C=%.2f Eg=%.2f | G A=%.2f Ec=%.2f Br=%.2f'%tuple(p1))
e1f,e2f=eps_param(p1,Ec)
print('  chain-match rms e1=%.3f e2=%.3f'%(np.sqrt(np.mean((e1f-e1c)**2)),np.sqrt(np.mean((e2f-e2c)**2))))

# ---- (2) joint NCS polish ----
def ncs_forward(p,dz,drg,dth):
    e1,e2=eps_param(p,E_wl)
    Nf=np.sqrt((e1+1j*e2).astype(complex))
    Nf=np.where(Nf.imag<0,np.conj(Nf),Nf)
    layers=[n_air, ef.bruggeman_ema(Nf,n_air,0.5), Nf, N_ox, N_si]
    dl=[drg,dz,3.0]
    out=[]
    for ang in ang0+dth:
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        t=np.abs(rho); p2=2*np.arctan(t); an=np.angle(rho)
        out.append((np.cos(p2),np.sin(p2)*np.cos(an),np.sin(p2)*np.sin(an)))
    return out,Nf

def r2(x):
    dz,drg,dth=x[0],x[1],x[2]; p=x[3:]
    mm,_=ncs_forward(p,dz,drg,dth)
    res=[]
    for j,(Np_,Cp,Sp) in enumerate(mm):
        res+=[W[:,j]*(Np_-Nm[:,j]),W[:,j]*(Cp-Cm[:,j]),W[:,j]*(Sp-Sm[:,j])]
    return np.concatenate(res)
x0=np.concatenate(([float(ch['d']),1.0,0.0],p1))
lo2=np.concatenate(([170,0,-0.3],lo)); hi2=np.concatenate(([205,5,0.3],hi))
r=least_squares(r2,x0,bounds=(lo2,hi2),xtol=1e-12,ftol=1e-12,max_nfev=1500)
x=r.x; dz,drg,dth=x[0],x[1],x[2]; p=x[3:]
mm,Nf=ncs_forward(p,dz,drg,dth)
P=[];L=[]
for (Np_,Cp,Sp) in mm:
    P.append(0.5*np.degrees(np.arccos(np.clip(Np_,-1,1))))
    L.append(np.degrees(np.arctan2(Sp,Cp)))
P=np.array(P).T; L=np.array(L).T
dd=(L-del_m+180)%360-180
ncs_rms=np.sqrt(2*r.cost/(3*NA*len(wl)))
print('\n[O2 FINAL] d=%.2fnm rough=%.2fnm dth=%.3f  NCS-rms=%.4f  psi/del RMSE %.2f/%.2f'%(
      dz,drg,dth,ncs_rms,np.sqrt(np.mean((P-psi_m)**2)),np.sqrt(np.mean(dd**2))))
print('C0=%.3f Drude A=%.3f g=%.4f | TL A=%.1f E0=%.2f C=%.2f Eg=%.2f | G A=%.2f Ec=%.2f Br=%.2f'%tuple(p))
Ep=np.sqrt(max(p[1],1e-9)); print('plasma-ish E=sqrt(A/C0)=%.3f eV (%.0f nm)'%(np.sqrt(p[1]/p[0]),1239.84/np.sqrt(p[1]/p[0])))
for tw in [300,350,400,550,800,1000,1100,1200,1400,1650]:
    i=int(np.argmin(abs(wl-tw)))
    print('  %4dnm n=%.3f k=%.4f'%(wl[i],Nf[i].real,Nf[i].imag))
np.savez(str(OUT_DIR / 'izo_final_O2.npz'),x=x,wl=wl,psi_m=psi_m,del_m=del_m,
         psi_c=P,del_c=L,n=Nf.real,k=Nf.imag,W=W,angles=ang0,d=dz,dr=drg,dth=dth,p=p)
print('saved izo_final_O2.npz')

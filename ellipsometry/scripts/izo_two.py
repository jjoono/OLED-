"""O2_IZO3 two-layer hypothesis: top oxidized IZO (no/weak Drude) over bottom
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

conductive IZO (Drude), shared UV interband. NCS staged fit."""
import numpy as np, ellipsometry_fit as ef
from scipy.optimize import least_squares
d=np.load(str(OUT_DIR / 'izo_data.npz'))
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

def Nlayer(C0,Ad,gd,tlp,gsp):
    e2i=np.maximum(tl_eps2(Eg,*tlp)+gsp[0]*(np.exp(-((Eg-gsp[1])/gsp[2])**2)-np.exp(-((Eg+gsp[1])/gsp[2])**2)),0)
    e1g=C0+kk(e2i)-Ad/(Eg**2+gd**2)
    e2g=e2i+Ad*gd/(Eg*(Eg**2+gd**2))
    e1=np.interp(E_wl,Eg,e1g); e2=np.interp(E_wl,Eg,np.maximum(e2g,0))
    Nf=np.sqrt((e1+1j*e2).astype(complex))
    return np.where(Nf.imag<0,np.conj(Nf),Nf)

def resid(x):
    dtop,dbot,drg,dth = x[0],x[1],x[2],x[3]
    C0t,Adt = x[4],x[5]
    C0b,Adb,gdb = x[6],x[7],x[8]
    tlp=(x[9],x[10],x[11],x[12]); gsp=(x[13],x[14],x[15]); gdt=x[16]
    Nt=Nlayer(C0t,Adt,gdt,tlp,gsp)
    Nb=Nlayer(C0b,Adb,gdb,tlp,gsp)
    layers=[n_air, ef.bruggeman_ema(Nt,n_air,0.5), Nt, Nb, N_ox, N_si]
    dl=[drg,dtop,dbot,3.0]
    res=[]
    for j,ang in enumerate(ang0+dth):
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        t=np.abs(rho); p2=2*np.arctan(t); an=np.angle(rho)
        res+=[W[:,j]*(np.cos(p2)-Nm[:,j]),
              W[:,j]*(np.sin(p2)*np.cos(an)-Cm[:,j]),
              W[:,j]*(np.sin(p2)*np.sin(an)-Sm[:,j])]
    return np.concatenate(res)

#      dtop dbot  dr  dth  C0t Adt  C0b  Adb  gdb  TLA  E0  C   Eg   GA  Ec  Br  gdt
x0=[  90.0, 96.0,1.0, 0.0, 4.0,0.2, 4.3, 2.2,0.10, 30., 4.6,1.0,3.30, 0.3,3.9,0.4,0.05]
lo=[  10.0, 10.0,0.0,-0.3, 3.0,0.0, 3.0, 0.0,0.004,0.0, 3.8,0.2,2.90, 0.0,3.3,0.05,0.004]
hi=[ 180.0,180.0,5.0, 0.3, 5.5,3.0, 5.5, 6.0,0.60, 200.,7.0,3.0,4.20, 5.0,5.0,1.2, 0.6]
best=None
for st in [(90,96),(50,140),(140,50),(30,160)]:
    x0[0],x0[1]=st
    r=least_squares(resid,x0,bounds=(lo,hi),max_nfev=600)
    tot=r.x[0]+r.x[1]
    print('start %s -> dtop=%.1f dbot=%.1f (tot %.1f) cost=%.3f'%(st,r.x[0],r.x[1],tot,r.cost))
    if best is None or r.cost<best.cost: best=r
r=least_squares(resid,best.x,bounds=(lo,hi),xtol=1e-12,ftol=1e-12,max_nfev=1200)
x=r.x
ncs_rms=np.sqrt(2*r.cost/(3*NA*len(wl)))
print('\n[O2 two-layer] dtop=%.1f dbot=%.1f rough=%.2f dth=%.3f  NCS-rms=%.4f'%(x[0],x[1],x[2],x[3],ncs_rms))
print('top: C0=%.3f Ad=%.3f g=%.3f | bottom: C0=%.3f Ad=%.3f g=%.3f'%(x[4],x[5],x[16],x[6],x[7],x[8]))
print('TL A=%.1f E0=%.2f C=%.2f Eg=%.2f | G A=%.2f Ec=%.2f Br=%.2f'%(x[9],x[10],x[11],x[12],x[13],x[14],x[15]))
# psi/del RMSE
Nt=Nlayer(x[4],x[5],x[16],(x[9],x[10],x[11],x[12]),(x[13],x[14],x[15]))
Nb=Nlayer(x[6],x[7],x[8],(x[9],x[10],x[11],x[12]),(x[13],x[14],x[15]))
layers=[n_air, ef.bruggeman_ema(Nt,n_air,0.5), Nt, Nb, N_ox, N_si]
dl=[x[2],x[0],x[1],3.0]
P=[];L=[]
for ang in ang0+x[3]:
    rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
    P.append(np.degrees(np.arctan(np.abs(rho)))); L.append(np.degrees(np.angle(rho)))
P=np.array(P).T; L=np.array(L).T
dd=(L-del_m+180)%360-180
print('psi/del RMSE: %.2f / %.2f'%(np.sqrt(np.mean((P-psi_m)**2)),np.sqrt(np.mean(dd**2))))
for tw in [350,550,800,1100,1400,1650]:
    i=int(np.argmin(abs(wl-tw)))
    print('  %4dnm  top n=%.3f k=%.4f | bot n=%.3f k=%.4f'%(wl[i],Nt[i].real,Nt[i].imag,Nb[i].real,Nb[i].imag))
np.savez(str(OUT_DIR / 'izo_two_O2.npz'),x=x,wl=wl,psi_c=P,del_c=L,psi_m=psi_m,del_m=del_m,
         n_top=Nt.real,k_top=Nt.imag,n_bot=Nb.real,k_bot=Nb.imag,W=W,angles=ang0)
print('saved izo_two_O2.npz')

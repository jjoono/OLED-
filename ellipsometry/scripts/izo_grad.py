"""O2: 5-slice linearly graded Drude (A,gamma) through depth; shared C0/UV."""
import numpy as np, ellipsometry_fit as ef
from scipy.optimize import least_squares
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

d=np.load(str(OUT_DIR / 'izo_data.npz'))
wl=d['o2wl']; psi_m=d['o2psi']; del_m=d['o2del']; dep=d['o2dep']
ang0=np.array([45.0,50.0,55.0,60.0,65.0]); NA=5; NS=5
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

def resid(x,ret=False):
    dz,drg,dth,C0,At,Ab,gt,gb=x[:8]
    tlp=(x[8],x[9],x[10],x[11]); gsp=(x[12],x[13],x[14])
    e2i=np.maximum(tl_eps2(Eg,*tlp)+gsp[0]*(np.exp(-((Eg-gsp[1])/gsp[2])**2)-np.exp(-((Eg+gsp[1])/gsp[2])**2)),0)
    kki=kk(e2i)
    Ns=[]
    for s in range(NS):
        f=(s+0.5)/NS
        Ad=At+(Ab-At)*f; gd=gt+(gb-gt)*f
        e1=C0+kki-Ad/(Eg**2+gd**2)
        e2=e2i+Ad*gd/(Eg*(Eg**2+gd**2))
        e1w=np.interp(E_wl,Eg,e1); e2w=np.interp(E_wl,Eg,np.maximum(e2,0))
        Nf=np.sqrt((e1w+1j*e2w).astype(complex))
        Ns.append(np.where(Nf.imag<0,np.conj(Nf),Nf))
    layers=[n_air, ef.bruggeman_ema(Ns[0],n_air,0.5)]+Ns+[N_ox,N_si]
    dl=[drg]+[dz/NS]*NS+[3.0]
    res=[];P=[];L=[]
    for j,ang in enumerate(ang0+dth):
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        t=np.abs(rho); p2=2*np.arctan(t); an=np.angle(rho)
        if ret:
            P.append(np.degrees(np.arctan(t))); L.append(np.degrees(an))
        res+=[W[:,j]*(np.cos(p2)-Nm[:,j]),W[:,j]*(np.sin(p2)*np.cos(an)-Cm[:,j]),
              W[:,j]*(np.sin(p2)*np.sin(an)-Sm[:,j])]
    if ret: return np.array(P).T,np.array(L).T,Ns
    return np.concatenate(res)

#    d    dr  dth  C0   At   Ab   gt   gb  TLA  E0  C    Eg   GA  Ec  Br
x0=[186., 1.5,0.0, 4.2, 3.0, 2.0,0.30,0.01,90., 5.3,0.3,3.60,0.8,3.85,0.2]
lo=[160., 0.0,-0.3,3.0, 0.0, 0.0,0.004,0.004,0., 3.8,0.1,2.90,0.0,3.3,0.05]
hi=[210., 6.0, 0.3,5.5, 6.0, 6.0,0.8, 0.8, 250.,7.0,3.0,4.20,5.0,5.0,1.2]
best=None
for At0,Ab0 in [(3.0,2.0),(2.0,2.0),(0.3,2.5),(2.5,0.3)]:
    x0[4],x0[5]=At0,Ab0
    r=least_squares(resid,x0,bounds=(lo,hi),max_nfev=500)
    print('seed At/Ab=%.1f/%.1f -> cost=%.3f d=%.1f At=%.2f Ab=%.2f gt=%.3f gb=%.3f'%(
          At0,Ab0,r.cost,r.x[0],r.x[4],r.x[5],r.x[6],r.x[7]))
    if best is None or r.cost<best.cost: best=r
r=least_squares(resid,best.x,bounds=(lo,hi),xtol=1e-12,ftol=1e-12,max_nfev=1200)
x=r.x
ncs=np.sqrt(2*r.cost/(3*NA*len(wl)))
P,L,Ns=resid(x,ret=True)
dd=(L-del_m+180)%360-180
print('\n[graded] d=%.1f rough=%.2f dth=%.3f C0=%.3f  NCS-rms=%.4f  psi/del %.2f/%.2f'%(
      x[0],x[1],x[2],x[3],ncs,np.sqrt(np.mean((P-psi_m)**2)),np.sqrt(np.mean(dd**2))))
print('Drude top A=%.2f g=%.3f  ->  bottom A=%.2f g=%.3f'%(x[4],x[6],x[5],x[7]))
print('TL A=%.1f E0=%.2f C=%.2f Eg=%.2f | G A=%.2f Ec=%.2f Br=%.2f'%tuple(x[8:15]))
np.savez(str(OUT_DIR / 'izo_grad_O2.npz'),x=x,wl=wl,psi_c=P,del_c=L,psi_m=psi_m,del_m=del_m,W=W)
print('saved')

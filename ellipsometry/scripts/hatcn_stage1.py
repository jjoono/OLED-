"""Stage 1: bare HATCN sheets. HATCN30 -> full (d + genosc eps) fit;
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

then HATCN2/4/6 -> thickness-only with fixed eps. NCS space, depol weights."""
import numpy as np, ellipsometry_fit as ef
from scipy.optimize import least_squares
d=np.load(str(OUT_DIR / 'aghatcn_data.npz'))
ang0=np.array([45.0,50.0,55.0,60.0,65.0]); NA=5

def unpack_sheet(A):
    wl=A[:,0]; ok=np.isfinite(wl)
    A=A[ok]; wl=A[:,0]
    P=A[:,1:11:2]; D=A[:,2:11:2]; dep=A[:,45:50]
    return wl,P,D,dep

wl,_,_,_=unpack_sheet(d['HATCN30'])
E_wl=1239.84193/wl
N_si=np.sqrt((np.interp(wl,*(lambda z: (z['si_wl'],z['si_e1']+1j*z['si_e2']))(np.load(str(OUT_DIR / 'izo_data.npz'))))).astype(complex)) if False else None
z=np.load(str(OUT_DIR / 'izo_data.npz'))
N_si=np.sqrt((np.interp(wl,z['si_wl'],z['si_e1'])+1j*np.interp(wl,z['si_wl'],z['si_e2'])).astype(complex))
N_ox=np.sqrt((np.interp(wl,z['ox_wl'],z['ox_e1'])+1j*np.interp(wl,z['ox_wl'],z['ox_e2'])).astype(complex))
n_air=np.ones(len(wl),dtype=complex)

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

def eps_hatcn(p):
    C0=p[0]
    e2=tl_eps2(Eg,p[1],p[2],p[3],p[4])
    e2=e2+p[5]*(np.exp(-((Eg-p[6])/p[7])**2)-np.exp(-((Eg+p[6])/p[7])**2))
    e2=np.maximum(e2,0)
    e1=C0+kk(e2)
    e1w=np.interp(E_wl,Eg,e1); e2w=np.interp(E_wl,Eg,np.maximum(e2,0))
    Nf=np.sqrt((e1w+1j*e2w).astype(complex))
    return np.where(Nf.imag<0,np.conj(Nf),Nf)

def ncs(P,D):
    pr=np.deg2rad(P); drr=np.deg2rad(D)
    return np.cos(2*pr),np.sin(2*pr)*np.cos(drr),np.sin(2*pr)*np.sin(drr)

def model_ncs(Nf,dz):
    layers=[n_air,Nf,N_ox,N_si]; dl=[dz,3.0]
    out=[]
    for ang in ang0:
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        t=np.abs(rho); p2=2*np.arctan(t); an=np.angle(rho)
        out.append((np.cos(p2),np.sin(p2)*np.cos(an),np.sin(p2)*np.sin(an)))
    return out

def fit_sheet_full(sheet,d0):
    wlS,P,D,dep=unpack_sheet(d[sheet])
    Nm,Cm,Sm=ncs(P,D)
    W=1.0/(1.0+(np.abs(dep)/3.0)**2)
    def res(x):
        Nf=eps_hatcn(x[1:])
        mm=model_ncs(Nf,x[0])
        r=[]
        for j,(a,b,c) in enumerate(mm):
            r+=[W[:,j]*(a-Nm[:,j]),W[:,j]*(b-Cm[:,j]),W[:,j]*(c-Sm[:,j])]
        return np.concatenate(r)
    #        d    C0   TLA  E0   C    Eg   GA   Ec   Br
    x0=[d0, 1.9, 6.0, 3.9, 0.5, 3.1, 0.5, 4.6, 0.4]
    lo=[2.0, 1.2, 0.0, 3.2, 0.1, 2.6, 0.0, 3.8, 0.1]
    hi=[40., 3.5, 60., 5.0, 2.0, 3.6, 8.0, 6.0, 1.5]
    best=None
    for dd0 in [d0*0.6,d0,d0*1.3]:
        x0[0]=dd0
        r=least_squares(res,x0,bounds=(lo,hi),max_nfev=400)
        if best is None or r.cost<best.cost: best=r
    r=least_squares(res,best.x,bounds=(lo,hi),xtol=1e-12,max_nfev=600)
    rms=np.sqrt(2*r.cost/(3*NA*len(wlS)))
    return r.x,rms

def fit_sheet_donly(sheet,epsp,d0):
    wlS,P,D,dep=unpack_sheet(d[sheet])
    Nm,Cm,Sm=ncs(P,D)
    W=1.0/(1.0+(np.abs(dep)/3.0)**2)
    Nf=eps_hatcn(epsp)
    def res(x):
        mm=model_ncs(Nf,x[0])
        r=[]
        for j,(a,b,c) in enumerate(mm):
            r+=[W[:,j]*(a-Nm[:,j]),W[:,j]*(b-Cm[:,j]),W[:,j]*(c-Sm[:,j])]
        return np.concatenate(r)
    best=None
    for dd0 in [d0*0.5,d0,d0*1.5,0.5]:
        r=least_squares(res,[max(dd0,0.3)],bounds=([0.0],[15.0]),xtol=1e-12)
        if best is None or r.cost<best.cost: best=r
    rms=np.sqrt(2*best.cost/(3*NA*len(wlS)))
    return float(best.x[0]),rms

x30,rms30=fit_sheet_full('HATCN30',18.0)
print('[HATCN30] d=%.2fnm rms=%.4f  C0=%.3f TL(A=%.1f E0=%.2f C=%.2f Eg=%.2f) G(A=%.2f Ec=%.2f Br=%.2f)'%(
      x30[0],rms30,*x30[1:]))
Nf=eps_hatcn(x30[1:])
for tw in [300,365,400,450,550,633,800,1200]:
    i=int(np.argmin(abs(wl-tw)))
    print('   %4dnm n=%.3f k=%.4f'%(wl[i],Nf[i].real,Nf[i].imag))
print('tooling: 30nm nominal -> %.2fnm = %.0f%%'%(x30[0],100*x30[0]/30))
res={}
for sheet,nom in [('HATCN2',2.0),('HATCN4',4.0),('HATCN6',6.0)]:
    dz,rms=fit_sheet_donly(sheet,x30[1:],nom*0.6)
    res[sheet]=(dz,rms)
    print('[%s] d=%.2fnm (nominal %.0f, %.0f%%) rms=%.4f'%(sheet,dz,nom,100*dz/nom,rms))
np.savez(str(OUT_DIR / 'hatcn_stage1.npz'),x30=x30,wl=wl,n=Nf.real,k=Nf.imag,
         d30=x30[0],d2=res['HATCN2'][0],d4=res['HATCN4'][0],d6=res['HATCN6'][0])
print('saved hatcn_stage1.npz')

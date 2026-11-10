"""Stage 1b: bare HATCN with hard cut wl<=1050nm (Si-backside excluded).
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

(1) HATCN30 full fit; (2) n-branch degeneracy scan (fix n633, fit d);
(3) HATCN2/4/6 d-only at each branch."""
import numpy as np, ellipsometry_fit as ef
from scipy.optimize import least_squares
d=np.load(str(OUT_DIR / 'aghatcn_data.npz'))
z=np.load(str(OUT_DIR / 'izo_data.npz'))
ang0=np.array([45.0,50.0,55.0,60.0,65.0]); NA=5

def unpack(A,wlmax=1050.0):
    wl=A[:,0]; ok=np.isfinite(wl)&(wl<=wlmax)
    A=A[ok]; wl=A[:,0]
    return wl,A[:,1:11:2],A[:,2:11:2],A[:,45:50]

wl,_,_,_=unpack(d['HATCN30'])
E_wl=1239.84193/wl
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
    e2=tl_eps2(Eg,p[1],p[2],p[3],p[4])+p[5]*(np.exp(-((Eg-p[6])/p[7])**2)-np.exp(-((Eg+p[6])/p[7])**2))
    e2=np.maximum(e2,0)
    e1=C0+kk(e2)
    e1w=np.interp(E_wl,Eg,e1); e2w=np.interp(E_wl,Eg,np.maximum(e2,0))
    Nf=np.sqrt((e1w+1j*e2w).astype(complex))
    return np.where(Nf.imag<0,np.conj(Nf),Nf)

def sheet_ncs(sheet):
    wlS,P,D,dep=unpack(d[sheet])
    pr=np.deg2rad(P); drr=np.deg2rad(D)
    W=1.0/(1.0+(np.abs(dep)/3.0)**2)
    return np.cos(2*pr),np.sin(2*pr)*np.cos(drr),np.sin(2*pr)*np.sin(drr),W

def model_ncs(Nf,dz):
    layers=[n_air,Nf,N_ox,N_si]; dl=[dz,3.0]
    out=[]
    for ang in ang0:
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        t=np.abs(rho); p2=2*np.arctan(t); an=np.angle(rho)
        out.append((np.cos(p2),np.sin(p2)*np.cos(an),np.sin(p2)*np.sin(an)))
    return out

def resfun(sheet):
    Nm,Cm,Sm,W=sheet_ncs(sheet)
    def res(x,epsp=None,dfix=None):
        if epsp is None:
            dz=x[0]; Nf=eps_hatcn(x[1:])
        else:
            dz=x[0]; Nf=eps_hatcn(epsp)
        mm=model_ncs(Nf,dz)
        r=[]
        for j,(a,b,c) in enumerate(mm):
            r+=[W[:,j]*(a-Nm[:,j]),W[:,j]*(b-Cm[:,j]),W[:,j]*(c-Sm[:,j])]
        return np.concatenate(r)
    return res,Nm.shape[0]

# (1) HATCN30 full
res30,npts=resfun('HATCN30')
x0=[18., 1.9, 6.0, 3.9, 0.5, 3.1, 0.5, 4.6, 0.4]
lo=[5.0, 1.2, 0.0, 3.2, 0.1, 2.6, 0.0, 3.8, 0.1]
hi=[45., 3.5, 60., 5.0, 2.0, 3.6, 8.0, 6.0, 1.5]
best=None
for d0 in [12,18,25,32]:
    x0[0]=d0
    r=least_squares(res30,x0,bounds=(lo,hi),max_nfev=400)
    if best is None or r.cost<best.cost: best=r
r=least_squares(res30,best.x,bounds=(lo,hi),xtol=1e-12,max_nfev=600)
x30=r.x
rms=np.sqrt(2*r.cost/(3*NA*npts))
Nf=eps_hatcn(x30[1:])
i633=int(np.argmin(abs(wl-633)))
print('[HATCN30 cut1050] d=%.2f rms=%.4f n633=%.3f  (tooling %.0f%%)'%(x30[0],rms,Nf[i633].real,100*x30[0]/30))

# (2) n-branch scan: set C0 analytically for target n633, fit d only
print('\n-- n-branch degeneracy scan (HATCN30) --')
e2_fixed=np.maximum(tl_eps2(Eg,*x30[1+1:1+5])+x30[6]*(np.exp(-((Eg-x30[7])/x30[8])**2)
          -np.exp(-((Eg+x30[7])/x30[8])**2)),0.0)
kk633=float(np.interp(1239.84193/633.0,Eg,kk(e2_fixed)))
branches={}
for n633 in [1.78,1.90,2.00,2.10,2.25]:
    C0b=n633**2-kk633
    pb=x30[1:].copy(); pb[0]=C0b
    b=None
    for d0 in [12,18,24,30,38]:
        r=least_squares(lambda xx: res30(np.concatenate((xx,pb))),[d0],
                        bounds=([3.0],[60.0]),xtol=1e-12)
        if b is None or r.cost<b.cost: b=r
    rm=np.sqrt(2*b.cost/(3*NA*npts))
    branches[n633]=float(b.x[0])
    print('  n633=%.2f (C0=%.3f): d=%.2fnm  rms=%.4f'%(n633,C0b,float(b.x[0]),rm))

# (3) HATCN2/4/6 d-only at main-fit eps and at n=1.90 branch C0
print('\n-- thin sheets d-only (cut 1050) --')
for sheet,nom in [('HATCN2',2.0),('HATCN4',4.0),('HATCN6',6.0)]:
    resS,nptsS=resfun(sheet)
    b=None
    for d0 in [0.5,1.5,3,6,10]:
        r=least_squares(lambda xx: resS(np.concatenate((xx,x30[1:]))),[d0],bounds=([0.0],[20.0]),xtol=1e-12)
        if b is None or r.cost<b.cost: b=r
    rmsS=np.sqrt(2*b.cost/(3*NA*nptsS))
    print('  %s: d=%.2fnm (nom %.0f) rms=%.4f'%(sheet,float(b.x[0]),nom,rmsS))
np.savez(str(OUT_DIR / 'hatcn_stage1b.npz'),x30=x30,wl=wl,n=Nf.real,k=Nf.imag)

"""Per-wavelength free (n,k) extraction for IZO: scan (d,rough) for 4-angle
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

consistency, then extract full nk(lambda) at the best geometry."""
import sys, numpy as np, ellipsometry_fit as ef
from scipy.optimize import least_squares
d=np.load(str(OUT_DIR / 'izo_data.npz'))
which=sys.argv[1] if len(sys.argv)>1 else 'N2'
if which=='N2':
    wl=d['n2wl']; psi_m=d['n2psi']; del_m=d['n2del']; dep=d['n2dep']
    angles=np.array([45.0,55.0,65.0,75.0])
else:
    wl=d['o2wl']; psi_m=d['o2psi']; del_m=d['o2del']; dep=d['o2dep']
    angles=np.array([45.0,50.0,55.0,60.0,65.0])
N_si=np.sqrt((np.interp(wl,d['si_wl'],d['si_e1'])+1j*np.interp(wl,d['si_wl'],d['si_e2'])).astype(complex))
N_ox=np.sqrt((np.interp(wl,d['ox_wl'],d['ox_e1'])+1j*np.interp(wl,d['ox_wl'],d['ox_e2'])).astype(complex))
W=1.0/(1.0+(np.abs(dep)/3.0)**2)

def pd_one(i,n,k,dz,dr):
    Nf=np.array([complex(n,max(k,0.0))])
    na=np.ones(1,dtype=complex)
    layers=[na, ef.bruggeman_ema(Nf,na,0.5), Nf, N_ox[i:i+1], N_si[i:i+1]]
    dl=[dr,dz,3.0]
    P=[];L=[]
    for ang in angles:
        rp,rs=ef._tmm(wl[i:i+1],layers,dl,ang); rho=rp/rs
        P.append(float(np.degrees(np.arctan(np.abs(rho[0])))))
        L.append(float(np.degrees(np.angle(rho[0]))))
    return np.array(P),np.array(L)

def extract_one(i,dz,dr,seeds):
    def f(x):
        P,L=pd_one(i,x[0],x[1],dz,dr)
        dd=(L-del_m[i]+180)%360-180
        return np.concatenate([W[i]*(P-psi_m[i]),W[i]*dd/3.0])
    best=None
    for n0,k0 in seeds:
        r=least_squares(f,[n0,k0],bounds=([0.3,0],[4.5,3.5]),xtol=1e-10,ftol=1e-10)
        if best is None or r.cost<best.cost: best=r
    return best

# stage A: (d, rough) scan on subsample
idx=np.arange(0,len(wl),8)
seeds=[(2.1,0.0),(1.6,0.3),(2.6,0.5),(1.2,1.0),(3.2,0.8)]
print('scan (d,rough) -> mean per-wl weighted resid')
best_geo=None
for dz in range(150,206,5):
    for dr in [0.0,3.0,6.0,10.0,15.0]:
        tot=0.0
        for i in idx:
            b=extract_one(i,dz,dr,seeds)
            tot+=np.sqrt(2*b.cost/ (2*len(angles)))
        m=tot/len(idx)
        if best_geo is None or m<best_geo[0]:
            best_geo=(m,dz,dr); print('  d=%d dr=%.0f: %.3f  *'%(dz,dr,m))
print('BEST geometry: d=%dnm rough=%.0fnm (resid %.3f)'%(best_geo[1],best_geo[2],best_geo[0]))

# stage B: full extraction at best geometry (seed-chain)
dz,dr=best_geo[1],best_geo[2]
ns=np.zeros(len(wl)); ks=np.zeros(len(wl)); rs=np.zeros(len(wl))
prev=(2.1,0.0)
for i in range(len(wl)):
    b=extract_one(i,dz,dr,[prev]+seeds)
    ns[i],ks[i]=b.x; rs[i]=np.sqrt(2*b.cost/(2*len(angles)))
    prev=(ns[i],ks[i])
print('mean resid=%.3f  median=%.3f'%(rs.mean(),np.median(rs)))
for tw in [300,350,400,500,600,800,1000,1200,1400,1600]:
    i=int(np.argmin(abs(wl-tw)))
    print('  %4dnm n=%.3f k=%.4f (res %.2f)'%(wl[i],ns[i],ks[i],rs[i]))
np.savez(str(OUT_DIR / 'izo_perwl_%s.npz')%which,wl=wl,n=ns,k=ks,res=rs,d=dz,dr=dr)
print('saved izo_perwl_%s.npz'%which)

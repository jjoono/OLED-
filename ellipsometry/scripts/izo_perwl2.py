"""Branch-disciplined per-wavelength nk extraction: descending-lambda chain with
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

continuity penalty; physical bounds exclude the film~substrate degenerate branch."""
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
NA=len(angles)

def pd_one(i,n,k,dz,dr):
    Nf=np.array([complex(n,max(k,0.0))])
    na=np.ones(1,dtype=complex)
    layers=[na, ef.bruggeman_ema(Nf,na,0.5), Nf, N_ox[i:i+1], N_si[i:i+1]]
    dl=[dr,dz,3.0]
    P=np.empty(NA);L=np.empty(NA)
    for j,ang in enumerate(angles):
        rp,rs=ef._tmm(wl[i:i+1],layers,dl,ang); rho=rp/rs
        P[j]=np.degrees(np.arctan(np.abs(rho[0]))); L[j]=np.degrees(np.angle(rho[0]))
    return P,L

def chain(dz,dr,il,alpha=8.0,n_start=2.0):
    ns=np.zeros(len(il)); ks=np.zeros(len(il)); rs=np.zeros(len(il))
    prev=(n_start,0.0)
    for m,i in enumerate(il):
        def f(x):
            P,L=pd_one(i,x[0],x[1],dz,dr)
            dd=(L-del_m[i]+180)%360-180
            cont=alpha*np.array([x[0]-prev[0],x[1]-prev[1]])
            return np.concatenate([W[i]*(P-psi_m[i]),W[i]*dd/3.0,cont])
        b=None
        for n0,k0 in [prev,(prev[0]+0.06,prev[1]),(prev[0]-0.06,prev[1]),(prev[0],prev[1]+0.05)]:
            r=least_squares(f,[n0,k0],bounds=([0.8,0.0],[3.4,2.0]),xtol=1e-10,ftol=1e-10)
            if b is None or r.cost<b.cost: b=r
        ns[m],ks[m]=b.x
        P,L=pd_one(i,ns[m],ks[m],dz,dr)
        dd=(L-del_m[i]+180)%360-180
        rs[m]=np.sqrt(np.mean(np.concatenate([W[i]*(P-psi_m[i]),W[i]*dd/3.0])**2))
        prev=(ns[m],ks[m])
    return ns,ks,rs

# geometry scan on subsample, descending lambda
il_sub=np.arange(len(wl)-1,-1,-10)
print('geometry scan (chain-resid mean):')
best=None
for dz in range(150,206,5):
    for dr in [0.0,2.0,4.0,7.0,10.0]:
        ns,ks,rs=chain(dz,dr,il_sub)
        m=rs.mean()
        if best is None or m<best[0]:
            best=(m,dz,dr); print('  d=%d dr=%.0f: %.3f *'%(dz,dr,m))
print('BEST: d=%d dr=%.0f (%.3f)'%(best[1],best[2],best[0]))
dz,dr=best[1],best[2]
il=np.arange(len(wl)-1,-1,-1)
ns,ks,rs=chain(dz,dr,il)
order=np.argsort(wl[il])
wlo=wl[il][order]; no=ns[order]; ko=ks[order]; ro=rs[order]
print('full chain: mean res=%.3f median=%.3f'%(ro.mean(),np.median(ro)))
for tw in [300,350,400,500,600,800,1000,1200,1400,1600]:
    i=int(np.argmin(abs(wlo-tw)))
    print('  %4dnm n=%.3f k=%.4f (res %.2f)'%(wlo[i],no[i],ko[i],ro[i]))
np.savez(str(OUT_DIR / 'izo_perwl_%s.npz')%which,wl=wlo,n=no,k=ko,res=ro,d=dz,dr=dr)
print('saved izo_perwl_%s.npz'%which)

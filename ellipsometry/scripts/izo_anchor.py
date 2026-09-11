"""O2: (1) window-cost vs fixed d (alias check), (2) chain nk extraction at
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

externally anchored d=179nm (IZO3 pre-anneal CE thickness)."""
import numpy as np, ellipsometry_fit as ef
from scipy.optimize import least_squares
d=np.load(str(OUT_DIR / 'izo_data.npz'))
wl=d['o2wl']; psi_m=d['o2psi']; del_m=d['o2del']; dep=d['o2dep']
ang0=np.array([45.0,50.0,55.0,60.0,65.0]); NA=5
N_si=np.sqrt((np.interp(wl,d['si_wl'],d['si_e1'])+1j*np.interp(wl,d['si_wl'],d['si_e2'])).astype(complex))
N_ox=np.sqrt((np.interp(wl,d['ox_wl'],d['ox_e1'])+1j*np.interp(wl,d['ox_wl'],d['ox_e2'])).astype(complex))
W=1.0/(1.0+(np.abs(dep)/3.0)**2)
pr=np.deg2rad(psi_m); dr_=np.deg2rad(del_m)
Nm=np.cos(2*pr); Cm=np.sin(2*pr)*np.cos(dr_); Sm=np.sin(2*pr)*np.sin(dr_)

def ncs_res_window(x,dfix,mask):
    drg,dth,A,B=x
    Nf=((A+B/wl**2)+0j).astype(complex)
    na=np.ones(len(wl),dtype=complex)
    layers=[na, ef.bruggeman_ema(Nf,na,0.5), Nf, N_ox, N_si]
    dl=[drg,dfix,3.0]
    res=[]
    for j,ang in enumerate(ang0+dth):
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        t=np.abs(rho); p2=2*np.arctan(t); an=np.angle(rho)
        res+=[ (W[:,j]*(np.cos(p2)-Nm[:,j]))[mask],
               (W[:,j]*(np.sin(p2)*np.cos(an)-Cm[:,j]))[mask],
               (W[:,j]*(np.sin(p2)*np.sin(an)-Sm[:,j]))[mask] ]
    return np.concatenate(res)

mask=(wl>=550)&(wl<=1100)
print('== window cost vs fixed d (alias flatness) ==')
for dfix in [167,171,175,179,183,187,191,195,200,206]:
    b=None
    for A0 in [1.9,2.05,2.2]:
        r=least_squares(ncs_res_window,[1.0,0.0,A0,2e4],args=(dfix,mask),
                        bounds=([0,-0.2,1.5,0],[4,0.2,2.7,2e5]))
        if b is None or r.cost<b.cost: b=r
    print('  d=%3d: cost=%.4f  n800=%.3f dth=%+.3f dr=%.1f'%(dfix,b.cost,b.x[2]+b.x[3]/800**2,b.x[1],b.x[0]))

# ---- chain extraction at d=179.1 ----
DZ,DR=179.1,1.0
def pd_one(i,n,k):
    Nf=np.array([complex(n,max(k,0.0))])
    na=np.ones(1,dtype=complex)
    layers=[na, ef.bruggeman_ema(Nf,na,0.5), Nf, N_ox[i:i+1], N_si[i:i+1]]
    dl=[DR,DZ,3.0]
    N_=np.empty(NA);C_=np.empty(NA);S_=np.empty(NA)
    for j,ang in enumerate(ang0):
        rp,rs=ef._tmm(wl[i:i+1],layers,dl,ang); rho=rp/rs
        t=abs(rho[0]); p2=2*np.arctan(t); an=np.angle(rho[0])
        N_[j]=np.cos(p2);C_[j]=np.sin(p2)*np.cos(an);S_[j]=np.sin(p2)*np.sin(an)
    return N_,C_,S_

alpha=2.0
il=np.arange(len(wl)-1,-1,-1)
ns=np.zeros(len(wl));ks=np.zeros(len(wl));rs=np.zeros(len(wl))
prev=(1.85,0.0)
for m,i in enumerate(il):
    def f(x):
        N_,C_,S_=pd_one(i,x[0],x[1])
        cont=alpha*np.array([x[0]-prev[0],x[1]-prev[1]])
        return np.concatenate([W[i]*(N_-Nm[i]),W[i]*(C_-Cm[i]),W[i]*(S_-Sm[i]),cont])
    b=None
    for s in [prev,(prev[0]+0.05,prev[1]),(prev[0]-0.05,prev[1]),(prev[0],prev[1]+0.04)]:
        r=least_squares(f,s,bounds=([1.1,0.0],[2.9,0.9]),xtol=1e-10,ftol=1e-10)
        if b is None or r.cost<b.cost: b=r
    ns[m],ks[m]=b.x
    N_,C_,S_=pd_one(i,ns[m],ks[m])
    rs[m]=np.sqrt(np.mean(np.concatenate([W[i]*(N_-Nm[i]),W[i]*(C_-Cm[i]),W[i]*(S_-Sm[i])])**2))
    prev=(ns[m],ks[m])
o=np.argsort(wl[il]); wlo=wl[il][o]; no=ns[o]; ko=ks[o]; ro=rs[o]
print('\nchain @ d=179.1: mean NCS-res=%.4f median=%.4f'%(ro.mean(),np.median(ro)))
for tw in [280,300,320,350,400,500,600,800,1000,1200,1400,1650]:
    i=int(np.argmin(abs(wlo-tw)))
    print('  %4dnm n=%.3f k=%.4f (res %.3f)'%(wlo[i],no[i],ko[i],ro[i]))
np.savez(str(OUT_DIR / 'izo_chain_O2.npz'),wl=wlo,n=no,k=ko,res=ro,d=DZ,dr=DR)
print('saved izo_chain_O2.npz')

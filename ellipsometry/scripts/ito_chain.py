"""ITO model-free chain nk at S1 geometry, both samples, + overlay comparison."""
import numpy as np, openpyxl, ellipsometry_fit as ef, csv
from scipy.optimize import least_squares
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

exec(open(str(OUT_DIR / 'ito_fit.py')).read().split('Eg=np.linspace')[0])  # loaders+data
# (re-uses wl1,P1,D1,dep1, wl2,P2,D2,dep2, z, ang0, NA)

def chain(tag,wl_,P_,D_,dep_,DZ,DR):
    m=(wl_>=250)&(wl_<=1080)
    wl=wl_[m]; P=P_[m]; D=D_[m]; dep=dep_[m]
    N_si=np.sqrt((np.interp(wl,z['si_wl'],z['si_e1'])+1j*np.interp(wl,z['si_wl'],z['si_e2'])).astype(complex))
    N_ox=np.sqrt((np.interp(wl,z['ox_wl'],z['ox_e1'])+1j*np.interp(wl,z['ox_wl'],z['ox_e2'])).astype(complex))
    pr=np.deg2rad(P); drr=np.deg2rad(D)
    Nm=np.cos(2*pr); Cm=np.sin(2*pr)*np.cos(drr); Sm=np.sin(2*pr)*np.sin(drr)
    W=1.0/(1.0+(np.abs(dep)/3.0)**2)
    def pd_one(i,n,k):
        Nf=np.array([complex(n,max(k,0.0))])
        na=np.ones(1,dtype=complex)
        layers=[na, ef.bruggeman_ema(Nf,na,0.5), Nf, N_ox[i:i+1], N_si[i:i+1]]
        dl=[DR,DZ,3.0]
        NN=np.empty(3);CC=np.empty(3);SS=np.empty(3)
        for j,ang in enumerate(ang0):
            rp,rs=ef._tmm(wl[i:i+1],layers,dl,ang); rho=rp/rs
            t=abs(rho[0]); p2=2*np.arctan(t); an=np.angle(rho[0])
            NN[j]=np.cos(p2);CC[j]=np.sin(p2)*np.cos(an);SS[j]=np.sin(p2)*np.sin(an)
        return NN,CC,SS
    il=np.arange(len(wl)-1,-1,-1)
    ns=np.zeros(len(wl));ks=np.zeros(len(wl));rs=np.zeros(len(wl))
    prev=np.array([1.75,0.02]); alpha=1.5
    LOB=np.array([1.0,0.0]); HIB=np.array([3.2,1.5])
    for mth,i in enumerate(il):
        def f(x):
            NN,CC,SS=pd_one(i,x[0],x[1])
            cont=alpha*np.array([x[0]-prev[0],x[1]-prev[1]])
            return np.concatenate([W[i]*(NN-Nm[i]),W[i]*(CC-Cm[i]),W[i]*(SS-Sm[i]),cont])
        b=None
        for s in [prev,prev+[0.05,0],prev-[0.05,0],prev+[0,0.03]]:
            s=np.clip(s,LOB+1e-6,HIB-1e-6)
            r=least_squares(f,s,bounds=(LOB,HIB),xtol=1e-10,ftol=1e-10)
            if b is None or r.cost<b.cost: b=r
        ns[mth],ks[mth]=b.x
        NN,CC,SS=pd_one(i,ns[mth],ks[mth])
        rs[mth]=np.sqrt(np.mean(np.concatenate([W[i]*(NN-Nm[i]),W[i]*(CC-Cm[i]),W[i]*(SS-Sm[i])])**2))
        prev=b.x.copy()
    o=np.argsort(wl[il])
    wlo=wl[il][o]; no=ns[o]; ko=ks[o]; ro=rs[o]
    print('%s chain @ d=%.2f/rough=%.2f: mean res=%.4f median=%.4f'%(tag,DZ,DR,ro.mean(),np.median(ro)))
    for tw in [280,300,330,360,400,450,550,633,800,1000]:
        i=int(np.argmin(abs(wlo-tw)))
        print('   %4dnm n=%.3f k=%.4f (res %.3f)'%(wlo[i],no[i],ko[i],ro[i]))
    return wlo,no,ko,ro

w1,n1,k1,r1=chain('#1',wl1,P1,D1,dep1,46.86,4.11)
w2,n2,k2,r2=chain('#2',wl2,P2,D2,dep2,49.86,5.11)
print('\n== batch overlay (interp #2 onto #1 grid) ==')
n2i=np.interp(w1,w2,n2); k2i=np.interp(w1,w2,k2)
m=(w1>320)&(w1<1050)
print('mean|dn|=%.4f  mean|dk|=%.4f  (320-1050nm)'%(np.mean(np.abs(n1[m]-n2i[m])),np.mean(np.abs(k1[m]-k2i[m]))))
with open(str(OUT_DIR / 'ITO_nk_chain.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['wl_nm','n_1','k_1','res_1','n_2','k_2'])
    for i in range(len(w1)):
        w.writerow(['%.2f'%w1[i],'%.4f'%n1[i],'%.4f'%k1[i],'%.4f'%r1[i],'%.4f'%n2i[i],'%.4f'%k2i[i]])
np.savez(str(OUT_DIR / 'ito_chain.npz'),w1=w1,n1=n1,k1=k1,r1=r1,w2=w2,n2=n2,k2=k2,r2=r2)
print('saved ITO_nk_chain.csv / ito_chain.npz')

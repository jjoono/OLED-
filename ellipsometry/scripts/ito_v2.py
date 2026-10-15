"""Re-verdict on modified #1 (now interleaved) + final comparison with #2."""
import numpy as np, openpyxl, ellipsometry_fit as ef, csv
from scipy.optimize import least_squares
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

fp=str(DATA_DIR / 'se추출.xlsx')
z=np.load(str(OUT_DIR / 'izo_data.npz'))
ang0=np.array([65.0,70.0,75.0]); NA=3
def tofloat(v):
    if isinstance(v,(int,float)): return float(v)
    try: return float(str(v).strip())
    except: return np.nan
wb=openpyxl.load_workbook(fp,data_only=True,read_only=True)
def read_sheet(name):
    ws=wb[name]
    hdr=None
    for i,row in enumerate(ws.iter_rows(min_row=3,max_row=3,values_only=True)):
        hdr=row
    depcols=[j for j,v in enumerate(hdr) if v and 'Depolarization' in str(v)]
    rows=[]
    for row in ws.iter_rows(min_row=4,values_only=True):
        rows.append([tofloat(v) for v in (list(row)+[None]*120)[:120]])
    A=np.array(rows); ok=np.isfinite(A[:,0]); A=A[ok]
    wl=A[:,0]; P=A[:,1:7:2]; D=A[:,2:7:2]
    dep=A[:,depcols] if len(depcols)==3 else np.zeros((len(wl),3))
    return wl,P,D,dep,depcols
wl1,P1,D1,dep1,dc1=read_sheet('#1')
wl2,P2,D2,dep2,dc2=read_sheet('Sheet2')
wb.close()
print('#1: %d pts %.0f-%.0f depcols=%s | #2: %d pts depcols=%s'%(len(wl1),wl1[0],wl1[-1],dc1,len(wl2),dc2))

def prep(wl_,P_,D_,dep_,lo,hi):
    m=(wl_>=lo)&(wl_<=hi)
    wl=wl_[m]; P=P_[m]; D=D_[m]; dep=dep_[m]
    N_si=np.sqrt((np.interp(wl,z['si_wl'],z['si_e1'])+1j*np.interp(wl,z['si_wl'],z['si_e2'])).astype(complex))
    N_ox=np.sqrt((np.interp(wl,z['ox_wl'],z['ox_e1'])+1j*np.interp(wl,z['ox_wl'],z['ox_e2'])).astype(complex))
    pr=np.deg2rad(P); drr=np.deg2rad(D)
    return (wl,np.cos(2*pr),np.sin(2*pr)*np.cos(drr),np.sin(2*pr)*np.sin(drr),
            1.0/(1.0+(np.abs(dep)/3.0)**2),N_si,N_ox)

def chain(tag,wl_,P_,D_,dep_,DZ,DR,scan_only=False,il_stride=1):
    wl,Nm,Cm,Sm,W,N_si,N_ox=prep(wl_,P_,D_,dep_,360 if scan_only else 250,1080)
    il=np.arange(len(wl)-1,-1,-il_stride)
    prev=np.array([1.9,0.01]); ns=[];ks=[];rs=[];wlv=[]
    for i in il:
        def f(x):
            Nf=np.array([complex(x[0],max(x[1],0.0))])
            na=np.ones(1,dtype=complex)
            layers=[na, ef.bruggeman_ema(Nf,na,0.5), Nf, N_ox[i:i+1], N_si[i:i+1]]
            dl=[DR,DZ,3.0]
            NN=np.empty(3);CC=np.empty(3);SS=np.empty(3)
            for j,ang in enumerate(ang0):
                rp,rs_=ef._tmm(wl[i:i+1],layers,dl,ang); rho=rp/rs_
                t=abs(rho[0]); p2=2*np.arctan(t); an=np.angle(rho[0])
                NN[j]=np.cos(p2);CC[j]=np.sin(p2)*np.cos(an);SS[j]=np.sin(p2)*np.sin(an)
            cont=1.5*np.array([x[0]-prev[0],x[1]-prev[1]])
            return np.concatenate([W[i]*(NN-Nm[i]),W[i]*(CC-Cm[i]),W[i]*(SS-Sm[i]),cont])
        b=None
        for s in [prev,prev+[0.06,0],prev-[0.06,0]]:
            s=np.clip(s,[1.0,0.0],[3.2,1.5])
            r=least_squares(f,s,bounds=([1.0,0.0],[3.2,1.5]),xtol=1e-9)
            if b is None or r.cost<b.cost: b=r
        prev=b.x.copy()
        ns.append(b.x[0]);ks.append(b.x[1]);rs.append(np.sqrt(2*b.cost/6));wlv.append(wl[i])
    o=np.argsort(wlv)
    return np.array(wlv)[o],np.array(ns)[o],np.array(ks)[o],np.array(rs)[o]

# geometry scan for new #1
print('-- #1 geometry scan (360-1080 score) --')
best=None
for DZ in [44,46,48,50,52,54]:
    for DR in [2.0,5.0]:
        _,_,_,r=chain('#1',wl1,P1,D1,dep1,DZ,DR,scan_only=True,il_stride=5)
        v=r.mean()
        if best is None or v<best[0]:
            best=(v,DZ,DR); print('  d=%.0f dr=%.0f: %.4f *'%(DZ,DR,v))
print('#1 BEST: d=%.1f rough=%.1f res=%.4f (수정 전 평탄바닥 0.0874)'%(best[1],best[2],best[0]))
DZ1,DR1=best[1],best[2]
w1,n1,k1,r1=chain('#1',wl1,P1,D1,dep1,DZ1,DR1)
print('#1 full chain: median res=%.4f'%np.median(r1))
f2=np.load(str(OUT_DIR / 'ito_final.npz'))
w2,n2,k2=f2['w2'],f2['n2'],f2['k2']
n2i=np.interp(w1,w2,n2); k2i=np.interp(w1,w2,k2)
m=(w1>380)&(w1<1000)
print('#1 vs #2 (380-1000): mean|dn|=%.4f mean|dk|=%.4f'%(np.mean(np.abs(n1[m]-n2i[m])),np.mean(np.abs(k1[m]-k2i[m]))))
for tw in [400,450,550,633,800,1000]:
    i=int(np.argmin(abs(w1-tw)))
    print('  %4dnm  #1 n=%.3f k=%.4f (res %.3f) | #2 n=%.3f k=%.4f'%(w1[i],n1[i],k1[i],r1[i],np.interp(w1[i],w2,n2),np.interp(w1[i],w2,k2)))
np.savez(str(OUT_DIR / 'ito_v2.npz'),w1=w1,n1=n1,k1=k1,r1=r1,d1=DZ1,dr1=DR1)
print('saved ito_v2.npz')

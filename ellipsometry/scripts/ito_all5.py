"""Fit all 5 sheets (#1,#2 vendor-A ITO; #3,#4 vendor-B IZO; #5 vendor-A ITO 2%O2).
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

Uniform pipeline: header-autodetect loader -> geometry scan (360-1080 score)
-> full chain extraction at best (d, rough)."""
import numpy as np, openpyxl, ellipsometry_fit as ef, csv
from scipy.optimize import least_squares
fp=str(DATA_DIR / 'se추출.xlsx')
z=np.load(str(OUT_DIR / 'izo_data.npz'))
ang0=np.array([65.0,70.0,75.0])
def tofloat(v):
    if isinstance(v,(int,float)): return float(v)
    try: return float(str(v).strip())
    except: return np.nan
wb=openpyxl.load_workbook(fp,data_only=True,read_only=True)
def read_sheet(name):
    ws=wb[name]
    allrows=[list(r) for r in ws.iter_rows(min_row=1,max_row=6,values_only=True)]
    hrow=None
    for i,r in enumerate(allrows):
        if any(v and 'Wavelength' in str(v) for v in r):
            hrow=i; break
    hdr=allrows[hrow]
    depc=[j for j,v in enumerate(hdr) if v and 'Depolarization' in str(v)]
    rows=[]
    for row in ws.iter_rows(min_row=hrow+2,values_only=True):
        rows.append([tofloat(v) for v in (list(row)+[None]*120)[:120]])
    A=np.array(rows); ok=np.isfinite(A[:,0]); A=A[ok]
    wl=A[:,0]; P=A[:,1:7:2]; D=A[:,2:7:2]
    dep=A[:,depc] if len(depc)==3 else np.zeros((len(wl),3))
    return wl,P,D,dep,len(depc)==3
SHEETS=['#1','#2','#3','#4','#5']
RAW={sh:read_sheet(sh) for sh in SHEETS}
wb.close()
for sh in SHEETS:
    wl,P,D,dep,hasdep=RAW[sh]
    print('%s: %d pts %.0f-%.0f depol=%s'%(sh,len(wl),wl[0],wl[-1],hasdep))

def prep(sh,lo,hi):
    wl_,P_,D_,dep_,_=RAW[sh]
    m=(wl_>=lo)&(wl_<=hi)
    wl=wl_[m]; P=P_[m]; D=D_[m]; dep=dep_[m]
    N_si=np.sqrt((np.interp(wl,z['si_wl'],z['si_e1'])+1j*np.interp(wl,z['si_wl'],z['si_e2'])).astype(complex))
    N_ox=np.sqrt((np.interp(wl,z['ox_wl'],z['ox_e1'])+1j*np.interp(wl,z['ox_wl'],z['ox_e2'])).astype(complex))
    pr=np.deg2rad(P); drr=np.deg2rad(D)
    return (wl,np.cos(2*pr),np.sin(2*pr)*np.cos(drr),np.sin(2*pr)*np.sin(drr),
            1.0/(1.0+(np.abs(dep)/3.0)**2),N_si,N_ox)

def chain(sh,DZ,DR,lo,hi,stride):
    wl,Nm,Cm,Sm,W,N_si,N_ox=prep(sh,lo,hi)
    il=np.arange(len(wl)-1,-1,-stride)
    prev=np.array([1.9,0.01]); out=[]
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
            s=np.clip(s,[1.0,0.0],[3.2,1.8])
            r=least_squares(f,s,bounds=([1.0,0.0],[3.2,1.8]),xtol=1e-9)
            if b is None or r.cost<b.cost: b=r
        prev=b.x.copy()
        out.append((wl[i],b.x[0],b.x[1],np.sqrt(2*b.cost/6)))
    out=np.array(out); o=np.argsort(out[:,0])
    return out[o,0],out[o,1],out[o,2],out[o,3]

RES={}
for sh in SHEETS:
    best=None
    for DZ in [42,44,46,48,50,52,54,56]:
        for DR in [2.0,5.0]:
            _,_,_,r=chain(sh,DZ,DR,360,1080,5)
            v=r.mean()
            if best is None or v<best[0]: best=(v,DZ,DR)
    v,DZ,DR=best
    w,n,k,r=chain(sh,DZ,DR,250,1080,1)
    RES[sh]=(w,n,k,r,DZ,DR)
    g=r<0.02
    print('%s BEST d=%.0f rough=%.0f scanres=%.4f | chain med=%.4f | n550=%.3f n633=%.3f n800=%.3f k633=%.4f k800=%.4f'%(
        sh,DZ,DR,v,np.median(r),np.interp(550,w,n),np.interp(633,w,n),np.interp(800,w,n),
        np.interp(633,w,k),np.interp(800,w,k)))

with open(str(OUT_DIR / 'ITO_IZO_all5_nk.csv'),'w',newline='') as f:
    w0=RES['#1'][0]
    wtr=csv.writer(f)
    hdr=['wl_nm']
    for sh in SHEETS: hdr+=['n_%s'%sh[1],'k_%s'%sh[1]]
    wtr.writerow(hdr)
    for i in range(len(w0)):
        row=['%.2f'%w0[i]]
        for sh in SHEETS:
            w,n,k,_,_,_=RES[sh]
            row+=['%.4f'%np.interp(w0[i],w,n),'%.4f'%np.interp(w0[i],w,k)]
        wtr.writerow(row)
np.savez(str(OUT_DIR / 'ito_izo_all5.npz'),
         **{('%s_%s'%(sh[1],t)):RES[sh][j] for sh in SHEETS for j,t in enumerate(['w','n','k','r'])},
         **{('%s_geo'%sh[1]):np.array([RES[sh][4],RES[sh][5]]) for sh in SHEETS})
print('saved ITO_IZO_all5_nk.csv / ito_izo_all5.npz')

"""Substrate-model sensitivity: re-extract nk of all 5 samples under
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

alternative Si/oxide references; quantify nk spread; T-check for ITO."""
import numpy as np, openpyxl, ellipsometry_fit as ef
from scipy.optimize import least_squares
fp=str(DATA_DIR / 'se추출.xlsx')
z=np.load(str(OUT_DIR / 'izo_data.npz'))
ang0=[65.0,70.0,75.0]
def tofloat(v):
    try: return float(str(v).strip())
    except: return np.nan
wb=openpyxl.load_workbook(fp,data_only=True,read_only=True)
def read_sheet(name):
    ws=wb[name]
    allrows=[list(r) for r in ws.iter_rows(min_row=1,max_row=6,values_only=True)]
    hrow=next(i for i,r in enumerate(allrows) if any(v and 'Wavelength' in str(v) for v in r))
    rows=[]
    for row in ws.iter_rows(min_row=hrow+2,values_only=True):
        rows.append([tofloat(v) for v in (list(row)+[None]*120)[:120]])
    A=np.array(rows); A=A[np.isfinite(A[:,0])]
    return A[:,0],A[:,1:7:2],A[:,2:7:2]
DATA={s:read_sheet('#'+s) for s in '12345'}
wb.close()

def subs(wl,variant):
    # Si
    if variant.startswith('aspnes'):
        N_si=ef.n_Si(wl)
    else:
        N_si=np.sqrt((np.interp(wl,z['si_wl'],z['si_e1'])+1j*np.interp(wl,z['si_wl'],z['si_e2'])).astype(complex))
    # oxide
    if 'ntve' in variant:
        N_ox=np.sqrt((np.interp(wl,z['ox_wl'],z['ox_e1'])+1j*np.interp(wl,z['ox_wl'],z['ox_e2'])).astype(complex))
        t_ox=3.0
    else:
        N_ox=ef.n_SiO2(wl)
        t_ox=float(variant.split('_')[-1])
    return N_si,N_ox,t_ox

def chain(s,variant,DZ,DR,lo,hi,stride):
    wl_,P_,D_=DATA[s]
    m=(wl_>=lo)&(wl_<=hi)
    wl=wl_[m]; P=P_[m]; D=D_[m]
    pr=np.deg2rad(P); drr=np.deg2rad(D)
    Nm=np.cos(2*pr); Cm=np.sin(2*pr)*np.cos(drr); Sm=np.sin(2*pr)*np.sin(drr)
    N_si,N_ox,t_ox=subs(wl,variant)
    il=np.arange(len(wl)-1,-1,-stride)
    prev=np.array([1.9,0.01]); out=[]
    for i in il:
        def f(x):
            Nf=np.array([complex(x[0],max(x[1],0.0))])
            na=np.ones(1,dtype=complex)
            layers=[na, ef.bruggeman_ema(Nf,na,0.5), Nf, N_ox[i:i+1], N_si[i:i+1]]
            dl=[DR,DZ,t_ox]
            NN=np.empty(3);CC=np.empty(3);SS=np.empty(3)
            for j,ang in enumerate(ang0):
                rp,rs=ef._tmm(wl[i:i+1],layers,dl,ang); rho=rp/rs
                t=abs(rho[0]); p2=2*np.arctan(t); an=np.angle(rho[0])
                NN[j]=np.cos(p2);CC[j]=np.sin(p2)*np.cos(an);SS[j]=np.sin(p2)*np.sin(an)
            cont=1.5*np.array([x[0]-prev[0],x[1]-prev[1]])
            return np.concatenate([NN-Nm[i],CC-Cm[i],SS-Sm[i],cont])
        b=None
        for sd in [prev,prev+[0.06,0],prev-[0.06,0]]:
            sd=np.clip(sd,[1.0,0.0],[3.2,1.8])
            r=least_squares(f,sd,bounds=([1.0,0.0],[3.2,1.8]),xtol=1e-9)
            if b is None or r.cost<b.cost: b=r
        prev=b.x.copy(); out.append((wl[i],b.x[0],b.x[1],np.sqrt(2*b.cost/6)))
    out=np.array(out); o=np.argsort(out[:,0])
    return out[o,0],out[o,1],out[o,2],out[o,3]

VARIANTS=['jaw_ntve','jaw_sio2_1','jaw_sio2_2','jaw_sio2_3','jaw_sio2_4','aspnes_sio2_2']
D0={'1':50,'2':50,'3':42,'4':42,'5':42}; R0={'1':5,'2':2,'3':2,'4':2,'5':2}
RES={}
for s in '12345':
    RES[s]={}
    for v in VARIANTS:
        best=None
        for dd in [D0[s]-2,D0[s],D0[s]+2]:
            _,_,_,r=chain(s,v,dd,R0[s],360,1080,6)
            m=r.mean()
            if best is None or m<best[0]: best=(m,dd)
        w,n,k,r=chain(s,v,best[1],R0[s],360,1080,2)
        RES[s][v]=(w,n,k,r,best[1])
        print('#%s %-13s d=%d res=%.4f  n550=%.3f n633=%.3f n800=%.3f k800=%.4f'%(
            s,v,best[1],best[0],np.interp(550,w,n),np.interp(633,w,n),np.interp(800,w,n),np.interp(800,w,k)))
    ns=[np.interp(633,RES[s][v][0],RES[s][v][1]) for v in VARIANTS]
    print('  -> #%s n633 spread over substrate models: %.3f - %.3f (Δ=%.3f)'%(s,min(ns),max(ns),max(ns)-min(ns)))
np.savez(str(OUT_DIR / 'substrate_sens.npz'),
         **{f'{s}_{v}_{t}':RES[s][v][j] for s in '12345' for v in VARIANTS for j,t in enumerate(['w','n','k','r'])})
print('saved substrate_sens.npz')

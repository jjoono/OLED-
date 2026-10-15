"""ITO 50nm on Si/NO — #1 & #2 (same batch). Staged NCS fit:
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

S1 window Cauchy -> S2 full-range Drude+TL(KK) -> S3 all free."""
import numpy as np, openpyxl, ellipsometry_fit as ef
from scipy.optimize import least_squares
fp=str(DATA_DIR / 'se추출.xlsx')
z=np.load(str(OUT_DIR / 'izo_data.npz'))
ang0=np.array([65.0,70.0,75.0]); NA=3

def tofloat(v):
    if isinstance(v,(int,float)): return float(v)
    try: return float(str(v).strip())
    except: return np.nan

wb=openpyxl.load_workbook(fp,data_only=True,read_only=True)
def read_sheet(name,mode):
    ws=wb[name]; rows=[]
    for row in ws.iter_rows(min_row=4,values_only=True):
        rows.append([tofloat(v) for v in (list(row)+[None]*120)[:120]])
    A=np.array(rows)
    ok=np.isfinite(A[:,0]); A=A[ok]
    wl=A[:,0]
    if mode=='grouped':
        P=np.column_stack([A[:,1],A[:,5],A[:,9]])
        D=np.column_stack([A[:,2],A[:,6],A[:,10]])
        dep=np.column_stack([A[:,97],A[:,103],A[:,109]])
    else:
        P=A[:,1:7:2]; D=A[:,2:7:2]
        dep=A[:,33:36]
    return wl,P,D,dep
wl1,P1,D1,dep1=read_sheet('#1','grouped')
wl2,P2,D2,dep2=read_sheet('Sheet2','interleaved')
wb.close()
print('#1: %d pts %.0f-%.0f | #2: %d pts %.0f-%.0f'%(len(wl1),wl1[0],wl1[-1],len(wl2),wl2[0],wl2[-1]))
for tag,wl_,dep_ in [('#1',wl1,dep1),('#2',wl2,dep2)]:
    for lo,hi in [(250,1100),(1100,1700)]:
        m=(wl_>=lo)&(wl_<hi)
        print('  %s depol %d-%d: mean|d|=%.2f max=%.2f'%(tag,lo,hi,np.nanmean(np.abs(dep_[m])),np.nanmax(np.abs(dep_[m]))))

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
def tl_eps2(E,A,E0,C,Egp):
    den=(E**2-E0**2)**2+C**2*E**2
    return np.where(E>Egp, A*E0*C*(E-Egp)**2/(den*np.maximum(E,1e-9)), 0.0)

def fit_sample(tag,wl_,P_,D_,dep_):
    m=(wl_>=250)&(wl_<=1080)
    wl=wl_[m]; P=P_[m]; D=D_[m]; dep=dep_[m]
    E_wl=1239.84193/wl
    N_si=np.sqrt((np.interp(wl,z['si_wl'],z['si_e1'])+1j*np.interp(wl,z['si_wl'],z['si_e2'])).astype(complex))
    N_ox=np.sqrt((np.interp(wl,z['ox_wl'],z['ox_e1'])+1j*np.interp(wl,z['ox_wl'],z['ox_e2'])).astype(complex))
    n_air=np.ones(len(wl),dtype=complex)
    pr=np.deg2rad(P); drr=np.deg2rad(D)
    Nm=np.cos(2*pr); Cm=np.sin(2*pr)*np.cos(drr); Sm=np.sin(2*pr)*np.sin(drr)
    W=1.0/(1.0+(np.abs(dep)/3.0)**2)
    def eps_ito(p):
        C0,Ad,gd,At,E0t,Ct,Egt=p
        e2i=np.maximum(tl_eps2(Eg,At,E0t,Ct,Egt),0)
        e1g=C0+kk(e2i)-Ad/(Eg**2+gd**2)
        e2g=e2i+Ad*gd/(Eg*(Eg**2+gd**2))
        e1=np.interp(E_wl,Eg,e1g); e2=np.interp(E_wl,Eg,np.maximum(e2g,0))
        Nf=np.sqrt((e1+1j*e2).astype(complex))
        return np.where(Nf.imag<0,np.conj(Nf),Nf)
    def fwd(Nf,dz,dr):
        layers=[n_air, ef.bruggeman_ema(Nf,n_air,0.5), Nf, N_ox, N_si]
        dl=[dr,dz,3.0]
        out=[]
        for ang in ang0:
            rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
            t=np.abs(rho); p2=2*np.arctan(t); an=np.angle(rho)
            out.append((np.cos(p2),np.sin(p2)*np.cos(an),np.sin(p2)*np.sin(an)))
        return out
    # S1: window Cauchy
    mw=(wl>=500)&(wl<=1050)
    def s1(x):
        dz,dr,A,B=x
        Nf=((A+B/wl**2)+0j).astype(complex)
        mm=fwd(Nf,dz,dr)
        r=[]
        for j,(a,b,c) in enumerate(mm):
            r+=[(W[:,j]*(a-Nm[:,j]))[mw],(W[:,j]*(b-Cm[:,j]))[mw],(W[:,j]*(c-Sm[:,j]))[mw]]
        return np.concatenate(r)
    b1=None
    for d0 in [40,50,60,70]:
        r=least_squares(s1,[d0,2.0,1.85,3e4],bounds=([25,0,1.5,0],[90,8,2.4,1.5e5]))
        if b1 is None or r.cost<b1.cost: b1=r
    dz1,dr1,A1,B1=b1.x
    n550=A1+B1/550**2
    print('%s [S1] d=%.2f rough=%.2f n550=%.3f n900=%.3f'%(tag,dz1,dr1,n550,A1+B1/900**2))
    # S2/S3
    def res_all(x):
        dz,dr=x[0],x[1]
        Nf=eps_ito(x[2:])
        mm=fwd(Nf,dz,dr)
        r=[]
        for j,(a,b,c) in enumerate(mm):
            r+=[W[:,j]*(a-Nm[:,j]),W[:,j]*(b-Cm[:,j]),W[:,j]*(c-Sm[:,j])]
        return np.concatenate(r)
    x0=[dz1,dr1, n550**2+0.3, 0.6,0.12, 40.,6.0,1.5,3.7]
    lo=[dz1-10,0.0, 1.5, 0.0,0.03, 0.0,4.5,0.3,3.0]
    hi=[dz1+10,8.0, 6.0, 8.0,1.0, 300.,8.5,4.0,4.3]
    r=least_squares(res_all,x0,bounds=(lo,hi),max_nfev=700)
    r=least_squares(res_all,r.x,bounds=(lo,hi),xtol=1e-12,max_nfev=700)
    x=r.x
    ncs=np.sqrt(2*r.cost/(3*NA*len(wl)))
    Nf=eps_ito(x[2:])
    print('%s [FINAL] d=%.2f rough=%.2f NCS-rms=%.4f'%(tag,x[0],x[1],ncs))
    print('   C0=%.3f Drude(A=%.3f g=%.3f) TL(A=%.1f E0=%.2f C=%.2f Eg=%.2f)'%tuple(x[2:]))
    for tw in [300,350,400,450,550,633,800,1000]:
        i=int(np.argmin(abs(wl-tw)))
        print('    %4dnm n=%.3f k=%.4f'%(wl[i],Nf[i].real,Nf[i].imag))
    return wl,x,Nf,ncs

wlA,xA,NA1,rA=fit_sample('#1',wl1,P1,D1,dep1)
wlB,xB,NB1,rB=fit_sample('#2',wl2,P2,D2,dep2)
print('\n== batch consistency ==')
print('d: %.2f vs %.2f nm | n550: %.3f vs %.3f | n633: %.3f vs %.3f'%(
      xA[0],xB[0],np.interp(550,wlA,NA1.real),np.interp(550,wlB,NB1.real),
      np.interp(633,wlA,NA1.real),np.interp(633,wlB,NB1.real)))
np.savez(str(OUT_DIR / 'ito_fit.npz'),wlA=wlA,xA=xA,nA=NA1.real,kA=NA1.imag,
         wlB=wlB,xB=xB,nB=NB1.real,kB=NB1.imag)
print('saved ito_fit.npz')

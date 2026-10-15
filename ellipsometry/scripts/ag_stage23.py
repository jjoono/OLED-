"""Stage 2 refined (Ag25, gamma bound, wl>=300) + Stage 3 (Ag3/Ag5 thin films:
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

Drude + island-plasmon Lorentzians, shared eps per nominal thickness)."""
import numpy as np, ellipsometry_fit as ef
from scipy.optimize import least_squares
d=np.load(str(OUT_DIR / 'aghatcn_data.npz'))
z=np.load(str(OUT_DIR / 'izo_data.npz'))
s1=np.load(str(OUT_DIR / 'hatcn_stage1b.npz'))
ang0=np.array([45.0,50.0,55.0,60.0,65.0]); NA=5
D_HATCN={'HATCN2':11.0,'HATCN4':11.8,'HATCN6':12.7,'HATCN30':float(s1['x30'][0])}

A0=d['HATCN2_Ag25']; wlfull=A0[np.isfinite(A0[:,0]),0]
mask=wlfull>=300.0
wl=wlfull[mask]
E_wl=1239.84193/wl
N_si=np.sqrt((np.interp(wl,z['si_wl'],z['si_e1'])+1j*np.interp(wl,z['si_wl'],z['si_e2'])).astype(complex))
N_ox=np.sqrt((np.interp(wl,z['ox_wl'],z['ox_e1'])+1j*np.interp(wl,z['ox_wl'],z['ox_e2'])).astype(complex))
n_air=np.ones(len(wl),dtype=complex)
N_hat=(np.interp(wl,s1['wl'],s1['n'])+1j*np.interp(wl,s1['wl'],s1['k'])).astype(complex)

def sheet_ncs(sh):
    A=d[sh]; ok=np.isfinite(A[:,0])
    A=A[ok][mask]
    P=A[:,1:11:2]; D_=A[:,2:11:2]; dep=A[:,45:50]
    pr=np.deg2rad(P); drr=np.deg2rad(D_)
    return (np.cos(2*pr),np.sin(2*pr)*np.cos(drr),np.sin(2*pr)*np.sin(drr),
            1.0/(1.0+(np.abs(dep)/3.0)**2))

def eps_metal(p):
    C0,Ad,gd,A1,E1,g1,A2,E2,g2=p
    eps=(C0+0j)*np.ones(len(E_wl),dtype=complex)
    eps-=Ad/(E_wl**2+1j*gd*E_wl)
    eps+=A1/(E1**2-E_wl**2-1j*g1*E_wl)
    eps+=A2/(E2**2-E_wl**2-1j*g2*E_wl)
    return eps

def model_sheet(pM,dM,drg,dHat):
    NM=np.sqrt(eps_metal(pM)); NM=np.where(NM.imag<0,np.conj(NM),NM)
    layers=[n_air, ef.bruggeman_ema(NM,n_air,0.5), NM, N_hat, N_ox, N_si]
    dl=[drg,dM,dHat,3.0]
    out=[]
    for ang in ang0:
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        t=np.abs(rho); p2=2*np.arctan(t); an=np.angle(rho)
        out.append((np.cos(p2),np.sin(p2)*np.cos(an),np.sin(p2)*np.sin(an)))
    return out

def global_fit(suffix,x0,lo,hi,dlo,dhi):
    SH=[f'HATCN{b}_{suffix}' for b in ['2','4','6','30']]
    DATA={sh:sheet_ncs(sh) for sh in SH}
    def resid(x):
        pM=x[:9]; res=[]
        for i,sh in enumerate(SH):
            dM=x[9+2*i]; drg=x[10+2*i]
            Nm,Cm,Sm,W=DATA[sh]
            mm=model_sheet(pM,dM,drg,D_HATCN[sh.split('_')[0]])
            for j,(a,b,c) in enumerate(mm):
                res+=[W[:,j]*(a-Nm[:,j]),W[:,j]*(b-Cm[:,j]),W[:,j]*(c-Sm[:,j])]
        return np.concatenate(res)
    LO=np.concatenate((lo,np.array(sum([[dlo,0.0] for _ in SH],[]))))
    HI=np.concatenate((hi,np.array(sum([[dhi,4.0] for _ in SH],[]))))
    X0=np.concatenate((x0,np.array(sum([[0.5*(dlo+dhi),1.0] for _ in SH],[]))))
    r=least_squares(resid,X0,bounds=(LO,HI),max_nfev=900)
    r=least_squares(resid,r.x,bounds=(LO,HI),xtol=1e-12,max_nfev=900)
    npts=DATA[SH[0]][0].shape[0]*len(SH)
    ncs=np.sqrt(2*r.cost/(3*NA*npts))
    x=r.x
    print('[%s] NCS-rms=%.4f'%(suffix,ncs))
    print('  eps: C0=%.2f Drude(A=%.1f g=%.3f) L1(A=%.1f E=%.2f g=%.2f) L2(A=%.1f E=%.2f g=%.2f)'%tuple(x[:9]))
    for i,sh in enumerate(SH):
        print('   %s: d=%.2fnm rough=%.2fnm'%(sh,x[9+2*i],x[10+2*i]))
    NM=np.sqrt(eps_metal(x[:9])); NM=np.where(NM.imag<0,np.conj(NM),NM)
    print('  nk: '+' '.join('%d:(%.3f,%.3f)'%(t,np.interp(t,wl,NM.real),np.interp(t,wl,NM.imag))
          for t in [320,400,500,633,800,1000,1400]))
    return x,NM,ncs

print('===== Stage 2 refined: Ag25 =====')
x25,N25,_=global_fit('Ag25',
    x0=np.array([3.7,75.,0.06,20.,4.3,0.8,30.,5.5,1.5]),
    lo=np.array([1.0,20.,0.02,0.0,3.9,0.1,0.0,4.8,0.2]),
    hi=np.array([9.0,150.,0.5,80.,4.7,3.0,200.,8.0,5.0]),
    dlo=15.,dhi=35.)
print('(lit Ag 633nm: n~0.14 k~3.99)')

print('\n===== Stage 3: Ag5 =====')
x5,N5,_=global_fit('Ag5',
    x0=np.array([3.0,20.,0.15,8.,2.9,0.8,15.,4.4,1.0]),
    lo=np.array([1.0,0.0,0.02,0.0,1.8,0.1,0.0,3.8,0.1]),
    hi=np.array([9.0,120.,1.0,60.,3.8,3.0,150.,6.5,5.0]),
    dlo=1.0,dhi=12.)

print('\n===== Stage 3: Ag3 =====')
x3,N3,_=global_fit('Ag3',
    x0=np.array([3.0,5.,0.25,10.,2.7,0.9,10.,4.4,1.2]),
    lo=np.array([1.0,0.0,0.02,0.0,1.8,0.1,0.0,3.8,0.1]),
    hi=np.array([9.0,120.,1.5,60.,3.8,3.0,150.,6.5,5.0]),
    dlo=0.5,dhi=9.)

np.savez(str(OUT_DIR / 'ag_stages.npz'),wl=wl,
         x25=x25,n25=N25.real,k25=N25.imag,
         x5=x5,n5=N5.real,k5=N5.imag,
         x3=x3,n3=N3.real,k3=N3.imag)
print('\nsaved ag_stages.npz')

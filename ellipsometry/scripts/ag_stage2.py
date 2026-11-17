"""Stage 2: global fit of 4x Ag25 sheets. Shared Ag eps (Drude + 2 Lorentz),
import os, pathlib
DATA_DIR = pathlib.Path(os.environ.get('ELLIPS_DATA', './data'))
OUT_DIR = pathlib.Path(os.environ.get('ELLIPS_OUT', './out'))
SCRATCH_DIR = pathlib.Path(os.environ.get('ELLIPS_SCRATCH', './out'))

per-sheet d_Ag and surface-rough EMA. HATCN underlayer fixed (stage1 eps, d from
arg or default 0.65x nominal). NCS residuals, depol weights, full range
(Ag opaque -> backside irrelevant)."""
import sys, numpy as np, ellipsometry_fit as ef
from scipy.optimize import least_squares
d=np.load(str(OUT_DIR / 'aghatcn_data.npz'))
z=np.load(str(OUT_DIR / 'izo_data.npz'))
s1=np.load(str(OUT_DIR / 'hatcn_stage1b.npz'))
ang0=np.array([45.0,50.0,55.0,60.0,65.0]); NA=5
SHEETS=['HATCN2_Ag25','HATCN4_Ag25','HATCN6_Ag25','HATCN30_Ag25']
D_HATCN=dict(zip(SHEETS,[11.0,11.8,12.7,float(s1['x30'][0])]))  # measured bare-sheet values (stage1b)

A0=d[SHEETS[0]]
wl=A0[np.isfinite(A0[:,0]),0]
E_wl=1239.84193/wl
N_si=np.sqrt((np.interp(wl,z['si_wl'],z['si_e1'])+1j*np.interp(wl,z['si_wl'],z['si_e2'])).astype(complex))
N_ox=np.sqrt((np.interp(wl,z['ox_wl'],z['ox_e1'])+1j*np.interp(wl,z['ox_wl'],z['ox_e2'])).astype(complex))
n_air=np.ones(len(wl),dtype=complex)
# HATCN eps from stage1b nk (interp onto this wl grid — same grid anyway)
N_hat=( np.interp(wl,s1['wl'],s1['n'])+1j*np.interp(wl,s1['wl'],s1['k']) ).astype(complex)

DATA={}
for sh in SHEETS:
    A=d[sh]; ok=np.isfinite(A[:,0])
    A=A[ok]
    P=A[:,1:11:2]; D_=A[:,2:11:2]; dep=A[:,45:50]
    pr=np.deg2rad(P); drr=np.deg2rad(D_)
    DATA[sh]=(np.cos(2*pr),np.sin(2*pr)*np.cos(drr),np.sin(2*pr)*np.sin(drr),
              1.0/(1.0+(np.abs(dep)/3.0)**2))

def eps_ag(p,E):
    C0,Ad,gd,A1,E1,g1,A2,E2,g2=p
    eps=(C0+0j)*np.ones(len(E),dtype=complex)
    eps-=Ad/(E**2+1j*gd*E)
    eps+=A1/(E1**2-E**2-1j*g1*E)
    eps+=A2/(E2**2-E**2-1j*g2*E)
    return eps

def model_sheet(pAg,dAg,drg,dHat):
    eps=eps_ag(pAg,E_wl)
    NAg=np.sqrt(eps); NAg=np.where(NAg.imag<0,np.conj(NAg),NAg)
    layers=[n_air, ef.bruggeman_ema(NAg,n_air,0.5), NAg, N_hat, N_ox, N_si]
    dl=[drg,dAg,dHat,3.0]
    out=[]
    for ang in ang0:
        rp,rs=ef._tmm(wl,layers,dl,ang); rho=rp/rs
        t=np.abs(rho); p2=2*np.arctan(t); an=np.angle(rho)
        out.append((np.cos(p2),np.sin(p2)*np.cos(an),np.sin(p2)*np.sin(an)))
    return out

def resid(x):
    pAg=x[:9]
    res=[]
    for i,sh in enumerate(SHEETS):
        dAg=x[9+2*i]; drg=x[10+2*i]
        Nm,Cm,Sm,W=DATA[sh]
        mm=model_sheet(pAg,dAg,drg,D_HATCN[sh])
        for j,(a,b,c) in enumerate(mm):
            res+=[W[:,j]*(a-Nm[:,j]),W[:,j]*(b-Cm[:,j]),W[:,j]*(c-Sm[:,j])]
    return np.concatenate(res)

# Ag literature seeds: eps_inf~3.7, Drude Ad~ (hbar wp)^2=9.2^2=84? for Ag wp~9.2eV -> Ad=wp^2~84, gd~0.02eV
# interband ~4.0eV (L1), plus ~ higher
x0=np.array([3.7, 75.0, 0.05, 20.0, 4.7, 0.6, 30.0, 5.5, 1.5]+
            sum([[25.0,1.0] for _ in SHEETS],[]))
lo=np.array([1.0, 20.0, 0.005, 0.0, 3.8, 0.05, 0.0, 4.8, 0.2]+
            sum([[15.0,0.0] for _ in SHEETS],[]))
hi=np.array([9.0, 150.0, 0.6, 80.0, 4.8, 3.0, 200.0, 8.0, 5.0]+
            sum([[35.0,4.0] for _ in SHEETS],[]))
r=least_squares(resid,x0,bounds=(lo,hi),max_nfev=800)
r=least_squares(resid,r.x,bounds=(lo,hi),xtol=1e-12,max_nfev=800)
x=r.x
npts=sum(DATA[sh][0].shape[0] for sh in SHEETS)
ncs=np.sqrt(2*r.cost/(3*NA*npts))
print('[Ag25 global] NCS-rms=%.4f'%ncs)
print('Ag eps: C0=%.2f Drude(A=%.1f g=%.4f) L1(A=%.1f E=%.2f g=%.2f) L2(A=%.1f E=%.2f g=%.2f)'%tuple(x[:9]))
for i,sh in enumerate(SHEETS):
    print('  %s: d_Ag=%.2fnm rough=%.2fnm'%(sh,x[9+2*i],x[10+2*i]))
eps=eps_ag(x[:9],E_wl)
NAg=np.sqrt(eps); NAg=np.where(NAg.imag<0,np.conj(NAg),NAg)
print('Ag nk: ',' '.join('%d:(%.3f,%.3f)'%(t,np.interp(t,wl,NAg.real),np.interp(t,wl,NAg.imag))
      for t in [350,400,500,633,800,1000,1400]))
print('(JAW/Palik Ag ref approx 633nm: n~0.14 k~3.99; 400nm: n~0.05 k~2.1)')
np.savez(str(OUT_DIR / 'ag_stage2.npz'),x=x,wl=wl,nAg=NAg.real,kAg=NAg.imag)
print('saved ag_stage2.npz')

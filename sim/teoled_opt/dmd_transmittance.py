"""Transparent-electrode check: absolute T/R/A of glass | ZnS d1 | HATCN 5 | Ag 8 (joint nk on HATCN) | ZnS d2 | air at normal
incidence (coherent film stack, incoherent 1.52 glass incl. its back face), bare film vs DMD optimised for T(550) and for
mean T(400-700).  "rel" = T / T(bare glass) (the "relative transmittance" convention, can exceed 1)."""
import sys, numpy as np
sys.path[:0]=['/home/user/OLED-/sim/teoled_opt']
import tmodel as T
from scipy.optimize import minimize
a=np.genfromtxt('/home/user/OLED-/sim/teoled_opt/nk_joint/Ag_HATCN_2-6_8nm_nk.csv',delimiter=',',names=True)
L=np.arange(400,701,5.)
AG=np.interp(L,a['wavelength_nm'],a['n_joint'])+1j*np.interp(L,a['wavelength_nm'],a['k_joint'])
ZN=np.array([T.ZNS[T.idx(l)] for l in L]); HT=np.array([T.HAT[T.idx(l)] for l in L]); NG=1.52
def rt(ns,ds,lam):           # Macleod characteristic matrix (N = n - ik convention)
    ns=[np.conj(complex(x)) for x in ns]
    M=np.eye(2,dtype=complex)
    for n,d in zip(ns[1:-1],ds):
        dl=2*np.pi*n*d/lam; M=M@np.array([[np.cos(dl),1j*np.sin(dl)/n],[1j*n*np.sin(dl),np.cos(dl)]])
    n0,ns_=ns[0],ns[-1]; B,C=M@np.array([1,ns_])
    r=(n0*B-C)/(n0*B+C); t=2*n0/(n0*B+C); return r,t
def TRA(d1,d2,dag=8.0,dhat=5.0):
    out=[]
    for i,l in enumerate(L):
        lay=[1.0,ZN[i],AG[i],HT[i],ZN[i],NG]; ds=[d2,dag,dhat,d1]
        lay=[x for x,dd in zip(lay,[1]+ds+[1]) if dd>0 or x in (1.0,NG)]; ds=[dd for dd in ds if dd>0]
        r,t=rt(lay,ds,l); Rf=abs(r)**2; Tc=NG.real*abs(t)**2
        rb,_=rt(lay[::-1],ds[::-1],l); Rb=abs(rb)**2      # coherent stack seen from glass
        Rg=((NG-1)/(NG+1))**2; Tg=1-Rg
        Ttot=Tc*Tg/(1-Rb*Rg); Rtot=Rf+Tc*Rg*Tc/(1-Rb*Rg)  # incoherent glass back face (Tc reused by reciprocity)
        out.append((Ttot,Rtot,1-Ttot-Rtot))
    return np.array(out)
glass=TRA(0,0,0,0); Tg=glass[:,0]
i550=list(L).index(550.)
for name,(d1,d2) in [('bare (glass/HATCN5/Ag8/air)',(0,0))]:
    o=TRA(d1,d2); print('%-34s T550 %.3f (rel %.3f)  A550 %.3f  Tavg %.3f'%(name,o[i550,0],o[i550,0]/Tg[i550],o[i550,2],o[:,0].mean()))
for tgt in ('550','avg'):
    f=(lambda x:-TRA(*np.abs(x))[i550,0]) if tgt=='550' else (lambda x:-TRA(*np.abs(x))[:,0].mean())
    best=min((minimize(f,x0,method='Nelder-Mead',options=dict(xatol=0.5,fatol=1e-5)) for x0 in ([30,30],[50,40],[20,60],[60,20])),key=lambda r:r.fun)
    d1,d2=np.abs(best.x); o=TRA(d1,d2)
    print('DMD opt %-4s ZnS %4.1f/HATCN5/Ag8/ZnS %4.1f: T550 %.3f (rel %.3f)  R550 %.3f A550 %.3f  Tavg %.3f (rel %.3f)'%(tgt,d1,d2,o[i550,0],o[i550,0]/Tg[i550],o[i550,1],o[i550,2],o[:,0].mean(),(o[:,0]/Tg).mean()))
print('bare glass T550 %.3f'%Tg[i550])

"""Fit of the measured 30-deg s-pol pixel reflectance (2026-10-27: 59/69/80/84/85/86/87 % at 400-700 nm) with single-parameter
absorption / degraded-Ag / beam-spill models (stack of refl_30s.py).  Result: no single uniform absorber reproduces the smooth
monotonic spectrum (best rms 4.8 %, kcs), spill alone 8.8 %."""
import sys, numpy as np
sys.path[:0] = ['/home/user/OLED-/sim/diag', '/home/user/OLED-/sim/fig3c']
import fit_ref_angles as F, cps2
from scipy.optimize import minimize_scalar, minimize
NS=1.52; A=30.0; TG=np.degrees(np.arcsin(np.sin(np.radians(A))/NS))
c,ct=np.cos(np.radians(A)),np.cos(np.radians(TG)); RF=((c-NS*ct)/(c+NS*ct))**2
L=np.array([400,450,500,550,600,650,700]); J=L-400; M=np.array([59,69,80,84,85,86,87])/100
D=F.D
print('organic k at 400/450/550:', [ (n, [round(float(np.imag(getattr(D,n)[j])),4) for j in (0,50,150)]) for n in ('B3o','EMo','TCo','TPo','HAT','AG')])
base=(172.5/185,160/210,1.0,12.7)
def Rs(j,kd=0,kh=0,kc=0,dag=0,geo=base):
    f,g,z,di=geo; S=F.stack(j,f,g,z,di,1.0,kd,kh,kc)
    if dag:
        S=F.stack(j,f,g,z,di,1.0,kd,kh,kc)
    k0=2*np.pi/S.lam; kx=NS*np.sin(np.radians(TG))*k0
    no=list(S.all_no); ne=list(S.all_ne)
    if dag:  # degrade Ag: add dag to real part of the 100 nm Ag (last layer of 'above' -> find by largest imag)
        i=int(np.argmax([np.imag(x) for x in no])); no[i]=no[i]+dag; ne[i]=ne[i]+dag
    rs,_=cps2.stack_rt(np.array(no),np.array(ne),S.all_d,np.array([kx]),k0,'s'); R=float(abs(rs[0])**2)
    return RF+(1-RF)**2*R/(1-RF*R)
def curve(**kw): return np.array([Rs(j,**kw) for j in J])
print('no extra', np.round(100*curve(),1))
for name,key,bd in (('kd','kd',(0,0.3)),('kh','kh',(0,1.5)),('kc','kc',(0,5)),('Ag dn','dag',(0,1.5))):
    r=minimize_scalar(lambda v: np.sum((curve(**{key:v})-M)**2),bounds=bd,method='bounded')
    print(name, round(r.x,4), 'rms %.2f%%'%(100*np.sqrt(r.fun/7)), np.round(100*curve(**{key:r.x}),1))
# spill: fraction s of beam on non-metal region approximated by bare glass+organics ~ (use R_out ~ 0.10 flat)
for Rout in (0.08,0.12):
    r=minimize_scalar(lambda s: np.sum(((1-s)*curve()+s*Rout-M)**2),bounds=(0,0.5),method='bounded')
    print('spill Rout',Rout, round(r.x,3),'rms %.2f%%'%(100*np.sqrt(r.fun/7)))

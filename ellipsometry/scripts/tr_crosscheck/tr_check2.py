"""Validate the incoherent-substrate T/R machinery on bare glass, then use the
measured T/R to determine d_Ag independently, with my SE n,k held fixed."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import numpy as np, json
from scipy.optimize import minimize_scalar
import ce_osc as osc, ellipsometry_fit as ef, ag_final as AG
import tr_check as TC

TR=_os.path.join(ELLIPS_OUT, r'TR260820')
V3=json.load(open('ag_v3_result.json'))
SEEDP=json.load(open('ag_seed_result.json'))
SP={k:[r for r in SEEDP[k] if abs(r['d_ox']-2.0)<1e-9][0]['p'][1:] for k in ('HATCN','MoOx')}
def seedN(p,wl):
    E=1239.841984/wl; e1,e2=osc.tauc_lorentz(E,p[1],p[2],p[3],p[4])
    N=np.sqrt((p[0]+e1+1j*e2).astype(complex)); return np.where(N.imag<0,-N,N)

def stack_TR(wl, layers, ds, ng):
    """coherent front stack on a thick transparent slab, incoherent back surface"""
    amb=np.ones(len(wl))
    rp,rs = ef._tmm(wl,[amb]+layers+[ng],ds,0.0)
    Rf=(np.abs(rp)**2+np.abs(rs)**2)/2
    tp,ts = TC._t_amp(wl,[amb]+layers+[ng],ds,0.0)
    Tf=(np.abs(tp)**2+np.abs(ts)**2)/2*ng.real
    rp2,rs2 = ef._tmm(wl,[ng]+layers[::-1]+[amb],ds[::-1],0.0)
    Rb=(np.abs(rp2)**2+np.abs(rs2)**2)/2
    tp2,ts2 = TC._t_amp(wl,[ng]+layers[::-1]+[amb],ds[::-1],0.0)
    Tb=(np.abs(tp2)**2+np.abs(ts2)**2)/2/ng.real
    Rg=((ng.real-1)/(ng.real+1))**2
    den=1-Rb*Rg
    return Tf*(1-Rg)/den, Rf + Tf*Tb*Rg/den

wl=np.arange(400.,801.,2.)
ng=np.full(len(wl),1.5230+0j)
T0,R0=stack_TR(wl,[],[],ng)
i=np.argmin(abs(wl-550))
print('bare glass check @550 nm:  T = %.2f %%  R = %.2f %%   (measured 91.66 %% / ~8.5 %%)'%(100*T0[i],100*R0[i]))

D=np.genfromtxt(TR+r'\ALL_SAMPLES_TRA.csv',delimiter=',',names=True)
o=np.argsort(D['wavelength_nm']); wm=D['wavelength_nm'][o]
MAP={4:'1-6',5:'1-7',6:'1-8',7:'2-5',8:'2-6',10:'2-7',12:'2-8'}
print('\nd_Ag from THEIR absolute T/R, with MY SE n,k fixed')
print('%-6s %-9s %-9s %-9s | %-8s %-8s'%('nom','d SE','d from T','d from T+R','res T %p','res R %p'))
for nom,sh in MAP.items():
    p=np.array(V3[sh]['p']); Ns=seedN(SP['HATCN'],wl); Na=AG.agN(p,wl)
    Tm=np.interp(wl,wm,D['HATCN5_Ag%d_T_pct'%nom][o])/100
    Rm=np.interp(wl,wm,D['HATCN5_Ag%d_Rcorr_pct'%nom][o])/100
    m=(wl>=430)&(wl<=780)
    def cost(d,use_r):
        T,R=stack_TR(wl,[Na,Ns],[d,5.0],ng)
        e=(T[m]-Tm[m])**2
        if use_r: e=e+(R[m]-Rm[m])**2
        return np.mean(e)
    dT=minimize_scalar(lambda d: cost(d,False),bounds=(2.,20.),method='bounded').x
    dTR=minimize_scalar(lambda d: cost(d,True),bounds=(2.,20.),method='bounded').x
    T,R=stack_TR(wl,[Na,Ns],[dTR,5.0],ng)
    print('%-6d %-9.2f %-9.2f %-9.2f | %-8.2f %-8.2f'%(
        nom,p[0],dT,dTR,100*np.mean(np.abs(T[m]-Tm[m])),100*np.mean(np.abs(R[m]-Rm[m]))))

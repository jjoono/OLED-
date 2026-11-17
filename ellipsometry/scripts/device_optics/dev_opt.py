"""Angle-resolved loss with MY SE n,k, in the realistic stack
organic n=1.8 / HATCN 7.5 nm / Ag d / capping n=2.1 (variable) / air.
Reports s, p and the Lambertian-weighted average, and scans (d_Ag, cap)."""
import numpy as np, json, ce_osc as osc, ellipsometry_fit as ef, ag_final as AG, tr_check as TC
V3=json.load(open('ag_v3_result.json')); SEEDP=json.load(open('ag_seed_result.json'))
SPp=[r for r in SEEDP['HATCN'] if abs(r['d_ox']-2.0)<1e-9][0]['p'][1:]
wl=np.array([550.0]); N_ORG=1.80
def seedN(p,w):
    E=1239.841984/w; e1,e2=osc.tauc_lorentz(E,p[1],p[2],p[3],p[4])
    N=np.sqrt((p[0]+e1+1j*e2).astype(complex)); return np.where(N.imag<0,-N,N)
NS=seedN(SPp,wl)
MAP={5:'1-7',7:'2-5',8:'2-6',10:'2-7',12:'2-8'}

def A_ang(N_ag,d_ag,dcap,th,pol,ncap=2.10):
    org=np.array([N_ORG+0j]); cap=np.array([ncap+0j]); air=np.array([1.0+0j])
    lay=[org,NS,N_ag,cap,air]; ds=[7.5,d_ag,dcap]
    rp,rs=ef._tmm(wl,lay,ds,np.rad2deg(th)); tp,ts=TC._t_amp(wl,lay,ds,np.rad2deg(th))
    s0=N_ORG*np.sin(th); qo=np.sqrt(N_ORG**2-s0**2+0j); qa=np.sqrt(1.0-s0**2+0j)
    T=0.0 if (np.real(qa)<=1e-9 or abs(np.imag(qa))>1e-9) else \
      (abs(tp[0])**2 if pol=='p' else abs(ts[0])**2)*np.real(qa)/np.real(qo)
    R=abs(rp[0])**2 if pol=='p' else abs(rs[0])**2
    return max(0.0,1-R-T)

th=np.deg2rad(np.arange(1.,89.,2.)); w=np.sin(th)*np.cos(th); w/=w.sum()
def summary(N,d,dcap,ncap=2.10):
    ap=np.array([A_ang(N,d,dcap,t,'p',ncap) for t in th])
    as_=np.array([A_ang(N,d,dcap,t,'s',ncap) for t in th])
    m=(ap+as_)/2
    return 100*np.sum(m*w),100*as_.max(),np.rad2deg(th[as_.argmax()]),100*ap.max(),100*m[0]

print('MY SE n,k, cap = 65 nm : where does >15%% come from?')
print('%-7s | %-9s %-9s %-11s %-9s %-9s'%('d_Ag','A(0 deg)','A avg','s peak','at deg','p peak'))
for nom,sh in MAP.items():
    N=AG.agN(np.array(V3[sh]['p']),wl)
    av,sp,ta,pp,a0=summary(N,V3[sh]['p'][0],65.0)
    print('%-7d | %-9.2f %-9.2f %-11.1f %-9.0f %-9.1f'%(nom,a0,av,sp,ta,pp))

print('\ncapping scan for the 8 nm film (angle-averaged is what matters)')
print('%-9s | %-9s %-9s %-9s'%('cap (nm)','A(0 deg)','A avg','s peak'))
N8=AG.agN(np.array(V3['2-6']['p']),wl); d8=V3['2-6']['p'][0]
bestc=None
for dc in (0,20,40,55,65,75,85,95,110,130,150):
    av,sp,ta,pp,a0=summary(N8,d8,float(dc))
    if bestc is None or av<bestc[0]: bestc=(av,dc)
    print('%-9d | %-9.2f %-9.2f %-9.1f'%(dc,a0,av,sp))
print('  -> best angle-averaged: cap %d nm (%.2f %%)'%(bestc[1],bestc[0]))

print('\ncapping INDEX scan at its own best thickness (8 nm Ag)')
print('%-8s | %-9s %-9s %-9s'%('n_cap','best cap','A avg','s peak'))
for nc in (1.5,1.8,2.1,2.4,2.7):
    best=None
    for dc in np.arange(0,200.,5.):
        av,sp,ta,pp,a0=summary(N8,d8,float(dc),nc)
        if best is None or av<best[0]: best=(av,dc,sp)
    print('%-8.1f | %-9.0f %-9.2f %-9.1f'%(nc,best[1],best[0],best[2]))

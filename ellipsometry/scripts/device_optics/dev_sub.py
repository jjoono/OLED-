"""Correct stack: organic n=1.8 / HATCN / Ag / capping / SUBSTRATE n=1.8 (semi-
infinite, its modes extracted by the MLA).  Because the substrate is index
matched to the organic there is no total internal reflection anywhere, so every
internal angle can leave - the picture is completely different from cap-to-air."""
import numpy as np, json, ce_osc as osc, ellipsometry_fit as ef, ag_final as AG, tr_check as TC
V3=json.load(open('ag_v3_result.json')); SEEDP=json.load(open('ag_seed_result.json'))
SPp=[r for r in SEEDP['HATCN'] if abs(r['d_ox']-2.0)<1e-9][0]['p'][1:]
wl=np.array([550.0]); N_ORG=1.80
def seedN(p,w):
    E=1239.841984/w; e1,e2=osc.tauc_lorentz(E,p[1],p[2],p[3],p[4])
    N=np.sqrt((p[0]+e1+1j*e2).astype(complex)); return np.where(N.imag<0,-N,N)
NS=seedN(SPp,wl); MAP={5:'1-7',7:'2-5',8:'2-6',10:'2-7',12:'2-8'}

def A_ang(N_ag,d_ag,dcap,th,pol,ncap=2.10,nout=1.80):
    org=np.array([N_ORG+0j]); cap=np.array([ncap+0j]); out=np.array([nout+0j])
    lay=[org,NS,N_ag,cap,out]; ds=[7.5,d_ag,dcap]
    rp,rs=ef._tmm(wl,lay,ds,np.rad2deg(th)); tp,ts=TC._t_amp(wl,lay,ds,np.rad2deg(th))
    s0=N_ORG*np.sin(th); qo=np.sqrt(N_ORG**2-s0**2+0j); qx=np.sqrt(nout**2-s0**2+0j)
    T=0.0 if (np.real(qx)<=1e-9 or abs(np.imag(qx))>1e-9) else \
      (abs(tp[0])**2 if pol=='p' else abs(ts[0])**2)*np.real(qx)/np.real(qo)
    R=abs(rp[0])**2 if pol=='p' else abs(rs[0])**2
    return max(0.0,1-R-T)

th=np.deg2rad(np.arange(1.,89.,2.)); w=np.sin(th)*np.cos(th); w/=w.sum()
def summ(N,d,dcap,ncap=2.10,nout=1.80):
    ap=np.array([A_ang(N,d,dcap,t,'p',ncap,nout) for t in th])
    as_=np.array([A_ang(N,d,dcap,t,'s',ncap,nout) for t in th])
    m=(ap+as_)/2
    return 100*np.sum(m*w),100*as_.max(),np.rad2deg(th[as_.argmax()]),100*ap.max(),100*m[0]

N8=AG.agN(np.array(V3['2-6']['p']),wl); d8=V3['2-6']['p'][0]
print('what changes when the top medium is a matched n=1.8 substrate, not air')
print('%-22s | %-9s %-9s %-10s %-9s'%('top medium','A(0 deg)','A avg','s peak','p peak'))
for lab,nout in [('air (my earlier calc)',1.0),('substrate n=1.8',1.80)]:
    av,sp,ta,pp,a0=summ(N8,d8,65.0,2.10,nout)
    print('%-22s | %-9.2f %-9.2f %-10.1f %-9.1f'%(lab,a0,av,sp,pp))

print('\nangle profile with the n=1.8 substrate (Ag 8 nm, cap 65 nm, n_cap 2.1)')
print('%-6s | %s'%('pol','  '.join('%5.0f'%t for t in (0,20,40,50,60,70,80))))
for pol in ('p','s'):
    print('%-6s |  %s'%(pol,' '.join('%5.1f'%(100*A_ang(N8,d8,65.,np.deg2rad(t),pol)) for t in (0,20,40,50,60,70,80))))

print('\nthickness series, matched substrate, cap 65 nm n=2.1')
print('%-7s | %-9s %-9s %-9s'%('d_Ag','A(0 deg)','A avg','s peak'))
for nom,sh in MAP.items():
    N=AG.agN(np.array(V3[sh]['p']),wl)
    av,sp,ta,pp,a0=summ(N,V3[sh]['p'][0],65.0)
    print('%-7d | %-9.2f %-9.2f %-9.1f'%(nom,a0,av,sp))

print('\ncapping now: does it still matter?  (Ag 8 nm, matched substrate)')
print('%-8s %-9s | %-9s %-9s'%('n_cap','cap nm','A avg','s peak'))
for nc in (1.8,2.1,2.4,2.7):
    best=None
    for dc in np.arange(0,200.,5.):
        av,sp,ta,pp,a0=summ(N8,d8,float(dc),nc)
        if best is None or av<best[0]: best=(av,dc,sp)
    print('%-8.1f %-9.0f | %-9.2f %-9.1f'%(nc,best[1],best[0],best[2]))

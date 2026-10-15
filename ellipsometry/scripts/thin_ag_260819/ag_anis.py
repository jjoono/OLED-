"""Three questions:
 (1) n,k are material constants, not angle-dependent - but an ISOTROPIC fit to
     oblique SE returns something different from what near-normal T/R sees if
     the film is uniaxial.  Test: invert (n,k) from EACH angle separately at
     fixed d.  Consistent angles => isotropic; a drift with angle => anisotropy.
 (2) how tightly is n pinned at 633 nm?  Map the residual over (n,k).
"""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef, ag_final as AG

V3=json.load(open('ag_v3_result.json')); SEEDP=json.load(open('ag_seed_result.json'))
SP={k:[r for r in SEEDP[k] if abs(r['d_ox']-2.0)<1e-9][0]['p'][1:] for k in ('HATCN','MoOx')}
DSEED={'HATCN':7.50,'MoOx':8.00}
def seedN(p,wl):
    E=1239.841984/wl; e1,e2=osc.tauc_lorentz(E,p[1],p[2],p[3],p[4])
    N=np.sqrt((p[0]+e1+1j*e2).astype(complex)); return np.where(N.imag<0,-N,N)

def per_angle(sh, wl_t):
    grp=L.SPEC[sh][0]; d_ag=V3[sh]['p'][0]; rough=V3[sh]['p'][1]
    wl,P,D,_=L.load(sh); i=np.argmin(abs(wl-wl_t))
    w1=wl[i:i+1]; ox=cf._mat('NTVE_JAW',w1); si=cf._mat('SI_JAW',w1); Ns=seedN(SP[grp],w1)
    out=[]
    for a,an in enumerate(L.ANG):
        meas=np.array(cf.ncs(P[i:i+1,a],D[i:i+1,a])).ravel()
        def r(v):
            Na=np.array([complex(v[0],max(v[1],0.0))]); Nr=ef.bruggeman_ema50(Na,np.ones_like(Na))
            rp,rs=ef._tmm(w1,[np.ones_like(Na),Nr,Na,Ns,ox,si],[rough,d_ag,DSEED[grp],2.0],an)
            rr=rp/rs; q,dl=np.arctan(np.abs(rr)),np.angle(rr)
            return np.array([np.cos(2*q)[0],(np.sin(2*q)*np.cos(dl))[0],(np.sin(2*q)*np.sin(dl))[0]])-meas
        b=least_squares(r,[0.15,4.0],bounds=([0.01,0.5],[3.0,8.0]),xtol=1e-13,ftol=1e-13)
        out.append((an,b.x[0],b.x[1],np.sqrt(2*b.cost/3)))
    return out

def nk_map(sh, wl_t):
    grp=L.SPEC[sh][0]; d_ag=V3[sh]['p'][0]; rough=V3[sh]['p'][1]
    wl,P,D,_=L.load(sh); i=np.argmin(abs(wl-wl_t))
    M=[cf.ncs(P[i:i+1,a],D[i:i+1,a]) for a in range(5)]
    ns=np.linspace(0.02,0.45,120); ks=np.linspace(3.6,4.8,90)
    NN,KK=np.meshgrid(ns,ks,indexing='ij'); Nf=(NN+1j*KK).ravel()
    Nr=ef.bruggeman_ema50(Nf,np.ones_like(Nf)); amb=np.ones_like(Nf)
    w=np.full(Nf.shape,wl[i]); ox=np.full(Nf.shape,cf._mat('NTVE_JAW',wl[i:i+1])[0])
    si=np.full(Nf.shape,cf._mat('SI_JAW',wl[i:i+1])[0]); Ns=np.full(Nf.shape,seedN(SP[grp],wl[i:i+1])[0])
    tot=np.zeros(Nf.shape)
    for a,an in enumerate(L.ANG):
        rp,rs=ef._tmm(w,[amb,Nr,Nf,Ns,ox,si],[rough,d_ag,DSEED[grp],2.0],an)
        rr=rp/rs; q,dl=np.arctan(np.abs(rr)),np.angle(rr)
        for c,v in enumerate([np.cos(2*q),np.sin(2*q)*np.cos(dl),np.sin(2*q)*np.sin(dl)]):
            tot+=(v-M[a][c][0])**2
    R=np.sqrt(tot/15).reshape(NN.shape)
    j=np.unravel_index(np.argmin(R),R.shape)
    return ns,ks,R,ns[j[0]],ks[j[1]],R[j]

for sh,lab in [('2-5','HATCN / Ag 7'),('1-7','HATCN / Ag 5')]:
    print('=== %s  (%s) ==='%(sh,lab))
    print('  (n,k) inverted from EACH angle separately, 633 nm:')
    for an,n,k,res in per_angle(sh,633.):
        print('     %2.0f deg   n=%.4f  k=%.4f   residual %.5f'%(an,n,k,res))
    ns,ks,R,n0,k0,r0=nk_map(sh,633.)
    print('  joint 5-angle best: n=%.4f k=%.4f  (residual %.5f)'%(n0,k0,r0))
    prof=R.min(axis=1)
    print('  residual vs n (k optimised at each n):')
    for t in (0.05,0.08,0.10,0.15,0.20,0.28,0.35):
        i=np.argmin(abs(ns-t)); print('     n=%.2f -> %.5f  (%.1fx the minimum)'%(ns[i],prof[i],prof[i]/r0))
    print()

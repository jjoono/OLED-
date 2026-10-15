import sys, numpy as np
sys.path[:0]=['/home/user/OLED-/sim/scatter']
import isl_model as I
from mla_plus_scatter import slab
from run_isl import top_layers, coh, MLA
w=np.cos(np.radians(I.TH))*np.sin(np.radians(I.TH)); w/=w.sum()
def isl_budget(tfe,nL,nm,S,alb,dE,dH,ag15=None):
    if ag15 is not None: I.AG15=ag15
    BT,BR=MLA[nL]; tl=top_layers(tfe,dH); bel=[I.iso(I.BEB,dE),I.iso(I.IZO,50.)]
    sD,sU,_,_=I.emission(bel,tl,nm,nL)
    st=coh(bel[::-1])+[(I.EML[0],25.)]+coh(tl)
    Rd,TDU=I.transfer(st,nm,nL); Ru,TUD=I.transfer(st[::-1],nL,nm); Mg=I.mirror_R(nm)
    Ts,Rs=(np.eye(90),np.zeros((90,90))) if S==0 else slab(S,0.6,N=6000,seed=3,albedo=alb)
    # explicit: downward flux at slab top b, inside slab with mirror
    X=np.linalg.solve(np.eye(90)-Rs*Mg[None,:],Ts)     # flux arriving at mirror per unit b
    Reff=Rs+Ts@(Mg[:,None]*X)
    Iden=np.eye(90)
    A=np.block([[Iden-Ru[:,None]*BR,-TDU@Reff],[-TUD@BR,Iden-Rd[:,None]*Reff]])
    x=np.linalg.solve(A,np.concatenate([sU,sD])); a,b=x[:90],x[90:]
    c=Reff@b; down_U=BR@a
    esc=BT@a
    abs_stack_below=((1-Rd-TDU.sum(0))*c).sum(); abs_stack_above=((1-Ru-TUD.sum(0))*down_U).sum()
    abs_mirror=((1-Mg)*(X@b)).sum(); abs_scat=b.sum()-Reff.sum(0)@b-abs_mirror
    born=sU.sum()+sD.sum()
    print(f"  emitted up {sU.sum():.3f} down {sD.sum():.3f} lost-at-birth {1-born:.3f} | escaped {esc:.3f} | abs: stack(from ISL) {abs_stack_below:.3f}, stack(from lens) {abs_stack_above:.3f}, Ag mirror {abs_mirror:.3f}, scatterer {abs_scat:.3f} | check {1-born+esc+abs_stack_below+abs_stack_above+abs_mirror+abs_scat:.3f}")
    print(f"  Lambertian per-hit: stack-from-ISL A={w@(1-Rd-TDU.sum(0)):.3f} T={w@TDU.sum(0):.3f}; Ag mirror in n{nm} A={w@(1-Mg):.3f}; stack-from-lens A'={w@(1-Ru-TUD.sum(0)):.3f}")
print('lens1.7, ISL n1.85 S3 lossless, dE80 dH120 (measured Ag15 on HATCN, n=0.146)'); isl_budget(False,1.7,1.85,3,1.0,80,120)
print('same, S10 albedo 0.995'); isl_budget(False,1.7,1.85,10,0.995,80,120)
print('same, S3, ideal top Ag (McPeak n=0.044)'); isl_budget(False,1.7,1.85,3,1.0,80,120,ag15=I.AG)

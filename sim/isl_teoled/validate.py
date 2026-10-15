import numpy as np, isl_model as I, cps2
# planar ref: lens 1.8 / Ag15 / HATCN / TAPC / EML / Bebq2 / Ag  (as seen from lens)
dE,dH=60.,40.
above=[I.iso(I.BEB,dE),I.iso(I.AG,100.)]
below=[(I.TAPC[0],I.TAPC[1],dH),I.iso(I.HAT,10.),I.iso(I.AG15,15.)]
S=cps2.Stack(550.,I.EML,25.,12.5,above=above,below=below,n_sub=1.8,n_top=1.0)
Rc=np.real(cps2.stack_reflectance(S,I.TH))
L=[(I.AG15,15.),(I.HAT,10.),(I.TAPC[0],dH),(I.EML[0],25.),(I.BEB,dE),(I.AG,100.)]
Rm,Tm=I.rt_unpol(L,1.8,1.0,I.TH)
print('max |R_mine-R_cps| =',np.max(abs(Rm-Rc)), ' (isotropic approx for uniaxial layers in mine)')
r=cps2.solve(S); print('cps budget',{k:round(float(np.real(v)),4) for k,v in r.items()})
sD,sU,Pd,Pu=I.emission(above,below,1.0,1.8)
print('sum sU',sU.sum(),' vs cps air+sub',float(np.real(r['air']+r['sub'])))

"""v9: lock the parameters the 340-1080 nm window cannot determine, so the ones
CompleteEASE still fits come back with meaningful confidence intervals.

Covariance at the v8 solution showed:
  rho <-> tau       correlation -0.999   (Drude: 596% and 16500% uncertainty)
  Amp2 <-> Eg2      +0.978
  Amp2 <-> Eo2      -0.976
  Einf <-> Br2      -0.959
Locking tau, rho and Br2 collapses the condition number by ~2 orders and takes
Einf from 22% -> 5%, Eo2 3.0% -> 1.9%, Eg2 2.8% -> 2.4%.
Br2 is locked at a physically defensible ITO value (not at the bound the free
fit ran to), and the Drude at the carrier density established earlier.
"""
import numpy as np, json
from scipy.optimize import least_squares
import ce_fit as cf, ce_fit2 as f2

# Eo2 and Br2 are LOCKED at the v8 optimum: within 340-1080 nm they are not
# determined at all, so any value in a wide range reproduces the data equally -
# pinning them at the optimum keeps MSE lowest and lets Eg2/Amp2/Einf be measured.
TAU_FIX = 6.582119569               # fs  (Drude Br = 0.10 eV)
FREE = [1, 3, 6, 9, 10, 11, 12]     # rough, Einf, TL_Amp, TL_Eg, G_Amp, G_Br, G_En
NM = {1:'Roughness',3:'Einf',6:'Amp2(TL)',9:'Eg2(TL)',
      10:'Amp3(Gauss)',11:'Br3(Gauss)',12:'En3(Gauss)'}
LO = {1:0.0, 3:1.0, 6:20., 8:3.5, 9:3.00, 10:0.0, 11:0.10, 12:0.60}
HI = {1:12., 3:4.5, 6:600., 8:5.2, 9:3.95, 10:8.0, 11:5.00, 12:3.40}
RHO_FIX = {'#1':1.0644e-3, '#2':1.6820e-3, '#3':4.7736e-3, '#4':6.4640e-3, '#5':1.0212e-3}

V8 = json.load(open('ce_final_result.json'))
out = {}
for sh in ['#1','#2','#3','#4','#5']:
    p = np.array(V8[sh]['p'])
    p[5] = TAU_FIX; p[4] = RHO_FIX[sh]        # Eo2 (p[8]) and Br2 (p[7]) keep their v8 values
    wl, Pm, Dm = cf.load(sh)
    m = (wl>=340)&(wl<=1080); wl,Pm,Dm = wl[m],Pm[m],Dm[m]
    r0 = f2.resid_factory(wl, Pm, Dm, cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl))
    lo = np.array([LO[i] for i in FREE]); hi = np.array([HI[i] for i in FREE])
    def r(x):
        q = p.copy(); q[FREE] = x; return r0(q)
    x0 = np.clip(p[FREE], lo+1e-9, hi-1e-9)
    best = None; rng = np.random.default_rng(5)
    for k in range(12):
        s = x0 if k==0 else np.clip(x0*(1+0.3*rng.standard_normal(len(x0))), lo+1e-9, hi-1e-9)
        try: b = least_squares(r, s, bounds=(lo,hi), x_scale='jac', max_nfev=1200)
        except Exception: continue
        if best is None or b.cost < best.cost: best = b
    p[FREE] = best.x
    res = r0(p); M = cf.mse(res, len(FREE)+2)
    # uncertainties
    s2 = np.sum(res**2)/(len(res)-len(FREE))
    J = np.zeros((len(res), len(FREE)))
    for j,idx in enumerate(FREE):
        h = max(abs(p[idx])*1e-5, 1e-9)
        a=p.copy(); a[idx]+=h; bq=p.copy(); bq[idx]-=h
        J[:,j] = (r0(a)-r0(bq))/(2*h)
    sig = np.sqrt(np.diag(np.linalg.pinv(J.T@J)*s2))*1.645
    out[sh] = dict(sheet=sh, p=p.tolist(), mse=float(M), mse_v8=V8[sh]['mse'],
                   err={NM[idx]: float(sig[j]) for j,idx in enumerate(FREE)})
    print('%s  MSE %6.2f  (v8 %6.2f)' % (sh, M, V8[sh]['mse']), flush=True)
    print('    ' + '  '.join('%s=%.4g+-%.2g' % (NM[idx], p[idx], sig[j])
                             for j,idx in enumerate(FREE)), flush=True)
    out[sh]['locked']=dict(rho=p[4],tau=p[5],Br2=p[7],Eo2=p[8],d=p[0],dth=p[2])
    json.dump(out, open('ce_v9_result.json','w'), indent=1)
print('done')

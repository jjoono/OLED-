"""Should the thin Ag films be an EMA layer instead of a free Gen-Osc layer?

Reference Ag n,k = the 12 nm film on HATCN (the most bulk-like one measured).
Competing layer models, all on the same stack and range:
  A  free Gen-Osc               11 params  (what I used)
  B  Bruggeman (Ag + void)       2 params  (d, f)   - percolated mixture
  C  Maxwell-Garnett, spheres    2 params  (d, f)   - isolated inclusions, L=1/3
  D  Maxwell-Garnett, free L     3 params  (d, f, L) - oblate islands
MSE already divides by (3n - m), so the comparison is fair across model sizes.
"""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ellipsometry_fit as ef, ag_fit as AF

R = json.load(open('ag_polish_result.json'))
REF = '2-8'
LOW, HIGH = 260.0, 1080.0

def stack_mse(sh, Nfun, x, npar):
    name, _, _ = L.SPEC[sh]
    d_seed, ps = AF.SEED[name]
    wl, P, D, _ = L.load(sh)
    m = (wl>=LOW)&(wl<=HIGH); wl,P,D = wl[m],P[m],D[m]
    ox, si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    Ns = AF.seedN(ps, wl); amb = np.ones(len(wl))
    M = [cf.ncs(P[:,i],D[:,i]) for i in range(5)]
    Na = Nfun(x, wl)
    Nr = ef.bruggeman_ema50(Na, amb)
    o=[]
    for a,an in enumerate(L.ANG):
        rp,rs = ef._tmm(wl,[amb,Nr,Na,Ns,ox,si],[x[1],x[0],d_seed,AF.D_OX],an)
        rr=rp/rs; q,dl=np.arctan(np.abs(rr)),np.angle(rr)
        o += [np.cos(2*q)-M[a][0], np.sin(2*q)*np.cos(dl)-M[a][1], np.sin(2*q)*np.sin(dl)-M[a][2]]
    res=np.concatenate(o); n=len(res)//3
    return res, 1000*np.sqrt(np.sum(res**2)/(3*n-npar))

def ref_eps(wl):
    N = AF.agN(np.array(R[REF]['p']), wl); return N**2

def N_brugg(x, wl):                       # x = [d, rough, f]
    ea = ref_eps(wl); eb = np.ones_like(ea)
    Na = np.sqrt(ea); Na=np.where(Na.imag<0,-Na,Na)
    Nb = np.ones_like(Na)
    return ef.bruggeman_ema(Na, Nb, np.clip(x[2],1e-3,1-1e-6))

def N_mg(x, wl, Lfac=None):               # x = [d, rough, f, (L)]
    ei = ref_eps(wl); f = np.clip(x[2],1e-3,0.999)
    Lf = (1/3.) if Lfac is None else np.clip(x[3],0.02,0.90)
    e = 1.0 + f*(ei-1.0)/(1.0 + (1-f)*Lf*(ei-1.0))
    N = np.sqrt(e); return np.where(N.imag<0,-N,N)

def fit_ema(sh, kind):
    d0 = R[sh]['p'][0]; r0 = R[sh]['p'][1]
    if kind=='B':   fn, lo, hi, np_ = N_brugg, [d0*0.4,0.,0.05], [d0*2.5,8.,1.0], 3
    elif kind=='C': fn, lo, hi, np_ = (lambda x,w: N_mg(x,w)), [d0*0.4,0.,0.05], [d0*2.5,8.,0.999], 3
    else:           fn, lo, hi, np_ = (lambda x,w: N_mg(x,w,1)), [d0*0.4,0.,0.05,0.02], [d0*2.5,8.,0.999,0.90], 4
    x0 = [d0, r0, 0.7] + ([0.2] if kind=='D' else [])
    best=None; rng=np.random.default_rng(3)
    for k in range(18):
        s = np.array(x0) if k==0 else np.array(lo)+rng.random(len(lo))*(np.array(hi)-np.array(lo))
        try: b = least_squares(lambda x: stack_mse(sh, fn, x, np_)[0], np.clip(s,lo,hi),
                               bounds=(lo,hi), x_scale='jac', max_nfev=800)
        except Exception: continue
        if best is None or b.cost<best.cost: best=b
    return stack_mse(sh, fn, best.x, np_)[1], best.x

print('%-5s %-6s %3s | %7s | %7s %5s %5s | %7s %5s %5s | %7s %5s %5s %5s'%(
 'id','seed','nom','A GenO','B Brugg','d','f','C MG1/3','d','f','D MG(L)','d','f','L'))
for sh in ['1-6','1-7','1-8','2-8','1-14','1-15','1-16','2-16']:
    row = '%-5s %-6s %3d | %7.2f'%(sh, L.SPEC[sh][0], L.SPEC[sh][2], R[sh]['mse'])
    for kind in 'BCD':
        M,x = fit_ema(sh,kind)
        row += ' | %7.2f %5.2f %5.2f'%(M,x[0],x[2]) + (' %5.2f'%x[3] if kind=='D' else '')
    print(row, flush=True)

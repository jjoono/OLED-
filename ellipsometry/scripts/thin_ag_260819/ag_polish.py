"""Local-minimum check: restart every sample from every OTHER sample's solution
(thickness kept at its own value) and keep whichever lands lowest."""
import numpy as np, json
from scipy.optimize import least_squares
import ag_load as L, ce_fit as cf, ellipsometry_fit as ef, ag_fit as AF

R = json.load(open('ag_fit_result.json'))
SH = [s for s in L.ORDER if L.SPEC[s][2] > 0]

def resid(sh):
    name, _, _ = L.SPEC[sh]
    d_seed, ps = AF.SEED[name]
    wl, P, D, _ = L.load(sh)
    m = (wl>=AF.LOW)&(wl<=AF.HIGH); wl,P,D = wl[m],P[m],D[m]
    ox, si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    Ns = AF.seedN(ps, wl)
    M = [cf.ncs(P[:,i],D[:,i]) for i in range(5)]; amb = np.ones(len(wl))
    def r(p):
        Na = AF.agN(p, wl); Nr = ef.bruggeman_ema50(Na, amb); o=[]
        for a,an in enumerate(L.ANG):
            rp,rs = ef._tmm(wl,[amb,Nr,Na,Ns,ox,si],[p[1],p[0],d_seed,AF.D_OX],an)
            rr=rp/rs; q,dl=np.arctan(np.abs(rr)),np.angle(rr)
            o += [np.cos(2*q)-M[a][0], np.sin(2*q)*np.cos(dl)-M[a][1],
                  np.sin(2*q)*np.sin(dl)-M[a][2]]
        return np.concatenate(o)
    return r, len(M[0][0])

out = {}
wp = np.array([450.,550.,633.,800.,1000.])
for sh in SH:
    r, nwl = resid(sh)
    best = None
    starts = [np.array(R[sh]['p'])]
    for o in SH:
        if o == sh: continue
        q = np.array(R[o]['p']).copy(); q[0] = R[sh]['p'][0]; q[1] = R[sh]['p'][1]
        starts.append(q)
    for s in starts:
        s = np.clip(s, AF.LO+1e-9, AF.HI-1e-9)
        try: b = least_squares(r, s, bounds=(AF.LO, AF.HI), x_scale='jac', max_nfev=1500)
        except Exception: continue
        if best is None or b.cost < best.cost: best = b
    M = 1000*np.sqrt(np.sum(best.fun**2)/(3*nwl*5 - len(AF.LO)))
    N = AF.agN(best.x, wp)
    out[sh] = dict(sheet=sh, seed=L.SPEC[sh][0], ag_nom=L.SPEC[sh][2], mse=float(M),
                   p=[float(v) for v in best.x], n=[float(v) for v in N.real],
                   k=[float(v) for v in N.imag])
    old = R[sh]
    flag = '  <-- improved' if M < old['mse']-0.05 else ''
    print('%-5s %-6s nom=%2d  MSE %6.2f -> %6.2f   n633 %.3f -> %.3f   tau %.1f -> %.1f%s'
          %(sh, L.SPEC[sh][0], L.SPEC[sh][2], old['mse'], M, old['n'][2], N.real[2],
            old['p'][4], best.x[4], flag), flush=True)
    json.dump(out, open('ag_polish_result.json','w'), indent=1)
print('done')

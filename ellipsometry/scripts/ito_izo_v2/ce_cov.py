"""Parameter uncertainties and correlations at the v8 solution.
Large error bars on Eo2/Eg2/Br2 are not noise - they are the fingerprint of a
degenerate direction, so quantify it and find what to lock."""
import numpy as np, json
import ce_fit as cf, ce_fit2 as f2, ellipsometry_fit as ef

R = json.load(open('ce_final_result.json'))
LOW, HIGH = 340.0, 1080.0
FREE = [1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]        # rough,Einf,rho,tau,TL*,G*  (d,dth locked)
NM = {1:'rough',3:'Einf',4:'rho(Drude)',5:'tau(Drude)',6:'Amp2(TL)',7:'Br2(TL)',
      8:'Eo2(TL)',9:'Eg2(TL)',10:'Amp3(G)',11:'Br3(G)',12:'En3(G)'}

def analyse(sh):
    p = np.array(R[sh]['p'])
    wl, Pm, Dm = cf.load(sh)
    m = (wl >= LOW) & (wl <= HIGH); wl, Pm, Dm = wl[m], Pm[m], Dm[m]
    ox, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
    r0 = f2.resid_factory(wl, Pm, Dm, ox, si)
    res0 = r0(p); n = len(res0); k = len(FREE)
    s2 = np.sum(res0**2) / (n - k)
    J = np.zeros((n, k))
    for j, idx in enumerate(FREE):
        h = max(abs(p[idx]) * 1e-5, 1e-9)
        pp = p.copy(); pp[idx] += h
        pm2 = p.copy(); pm2[idx] -= h
        J[:, j] = (r0(pp) - r0(pm2)) / (2 * h)
    JTJ = J.T @ J
    C = np.linalg.pinv(JTJ) * s2
    sig = np.sqrt(np.diag(C))
    corr = C / np.outer(sig, sig)
    ev = np.linalg.eigvalsh(JTJ)
    print('=== %s   MSE=%.1f ===' % (sh, R[sh]['mse']))
    print('  parameter        value      90%% conf      rel.')
    for j, idx in enumerate(FREE):
        rel = 100 * 1.645 * sig[j] / abs(p[idx]) if p[idx] else np.inf
        flag = '   <-- UNIDENTIFIABLE' if rel > 60 else ('   <-- weak' if rel > 20 else '')
        print('  %-14s %10.4g  +-%9.3g  %7.1f%%%s' % (NM[idx], p[idx], 1.645*sig[j], rel, flag))
    print('  condition number of J^T J = %.2e  (>1e8 => degenerate directions)' % (ev[-1]/max(ev[0], 1e-30)))
    print('  strongest correlations:')
    pr = [(abs(corr[a, b]), NM[FREE[a]], NM[FREE[b]], corr[a, b])
          for a in range(k) for b in range(a+1, k)]
    pr.sort(reverse=True)
    for v, x, y, c in pr[:6]:
        print('     %-14s %-14s  %+.4f' % (x, y, c))
    # which direction is flattest?
    w, V = np.linalg.eigh(JTJ)
    v0 = V[:, 0]
    print('  flattest direction (eigval %.2e):' % w[0],
          ', '.join('%s %+.2f' % (NM[FREE[i]], v0[i]) for i in np.argsort(-abs(v0))[:4]))

for sh in ['#1', '#3']:
    analyse(sh); print()

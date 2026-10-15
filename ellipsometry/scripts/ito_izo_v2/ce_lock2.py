import numpy as np, json, ce_fit as cf, ce_fit2 as f2, ellipsometry_fit as ef

# ---- (1) is my TMM's ABSOLUTE reflectance trustworthy?  bare Si, analytic check
print('TMM absolute-reflectance check: bare Si (SI_JAW) at 65 deg')
wl = np.array([250., 400., 633., 800.])
si = cf._mat('SI_JAW', wl); amb = np.ones_like(si)
rp, rs = ef._tmm(wl, [amb, si], [], 65.0)
th = np.deg2rad(65.0)
c1 = np.cos(th); c2 = np.sqrt(1 - (np.sin(th)/si)**2)
fr_s = (c1 - si*c2)/(c1 + si*c2)
fr_p = (si*c1 - c2)/(si*c1 + c2)
for i, w in enumerate(wl):
    print('  %4.0f nm  |rs| tmm %.4f analytic %.4f   |rp| tmm %.4f analytic %.4f'
          % (w, abs(rs[i]), abs(fr_s[i]), abs(rp[i]), abs(fr_p[i])))

# ---- (2) which parameters to LOCK so the free ones get meaningful error bars
R = json.load(open('ce_final_result.json'))
NM = {1:'rough',3:'Einf',4:'rho',5:'tau',6:'Amp2',7:'Br2',8:'Eo2',9:'Eg2',10:'Amp3',11:'Br3',12:'En3'}
SCHEMES = [('A none',            [1,3,4,5,6,7,8,9,10,11,12]),
           ('B lock tau',        [1,3,4,  6,7,8,9,10,11,12]),
           ('C lock tau,Br2',    [1,3,4,  6,  8,9,10,11,12]),
           ('D lock rho,tau,Br2',[1,3,    6,  8,9,10,11,12])]

def unc(sh, free):
    p = np.array(R[sh]['p'])
    wl, Pm, Dm = cf.load(sh)
    m = (wl >= 340) & (wl <= 1080); wl, Pm, Dm = wl[m], Pm[m], Dm[m]
    r0 = f2.resid_factory(wl, Pm, Dm, cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl))
    res0 = r0(p); s2 = np.sum(res0**2)/(len(res0)-len(free))
    J = np.zeros((len(res0), len(free)))
    for j, idx in enumerate(free):
        h = max(abs(p[idx])*1e-5, 1e-9)
        a = p.copy(); a[idx] += h; b = p.copy(); b[idx] -= h
        J[:, j] = (r0(a) - r0(b))/(2*h)
    JTJ = J.T @ J
    sig = np.sqrt(np.diag(np.linalg.pinv(JTJ)*s2))
    ev = np.linalg.eigvalsh(JTJ)
    return {NM[idx]: 100*1.645*sig[j]/max(abs(p[idx]), 1e-12) for j, idx in enumerate(free)}, ev[-1]/max(ev[0],1e-30)

for sh in ['#1', '#3']:
    print('\n=== %s : 90%% confidence interval as %% of value ===' % sh)
    hdr = ['rough','Einf','rho','tau','Amp2','Br2','Eo2','Eg2','Amp3','Br3','En3']
    print('  %-20s' % 'scheme' + ''.join('%8s' % h for h in hdr) + '   cond')
    for tag, free in SCHEMES:
        u, cond = unc(sh, free)
        print('  %-20s' % tag + ''.join(('%7.1f ' % u[h]) if h in u else '   fix  ' for h in hdr)
              + '  %.1e' % cond)

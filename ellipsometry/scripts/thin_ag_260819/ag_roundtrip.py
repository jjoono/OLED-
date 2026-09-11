"""Round-trip check: read the generated .mod back, rebuild the stack from what
the FILE says, and score it against the measured data.  If the MSE matches the
fit, the file really encodes the model."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import re, os, json, numpy as np
import ag_load as L, ce_fit as cf, ce_osc as osc, ellipsometry_fit as ef

D = _os.path.join(ELLIPS_OUT, r'mod_ag260819_v2')
def parse(path):
    t = open(path,'rb').read().decode('utf-8-sig')
    rough = float(re.search(r"^\t([\d.eE+-]+)\t[TF]\t0\.0\t150\.0\tF\t'Roughness'", t, re.M).group(1))/10
    dox   = float(re.search(r"^\t\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tT\t'Native Oxide'", t, re.M).group(1))/10
    dth   = float(re.search(r"^\t([-\d.eE+]+)\t[TF]\t-2\.0\t2\.0\tF\t'Angle Offset'", t, re.M).group(1))
    thick = [float(x)/10 for x in re.findall(r"^\t\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'Thickness # \d+'", t, re.M)]
    layers = []
    for fb in re.findall(r"start_Gen-Osc Fit Parms\r\n(.*?)\r\n\t\tend_Gen-Osc Fit Parms", t, re.S):
        Ls = [x.strip('\t') for x in fb.split('\r\n')]
        einf = float(Ls[0].split('\t')[1]); cur=None; oscs=[]
        for ln in Ls[1:]:
            s = ln.strip()
            if s.startswith("'") and s.endswith("'"):
                cur = (s.strip("'"), []); oscs.append(cur)
            elif "\t'" in ln and cur is not None:
                f = ln.split('\t'); cur[1].append((f[5].strip("'"), float(f[0])))
        layers.append((einf, oscs))
    return rough, dox, dth, thick, layers

def N_of(einf, oscs, wl):
    E = 1239.841984/wl; e1 = np.zeros_like(wl); e2 = np.zeros_like(wl)
    for typ, ps in oscs:
        d = dict(ps)
        if typ == 'Drude(RT)':
            a,b = osc.drude_rt(E, d['Resistivity (Ohm\u00b7cm)'], d['Scat. Time (fs)'])
        elif typ == 'Tauc-Lorentz':
            a,b = osc.tauc_lorentz(E, d['Amp'], d['Br'], d['Eo'], d['Eg'])
        elif typ == 'Gaussian':
            a,b = osc.gaussian(E, d['Amp'], d['Br'], d['En'])
        else:
            raise ValueError(typ)
        e1 += a; e2 += b
    N = np.sqrt((einf + e1 + 1j*e2).astype(complex))
    return np.where(N.imag<0,-N,N)

R = json.load(open('ag_final_result.json'))
print('%-24s %-6s %8s %8s   %s'%('file','sheet','MSE(mod)','MSE(fit)','match'))
for f in sorted(os.listdir(D)):
    if not f.endswith('.mod'): continue
    sh = f.split('_')[0]
    rough, dox, dth, thick, layers = parse(os.path.join(D,f))
    wl,P,Dl,_ = L.load(sh)
    m = (wl>=260)&(wl<=1080); wl,P,Dl = wl[m],P[m],Dl[m]
    ox, si = cf._mat('NTVE_JAW',wl), cf._mat('SI_JAW',wl)
    Ns = [N_of(e,o,wl) for e,o in layers]
    amb = np.ones(len(wl))
    stack = [amb] + ([ef.bruggeman_ema50(Ns[-1],amb)] if rough>0 else []) + Ns[::-1] + [ox, si]
    ds    = ([rough] if rough>0 else []) + thick[::-1] + [dox]
    M = [cf.ncs(P[:,i],Dl[:,i]) for i in range(5)]; o=[]
    for a,an in enumerate(L.ANG+dth):
        rp,rs = ef._tmm(wl, stack, ds, an)
        rr=rp/rs; q,dl=np.arctan(np.abs(rr)),np.angle(rr)
        o += [np.cos(2*q)-M[a][0], np.sin(2*q)*np.cos(dl)-M[a][1], np.sin(2*q)*np.sin(dl)-M[a][2]]
    res=np.concatenate(o); n=len(res)//3
    npar = 14 if len(layers)>1 else 6
    mse = 1000*np.sqrt(np.sum(res**2)/(3*n-npar))
    ref = R[sh]['mse'] if sh in R else None
    print('%-24s %-6s %8.2f %8s   %s'%(f[:24],sh,mse,('%.2f'%ref) if ref else '   -   ',
          'OK' if (ref is None or abs(mse-ref)<0.3) else 'MISMATCH'))

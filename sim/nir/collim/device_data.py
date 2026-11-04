"""Per-wavelength device data for a fixed NIR stack (nir_opt model): Purcell factor F, eta_sub, substrate angular power
P(theta) (1-deg bins) and stack reflectance R(theta) for light returning from the substrate.  600-900 nm, 2-nm steps.
Usage: python3 device_data.py <emitter-json from ../results_emitter> <Theta> <out.npz>"""
import sys, json, numpy as np
sys.path.insert(0, '/home/user/OLED-/sim/nir')
import nir_opt as N, cps2
from multiprocessing import Pool
def one(a):
    lam, x, te, nd, nsub, th = a; N.THETA = th
    S = cps2.Stack(float(lam), (N.NORG, N.NORG), N.DEML, N.DEML / 2, above=[(N.NORG, N.NORG, x[0]), (N.AG[int(lam - 600)],) * 2 + (100.0,)],
                   below=N.stack(x, te, nd, lam), n_sub=nsub, h=th)
    r = cps2.solve_pol(S, npts=3000)
    return float(np.real(r['P_tot'])), float(np.real(r['air'] + r['sub'])), np.real(cps2.sub_angular(S, N.TH)), np.real(cps2.stack_reflectance(S, N.TH))
if __name__ == '__main__':
    d = json.load(open(sys.argv[1])); th = float(sys.argv[2]); x = np.array(list(d['structure'].values()), float)
    lams = np.arange(600, 901, 2)
    with Pool(4) as p: res = p.map(one, [(l, x, d['TE'], d['ndbr'], d['n_sub'], th) for l in lams])
    np.savez(sys.argv[3], lam=lams, F=np.array([r[0] for r in res]), es=np.array([r[1] for r in res]),
             P=np.array([r[2] for r in res]), R=np.array([r[3] for r in res]), x=x, n_sub=d['n_sub'], TE=d['TE'], ndbr=d['ndbr'])
    print('saved', sys.argv[3], flush=True)

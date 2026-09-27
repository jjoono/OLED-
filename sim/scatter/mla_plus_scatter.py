"""Hemisphere MLA (n 1.8) on top of an index-matched scattering layer (n 1.8, HG) on the n 1.8 substrate.
The slab is index-matched on both sides (no Fresnel), so its up/down responses are symmetric: T[out,in], R[out,in].
Combined structure (adding method on 1-deg angular bins):
  u   = (I - R B_R^MLA)^-1 T v          light going up into the MLA, all internal bounces summed
  B_T^eff = B_T^MLA . (I - R B_R^MLA)^-1 T
  B_R^eff = R + T B_R^MLA (I - R B_R^MLA)^-1 T
then the device series as before."""
import numpy as np, scipy.io as sio, sys, csv
from multiprocessing import Pool
sys.path[:0] = ['/home/user/OLED-/sim/fig3c']
import cps2, materials as M
TH = np.arange(90) + 0.5; g0 = np.cos(np.radians(TH))*np.sin(np.radians(TH)); g0 /= g0.sum()
def hg_cos(g, u):
    if abs(g) < 1e-6: return 2*u - 1
    return (1 + g*g - ((1 - g*g)/(1 - g + 2*g*u))**2)/(2*g)
def slab(S, g, N=6000, seed=1, albedo=1.0):
    rng = np.random.default_rng(seed); T = np.zeros((90, 90)); R = np.zeros((90, 90))
    for a in range(90):
        th = np.radians(a + 0.5); ph = rng.uniform(0, 2*np.pi, N); st = np.sin(th)
        d = np.stack([st*np.cos(ph), st*np.sin(ph), np.full(N, np.cos(th))], 1); z = np.zeros(N); alive = np.ones(N, bool)
        while alive.any():
            idx = np.where(alive)[0]; znew = z[idx] - np.log(rng.random(idx.size))/S*d[idx, 2]
            top = znew >= 1; bot = znew <= 0; mid = ~(top | bot)
            for m, arr, sgn in ((top, T, 1), (bot, R, -1)):
                ii = idx[m]
                if ii.size:
                    ang = np.degrees(np.arccos(np.clip(sgn*d[ii, 2], 0, 1))).astype(int).clip(0, 89)
                    np.add.at(arr[:, a], ang, 1.0/N); alive[ii] = False
            ii = idx[mid]; z[ii] = znew[mid]
            if albedo < 1:
                dead = rng.random(ii.size) > albedo; alive[ii[dead]] = False; ii = ii[~dead]
            ct = hg_cos(g, rng.random(ii.size)); stt = np.sqrt(np.clip(1 - ct**2, 0, 1)); p2 = rng.uniform(0, 2*np.pi, ii.size)
            dd = d[ii]; ux, uy, uz = dd[:, 0], dd[:, 1], dd[:, 2]; den = np.sqrt(np.clip(1 - uz**2, 1e-12, None))
            nx = stt*(ux*uz*np.cos(p2) - uy*np.sin(p2))/den + ux*ct
            ny = stt*(uy*uz*np.cos(p2) + ux*np.sin(p2))/den + uy*ct
            nz = -stt*np.cos(p2)*den + uz*ct
            d[ii] = np.stack([nx, ny, nz], 1)
    return T, R
mla = sio.loadmat('/home/user/OLED-/sim/mla/lt_hemisphere_bsdf.mat')['BSDF_MLA'][:, :, 10]
BTm, BRm = mla[:90].sum(0), mla[90:][::-1]
EML = (1.83273, 1.67122); TAPC = (1.69075, 1.66375); B3 = M.ETL['B3PyMPM']
St = cps2.Stack(550.0, EML, 25.0, 12.5, above=[(B3[0], B3[1], 300.0), (M.AG, M.AG, 100.0)],
                below=[(TAPC[0], TAPC[1], 300.0), (M.ITO, M.ITO, 50.0)], n_sub=1.8)
Rd = np.real(cps2.stack_reflectance(St, TH)); P = np.real(cps2.sub_angular(St, TH)); P /= P.sum()
def series(BT, BR, K=400):
    v = P.copy(); e = 0; ps = []
    for k in range(K):
        U = v.sum(); t = BT @ v; e += t
        if k < 8: ps.append(t/U)
        v = Rd*(BR @ v)
    return e, ps
def combo(a):
    S, g, alb = a; T, R = slab(S, g, seed=int(S*97 + g*13), albedo=alb)
    X = np.linalg.solve(np.eye(90) - R @ BRm, T)
    BT = BTm @ X; BR = R + T @ BRm @ X
    e, ps = series(BT, BR); pl = BT @ g0
    return (S, g, alb, e, pl, ps)
if __name__ == '__main__':
    e0, ps0 = series(BTm, BRm)
    print('MLA alone      eta_ext %.4f  p_Lamb %.3f  p_1..p_6 %s' % (e0, BTm @ g0, ' '.join('%.3f' % x for x in ps0[:6])))
    jobs = [(S, g, 1.0) for S in (0.3, 1, 3, 10) for g in (0.0, 0.5, 0.9)] + [(S, g, a) for S, g in ((1, 0.5), (3, 0.9)) for a in (0.999, 0.995)]
    with Pool(4) as p: res = p.map(combo, jobs)
    rows = [['S', 'g', 'albedo', 'eta_ext', 'p_Lamb', 'p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6']]
    for S, g, a, e, pl, ps in res:
        rows.append([S, g, a, e, pl] + ps[:6])
        print('MLA+scat S%-4s g%-4s alb %-6s eta_ext %.4f  p_Lamb %.3f  p_1..p_6 %s' % (S, g, a, e, pl, ' '.join('%.3f' % x for x in ps[:6])))
    with open('mla_plus_scatter.csv', 'w', newline='') as f: csv.writer(f).writerows(rows)

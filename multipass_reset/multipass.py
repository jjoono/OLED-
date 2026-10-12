"""Multi-pass extraction with an angle-preserving rear (specular reflector) versus a
rear that resets the returned light to a Lambertian distribution on every pass.

BSDFs: hftrace.py (validated against LightTools, p = 0.4275 vs 0.4274).
Device: CPS dipole model of the preprint (cps2), 550 nm, Al or Ag cathode,
ITO 150 nm, n_sub = n_lens = 1.51, Theta = 0.865, cavity (dETL, dHTL) of each design.
"""
import sys, numpy as np, csv, json
SP = '/tmp/claude-0/-home-user-OLED-/a92e5b80-064b-5ded-9e05-afa68bed9893/scratchpad'
sys.path.insert(0, SP + '/pre/sim/fig3c')
import cps2, materials as M

AL = 0.958336812040867 + 6.68678039868499j
TH = np.arange(90) + 0.5
W = np.cos(np.radians(TH)) * np.sin(np.radians(TH)); Wn = W / W.sum()
BANDS = [(0, 20), (20, 40), (40, 60), (60, 80)]
LAMB = np.array([np.sin(np.radians(b)) ** 2 - np.sin(np.radians(a)) ** 2 for a, b in BANDS])
D = np.load(SP + '/mp/bsdfs.npz')
NAMES = ['0-20', '20-40', '40-60', '60-80', 'total']


def device(refl, detl, dhtl):
    S = cps2.Stack(550.0, (1.8, 1.8), 20.0, 10.0, above=[(1.8, 1.8, detl), (refl, refl, 100.0)],
                   below=[(1.8, 1.8, dhtl), (M.ITO, M.ITO, 150.0)], n_sub=1.51, h=0.865)
    r = cps2.solve_pol(S, npts=6000)
    R = np.real(cps2.stack_reflectance(S, TH)); P = np.real(cps2.sub_angular(S, TH))
    return r['air'] + r['sub'], R, P / P.sum()


def run(BT, BR, R, P, mode='spec', m=1.0, dA=0.0, K=400):
    """mode 'spec': v' = R * w ; 'reset': randomized fraction m returns as Lambertian
    with round-trip loss A'_lamb + dA. Returns eta_ext, air distribution, p_k list."""
    BTm, BRm = BT.T, BR.T                       # [out, in]
    Alamb = 1 - (R * W).sum() / W.sum()
    v = P.copy(); air = np.zeros(90); pk = []
    for k in range(K):
        U = v.sum()
        if U < 1e-9:
            break
        a = BTm @ v; air += a; pk.append(a.sum() / U)
        w = BRm @ v
        if mode == 'spec':
            v = R * w
        else:
            v = (1 - m) * (R * w) + m * (1 - Alamb - dA) * w.sum() * Wn
    return air.sum(), air, np.array(pk)


def comp(air):
    return np.array([air[a:b].sum() for a, b in BANDS]) / air.sum()


if __name__ == '__main__':
    designs = sorted({k.rsplit('_', 1)[0] for k in D.files if k.endswith('_BT')} - {'free_0-20'})  # vertical rim wall: not traceable as a height field
    rows = []; res = {}
    for d in designs:
        BT, BR, xv = D[d + '_BT'], D[d + '_BR'], D[d + '_x']
        for cat, refl in (('Al', AL), ('Ag', M.AG)):
            es, R, P = device(refl, xv[-3], xv[-2])
            Al_ = 1 - (R * W).sum() / W.sum()
            p_l = (BT.sum(1) * Wn).sum()
            for mode, dA in (('spec', 0.0), ('reset', 0.0), ('reset', 0.01), ('reset', 0.02), ('reset', 0.03)):
                e, air, pk = run(BT, BR, R, P, mode, 1.0, dA)
                S = comp(air); bandEQE = es * e * S
                key = f'{d}|{cat}|{mode}|{dA}'
                res[key] = dict(eta_sub=es, eta_ext=e, EQE=es * e, S=S.tolist(), band=bandEQE.tolist(),
                                pk=pk[:8].tolist(), Aprime=Al_, p_lamb=p_l)
                rows.append([d, cat, mode, dA, round(es, 4), round(Al_, 4), round(p_l, 4), round(e, 4), round(es * e, 4),
                             *np.round(S, 4), *np.round(bandEQE, 5), *np.round(pk[:6], 4)])
            print(d, cat, 'spec EQE %.4f  reset EQE %.4f' % (res[f'{d}|{cat}|spec|0.0']['EQE'], res[f'{d}|{cat}|reset|0.0']['EQE']), flush=True)
    # mixing-fraction sweep on the total-EQE designs
    msweep = []
    for d in ('hemi_total', 'free_total'):
        BT, BR, xv = D[d + '_BT'], D[d + '_BR'], D[d + '_x']
        for cat, refl in (('Al', AL), ('Ag', M.AG)):
            es, R, P = device(refl, xv[-3], xv[-2])
            for m in np.linspace(0, 1, 11):
                e, air, pk = run(BT, BR, R, P, 'reset', m, 0.0)
                msweep.append([d, cat, round(m, 2), round(e, 4), round(es * e, 4), *np.round(comp(air), 4)])
    hdr = ['design', 'cathode', 'rear', 'dA_extra', 'eta_sub', "A'(Lamb)", 'p(Lamb)', 'eta_ext', 'EQE',
           'S_0-20', 'S_20-40', 'S_40-60', 'S_60-80', 'EQE_0-20', 'EQE_20-40', 'EQE_40-60', 'EQE_60-80',
           'p_1', 'p_2', 'p_3', 'p_4', 'p_5', 'p_6']
    with open(SP + '/mp/multipass_results.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(hdr); w.writerows(rows)
    with open(SP + '/mp/mixing_sweep.csv', 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['design', 'cathode', 'm', 'eta_ext', 'EQE', 'S_0-20', 'S_20-40', 'S_40-60', 'S_60-80']); w.writerows(msweep)
    json.dump(res, open(SP + '/mp/multipass_results.json', 'w'))
    print('Lambertian partition', np.round(LAMB, 3))

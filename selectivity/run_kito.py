"""Round-trip loss set by the transparent electrode (design rule 2 of the round-trip-loss paper):
low-loss n_sub = 1.8 stack with the ITO extinction spectrum scaled to k(550) = 0 ... 0.01.
For each, no filter / ideal 10-20-30 deg cones / best TiO2-SiO2 16-pair edge DBR (full randomizer)."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, json
from multiprocessing import Pool
import sel, common

K550 = sel.D['n_ITO'][150].imag
SW = common.Sweep()


def one(k):
    base = sel.D['n_ITO'].copy()
    sel.D['n_ITO'] = base.real + 1j * base.imag * (k / K550)
    dev = sel.device('LL18'); sel.D['n_ITO'] = base
    A = 1 - (dev['R'] * sel.WL).sum(1)
    Ts, Rs = sel.mixer(1.0); out = {'k': k}
    for sp in ('G30', 'v-DABNA (HF)'):
        spec = common.SPECTRA[sp]; o = {'A_prime': float(spec @ A), 'eta_sub': float(spec @ dev['eta_sub'])}
        Rf, Tf = sel.filter_RT('fresnel', 1.8)
        o['nofilter'] = sel.metrics(dev['eta_sub'], sel.cones(sel.solve(dev, Rf, Tf, Ts, Rs), 1.8), spec)
        for b in (10.0, 20.0, 30.0):
            Rf, Tf = sel.filter_RT(('ideal', b), 1.8)
            o[f'ideal{b:g}'] = sel.metrics(dev['eta_sub'], sel.cones(sel.solve(dev, Rf, Tf, Ts, Rs), 1.8), spec)
            d, _ = SW.best_edge('LL18', 'full', spec, b, mat='TiO2/SiO2', N=16)
            Rf, Tf = sel.filter_RT(d, 1.8)
            o[f'dbr{b:g}'] = sel.metrics(dev['eta_sub'], sel.cones(sel.solve(dev, Rf, Tf, Ts, Rs), 1.8), spec)
            o[f'dbr{b:g}_design'] = str(d)
        out[sp] = o
    return out


if __name__ == '__main__':
    with Pool(4) as p:
        res = p.map(one, [0.0, 0.001, 0.002, 0.0032, 0.005, 0.0075, 0.01])
    json.dump(res, open('data/kito.json', 'w'))
    for r in res:
        o = r['G30']; print(r['k'], round(o['A_prime'], 4), round(o['nofilter']['EQE_total'], 3), round(o['ideal20']['EQE_20'], 3), round(o['dbr20']['EQE_20'], 3))

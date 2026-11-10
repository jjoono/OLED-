"""Thin-Ag (20 nm) microcavity bottom emitter + randomizing scattering layer on the glass exit (no filter / edge DBR).
Cavity scanned (ETL/HTL 30 nm steps); full-randomizer closed form for the scan, then matrix model with the
physical slab (g 0.8, tau* 1) for the best cavities.  Reports total EQE, on-axis cd/A, du'v'."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
import json, numpy as np
from multiprocessing import Pool
import sel, common, run_display as RD, display_refs as DR

E = np.arange(30.0, 301.0, 30.0)


def one(args):
    cath, e, h = args
    d = DR.thin_ag_cavity(cath, e, h)
    return args, d


if __name__ == '__main__':
    jobs = [(c, e, h) for c in ('Ag', 'Al') for e in E for h in E]
    with Pool(3) as p:
        devs = dict(p.map(one, jobs, chunksize=4))
    Rf, Tf = sel.filter_RT('fresnel', 1.51); F = sel.cone_frac(1.51)
    pt = (Tf * sel.WL).sum(1); p5 = (Tf * sel.WL * F[0]).sum(1); r = (Rf * sel.WL).sum(1)
    slab = tuple(RD.SL['g0.8_t1.0_w1']); out = {}
    for sp in ('G30', 'Ir(ppy)2acac'):
        spec = common.SPECTRA[sp]
        for cath in ('Ag', 'Al'):
            sc = {}
            for (c, e, h), d in devs.items():
                if c != cath: continue
                rho = (d['R'] * sel.WL).sum(1); den = 1 - r * rho
                sc[(e, h)] = (spec @ (d['eta_sub'] * pt / den), spec @ (d['eta_sub'] * p5 / den), spec @ (1 - rho))
            for obj, k in (('max EQE', 0), ('max on-axis', 1)):
                e, h = max(sc, key=lambda q: sc[q][k]); d = devs[(cath, e, h)]
                for nm, filt in (('scatterer', 'fresnel'), ('planar (no scatterer)', None)):
                    if filt is None:
                        I = RD.to_air_angle(RD.planar_air(d), 1.51)
                    else:
                        I = RD.to_air_angle(sel.solve(d, *sel.filter_RT(filt, 1.51), *slab), 1.51)
                    m = RD.display_metrics(I, d['eta_sub'], spec)
                    key = f'{sp}|{cath} cathode, 20nm Ag anode, ETL {e:.0f} HTL {h:.0f} ({obj})|{nm}'
                    out[key] = dict(cdA0=m['cdA0'], EQE=m['EQE'], duv30=m['duv30'], duv60=m['duv60'], half=m['half_lum_angle'], Aprime=float(sc[(e, h)][2]))
                    print(f"{key:85s} cd/A {m['cdA0']:5.0f} EQE {m['EQE']:.3f} du'v'30 {m['duv30']:.3f} A' {sc[(e, h)][2]:.3f}", flush=True)
    json.dump(out, open('data/cavscat.json', 'w'), indent=1)

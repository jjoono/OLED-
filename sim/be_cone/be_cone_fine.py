"""Fine-spectrum (2 nm) confirmation of the be_cone_opt.py optima.  The 10 nm spectral grid overrates some thick
resonant structures, so for every (IZO, extractor, bound) the 8 best mutually separated (>= 30 nm) grid points are
re-evaluated on the 2 nm grid, the best is refined locally (+-10 nm, step 5), and results are made monotonic in the
bound (an optimum under a smaller bound is admissible under a larger one).  -> be_cone_best.csv"""
import csv, itertools, numpy as np
from multiprocessing import Pool
import be_cone_opt as B

G = {}
for r in csv.DictReader(open('be_cone_grid.csv')):
    G[(int(r['d_ETL']), int(r['d_HTL']), r['IZO'])] = {k: float(v) for k, v in r.items() if k.startswith('EQE_')}
BOUNDS = (100, 200, 300, 500); EXTS = list(B.EXT)


def top(izo, ext, bound, n=8, sep=30):
    c = sorted(((v[f'EQE_{ext}'], k) for k, v in G.items() if k[2] == izo and k[0] <= bound and k[1] <= bound), reverse=True)
    out = []
    for _, k in c:
        if all(abs(k[0] - o[0]) >= sep or abs(k[1] - o[1]) >= sep for o in out): out.append(k)
        if len(out) == n: break
    return out


if __name__ == '__main__':
    combos = list(itertools.product(B.IZOS, EXTS, BOUNDS))
    cand = {k for cb in combos for k in top(*cb)}
    with Pool(4) as pool:
        F = B.evaluate(sorted(cand), 2, pool)
        loc = set()
        for izo, ext, bd in combos:
            (e, h, _), _ = max(((k, F[k]) for k in top(izo, ext, bd)), key=lambda kv: kv[1][f'EQE_{ext}'])
            loc |= {(x, y, izo) for x in range(max(20, e - 10), min(bd, e + 10) + 1, 5) for y in range(max(20, h - 10), min(bd, h + 10) + 1, 5)}
        F.update(B.evaluate(sorted(loc - set(F)), 2, pool))
    rows = []
    for izo, ext in itertools.product(B.IZOS, EXTS):
        for bd in BOUNDS:
            k, v = max(((k, v) for k, v in F.items() if k[2] == izo and k[0] <= bd and k[1] <= bd), key=lambda kv: kv[1][f'EQE_{ext}'])
            rows.append(dict(IZO=izo, extractor=ext, bound=bd, d_ETL=k[0], d_HTL=k[1], EQE=v[f'EQE_{ext}'], eta_sub=v['eta_sub'],
                             eta_ext=v[f'EQE_{ext}'] / v['eta_sub'], spp=v['spp'], EQE_grid10=G.get(k, {}).get(f'EQE_{ext}', float('nan'))))
            print({a: (round(b, 4) if isinstance(b, float) else b) for a, b in rows[-1].items()}, flush=True)
    with open('be_cone_best.csv', 'w', newline='') as f:
        wr = csv.DictWriter(f, list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)

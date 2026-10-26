"""Tabulate results(_refined)/*.json -> summary.csv / summary.md.  EQE(q) = q F EQE1 / (1 - q + q F)."""
import json, glob, os, csv
rows = []
for f in sorted(glob.glob('results/*.json')):
    g = f.replace('results/', 'results_refined/'); d = json.load(open(g if os.path.exists(g) else f))
    F, e1 = d['F'], d['EQE_q1.0']; eq = lambda q: q * F * e1 / (1 - q + q * F)
    s = d['structure']; dbr = '/'.join('%.0f' % v for k, v in s.items() if k[:4] in ('SiO2', 'TiO2'))
    rows.append(dict(TE=d['TE'], n_sub=d['n_sub'], ext=d['extractor'], ndbr=d['ndbr'], ETL=s['d_ETL'], HTL=s['d_HTL'], d_TE=s['d_TE'],
                     DBR=dbr or '-', F=round(F, 3), eta_sub=round(d['eta_sub'], 3), eta_ext=round(d['eta_ext'], 3), A_round=round(d['A_round'], 3),
                     spp=round(d['spp'], 3), absorb=round(d['absorb'], 3),
                     **{'EQE_q%.1f' % q: round(100 * eq(q), 1) for q in (0.1, 0.3, 0.5, 0.8, 1.0)}))
with open('summary.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, rows[0].keys()); w.writeheader(); w.writerows(rows)
with open('summary.md', 'w') as fh:
    fh.write('| ' + ' | '.join(rows[0]) + ' |\n|' + '---|' * len(rows[0]) + '\n')
    for r in rows: fh.write('| ' + ' | '.join(str(v) for v in r.values()) + ' |\n')
print(open('summary.md').read())

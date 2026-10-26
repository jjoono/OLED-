"""Tabulate results_emitter/*.json -> summary_emitters.csv/.md"""
import json, glob, csv
rows = []
for f in sorted(glob.glob('results_emitter/*.json')):
    d = json.load(open(f)); s = d['structure']
    rows.append(dict(emitter=d['emitter'], peak=d['lam0'], PLQY=d['PLQY'], Theta=d['Theta'], TE=d['TE'], n_sub=d['n_sub'], ext=d['extractor'],
                     ndbr=d['ndbr'], ETL=s['d_ETL'], HTL=s['d_HTL'], d_TE=s['d_TE'],
                     DBR='/'.join('%.0f' % v for k, v in s.items() if k[:4] in ('SiO2', 'TiO2')) or '-',
                     F=round(d['F'], 3), eta_sub=round(d['eta_sub'], 3), eta_ext=round(d['eta_ext'], 3), A_round=round(d['A_round'], 3),
                     EQE=round(100 * d['EQE'], 1), EQE_q1=round(100 * d['EQE_q1'], 1)))
with open('summary_emitters.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, rows[0].keys()); w.writeheader(); w.writerows(rows)
with open('summary_emitters.md', 'w') as fh:
    fh.write('| ' + ' | '.join(rows[0]) + ' |\n|' + '---|' * len(rows[0]) + '\n')
    for r in rows: fh.write('| ' + ' | '.join(str(v) for v in r.values()) + ' |\n')
print(open('summary_emitters.md').read())

"""Score green_joint_fit raw results against tilt-corrected targets: common eta (LSQ on front EQE), per-device
profile RMS and z-scores. Usage: python3 green_joint_score.py raw.json targets.json out.json"""
import sys, json, numpy as np
R = json.load(open(sys.argv[1])); T = json.load(open(sys.argv[2])); DEV = ('bare', 'hemisphere', 'freeform')
SP, SE = 0.03, 1.0  # profile tolerance (I/I0), front-EQE tolerance (pp)
out = []
for r in R:
    fm = np.array([np.mean(T[d]['front']) for d in DEV]); fs = np.array([r[d]['front'] * 100 for d in DEV])
    eta = float(fs @ fm / (fs @ fs)); row = dict(etl=r['etl'], htl=r['htl'], eta=eta, cost=0.0)
    for k, d in enumerate(DEV):
        m = np.array([np.nan if v is None else v for v in T[d]['mean']]); s = np.array(r[d]['rel'])
        ok = ~np.isnan(m); ok[0] = False
        if d != 'bare': ok[7:] = False
        res = s[ok] - m[ok]; prms = float(np.sqrt(np.mean(res ** 2))); pmax = float(np.abs(res).max()); de = eta * fs[k] - fm[k]
        row[d] = dict(front_sim=eta * fs[k], front_meas=fm[k], dEQE=de, prof_rms=prms, prof_max=pmax, PE=r[d]['PE'])
        row['cost'] += (prms / SP) ** 2 + (de / SE) ** 2
    row['within'] = all(row[d]['prof_rms'] < SP and abs(row[d]['dEQE']) < SE for d in DEV)
    out.append(row)
out.sort(key=lambda x: x['cost']); json.dump(out, open(sys.argv[3], 'w'), indent=1)
print('ETL HTL  eta  cost | bare dE rms | hemi dE rms | ff dE rms')
for o in out[:15]:
    print(f"{o['etl']:3.0f} {o['htl']:3.0f} {o['eta']:.3f} {o['cost']:6.1f} | " + ' | '.join(f"{o[d]['dEQE']:+5.1f} {o[d]['prof_rms']:.3f}" for d in DEV) + (' *' if o['within'] else ''))

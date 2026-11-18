"""Prospective test of the reduced-model applicability conditions (main text Sec. 4.3) on twelve new structures traced
by multipass_reset/prospective.py. Step 1 writes data/prospective_predictions.json from M and the first-pass
distribution only (Al, 550 nm): 'valid' if lambda_2 <= 0.999 and lambda_2 - lambda_eff <= 0.35, otherwise 'outside'.
Step 2 evaluates the reduced model against the full series (three devices x three wavelengths) and writes
data/prospective_results.json."""
import json, glob, os, time, numpy as np, mixing
src = open('oe_figs.py').read(); exec(src[:src.index('# ------------------------------------------------------------------ Fig. 1')])
TH = np.arange(90) + .5
def modes(BR, dev='Al', i=75, cut=85):
    k = TH < cut; B = BR[np.ix_(k, k)]; cs = B.sum(0); M = B / np.maximum(cs, 1e-12)[None, :]
    w = W1[k]; pi = w * cs; pi /= pi.sum()
    S = np.diag(pi ** -0.5) @ M @ np.diag(pi ** 0.5); S = 0.5 * (S + S.T)
    lam, U = np.linalg.eigh(S); o = np.argsort(-lam); lam, U = lam[o], U[:, o]
    P = DV[f'{dev}_P'][i].reshape(90, 2).sum(1)[k]; u0 = B @ P; u0 /= u0.sum()
    a = U.T @ (u0 / np.sqrt(pi)); a2 = a[1:] ** 2; a2 /= a2.sum()
    return dict(lam_eff=float((a2 * lam[1:]).sum()), w_slow=float(a2[lam[1:] > 0.8].sum()))
S = {}
for f in sorted(glob.glob('../multipass_reset/prospective/*.npz')):
    z = np.load(f); S[os.path.basename(f)[:-4]] = (z['BT'], z['BR'].T)
pred = {'criterion': 'valid if lambda_2 <= 0.999 and lambda_2 - lambda_eff <= 0.35 (Al, 550 nm first pass)', 'written': time.strftime('%Y-%m-%d %H:%M:%S'), 'structures': {}}
for nm, (BT, BR) in S.items():
    l2 = mixing.lam2_cut(BR); m = modes(BR); gap = l2 - m['lam_eff']
    pred['structures'][nm] = dict(lam2=l2, lam_eff=m['lam_eff'], gap=gap, w_slow=m['w_slow'], prediction='valid' if (l2 <= 0.999 and gap <= 0.35) else 'outside')
json.dump(pred, open('data/prospective_predictions.json', 'w'), indent=1)
print('predictions written'); [print(f"  {k:16s} lam2 {v['lam2']:.3f} lam_eff {v['lam_eff']:.3f} gap {v['gap']:.3f} -> {v['prediction']}") for k, v in pred['structures'].items()]
res = {}
for nm, (BT, BR) in S.items():
    pT = BT.sum(1); l2 = pred['structures'][nm]['lam2']; errs = []; rand = []
    for dev in DEVS:
        for i, lam in WLS:
            P = DV[f'{dev}_P'][i].reshape(90, 2).sum(1); R = DV[f'{dev}_R'][i].reshape(90, 2).mean(1)
            full = mixing.full(pT, BR, R, P); red = mixing.reduced(pT, W1, R, P, max(1 - l2, 0.0))
            errs.append((red - full) * 100); rand.append((mixing.closed(pT, W1, R) - full) * 100)
    e = np.abs(errs); res[nm] = dict(**pred['structures'][nm], p=float(pT @ (W1 / W1.sum())), err_mean=float(e.mean()), err_max=float(e.max()), rand_mean=float(np.mean(np.abs(rand))), rand_max=float(np.max(np.abs(rand))))
json.dump(res, open('data/prospective_results.json', 'w'), indent=1)
print('results')
for k, v in res.items(): print(f"  {k:16s} {v['prediction']:8s} p {v['p']:.3f} reduced err mean {v['err_mean']:.2f} max {v['err_max']:.2f} | fully randomized mean {v['rand_mean']:.2f} max {v['rand_max']:.2f}")
val = [v for v in res.values() if v['prediction'] == 'valid']; out = [v for v in res.values() if v['prediction'] != 'valid']
print('valid: n', len(val), 'mean', np.mean([v['err_mean'] for v in val]).round(2), 'max', max(v['err_max'] for v in val) if val else None)
print('outside: n', len(out), 'max errors', [round(v['err_max'], 2) for v in out])

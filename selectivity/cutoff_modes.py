"""(1) Reduced-model error and metric values versus the angular cutoff (80, 85, 88, 90 deg), 342 library cases.
(2) Why a single rate fails for cone arrays: overlap of the light that actually circulates with the slow modes of M,
and the observable-weighted relaxation lam_eff (Rayleigh quotient of M for the deviation of the returned light from the
stationary distribution, etendue inner product; defined in mixing.mixing_defs before the blind test).
Output: data/cutoff_modes.json"""
import json, glob, os, numpy as np, mixing
src = open('oe_figs.py').read(); exec(src[:src.index('# ------------------------------------------------------------------ Fig. 1')])
TH = np.arange(90) + .5
def lam2_at(BR, cut):
    if cut >= 90:
        B = BR
    else:
        k = TH < cut; B = BR[np.ix_(k, k)]
    M = B / np.maximum(B.sum(0, keepdims=True), 1e-12); ev = np.linalg.eigvals(M); ev = ev[np.argsort(-np.abs(ev))]
    return float(np.abs(ev[1]))
out = {'cutoff': {}, 'modes': {}}
for cut in (80, 85, 88, 90):
    errs = []; byc = {}
    for nm, (c, BT, BR) in LIB.items():
        l2 = lam2_at(BR, cut); pT = BT.sum(1)
        for dev in DEVS:
            for i, lam in WLS:
                P = DV[f'{dev}_P'][i].reshape(90, 2).sum(1); R = DV[f'{dev}_R'][i].reshape(90, 2).mean(1)
                e = abs(mixing.reduced(pT, W1, R, P, max(1 - l2, 0)) - EV[(nm, dev, lam)]['eta']) * 100
                errs.append(e); byc.setdefault(c, []).append(e)
    out['cutoff'][cut] = dict(mean=float(np.mean(errs)), max=float(np.max(errs)), by_class={k: [float(np.mean(v)), float(np.max(v))] for k, v in byc.items()})
    print(f'cutoff {cut}: reduced-model error mean {np.mean(errs):.2f} max {np.max(errs):.2f} pp |', {k: round(float(np.mean(v)), 2) for k, v in byc.items()})

def modes(BR, BT, dev='Al', i=75, cut=85):
    """Overlap of the first returned light with the eigenmodes of M (etendue-weighted inner product)."""
    k = TH < cut; B = BR[np.ix_(k, k)]; cs = B.sum(0); M = B / np.maximum(cs, 1e-12)[None, :]
    w = W1[k]; pi = w * cs; pi /= pi.sum()                                     # stationary distribution of the truncated, renormalized M: w_i s_i
    S = np.diag(pi ** -0.5) @ M @ np.diag(pi ** 0.5); S = 0.5 * (S + S.T)           # symmetrized (reciprocal)
    lam, U = np.linalg.eigh(S); o = np.argsort(-lam); lam, U = lam[o], U[:, o]
    P = DV[f'{dev}_P'][i].reshape(90, 2).sum(1)[k]; u0 = B @ P; u0 /= u0.sum()    # light returned after the first pass
    a = U.T @ (u0 / np.sqrt(pi)); a2 = a[1:] ** 2; a2 /= a2.sum()                    # deviation from stationarity, by mode
    lam_eff = float((a2 * lam[1:]).sum())
    return dict(lam2=float(lam[1]), lam_eff=lam_eff, w_slow=float(a2[lam[1:] > 0.8].sum()), top=[[float(lam[j + 1]), float(a2[j])] for j in range(4)])

lib = {nm: (BT, BR) for nm, (c, BT, BR) in LIB.items()}
for f in sorted(glob.glob('../multipass_reset/blind/*.npz')):
    z = np.load(f); lib['blind:' + os.path.basename(f)[:-4]] = (z['BT'], z['BR'].T)
R_ = json.load(open('../multipass_reset/blind/results.json'))
for nm, (BT, BR) in lib.items():
    m = modes(BR, BT); pT = BT.sum(1)
    errs = {'lam2': [], 'lam_eff': []}
    for dev in DEVS:
        for i, lam in WLS:
            P = DV[f'{dev}_P'][i].reshape(90, 2).sum(1); R = DV[f'{dev}_R'][i].reshape(90, 2).mean(1)
            full = mixing.full(pT, BR, R, P)
            errs['lam2'].append(abs(mixing.reduced(pT, W1, R, P, max(1 - mixing.lam2_cut(BR), 0)) - full) * 100)
            errs['lam_eff'].append(abs(mixing.reduced(pT, W1, R, P, float(np.clip(1 - m['lam_eff'], 0, 1))) - full) * 100)
    m.update(err_lam2=[float(np.mean(errs['lam2'])), float(np.max(errs['lam2']))], err_lam_eff=[float(np.mean(errs['lam_eff'])), float(np.max(errs['lam_eff']))])
    out['modes'][nm] = m
json.dump(out, open('data/cutoff_modes.json', 'w'), indent=1)
for nm, m in out['modes'].items():
    if nm.startswith('blind') or nm in ('MLA hemi_total', 'MLA free_40-60', 'W_slope0.2', 'A_shape1.0', 'scatter g0.8_t1.0'):
        print(f"{nm:24s} lam2 {m['lam2']:.3f} lam_eff {m['lam_eff']:.3f} slow-mode weight {m['w_slow']:.3f} | err lam2 {m['err_lam2'][0]:.2f}/{m['err_lam2'][1]:.2f} lam_eff {m['err_lam_eff'][0]:.2f}/{m['err_lam_eff'][1]:.2f}")
lb = [m for nm, m in out['modes'].items() if not nm.startswith('blind') and nm != 'planar']
print('library (non-planar) lam_eff: mean err %.2f max %.2f | lam2 mean %.2f max %.2f' % (np.mean([m['err_lam_eff'][0] for m in lb]), max(m['err_lam_eff'][1] for m in lb), np.mean([m['err_lam2'][0] for m in lb]), max(m['err_lam2'][1] for m in lb)))

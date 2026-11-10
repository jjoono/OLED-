"""Ray-count and domain-size checks (review items 9, 16): compare metrics with the library traces."""
import json, glob, os, numpy as np
src = open('oe_figs.py').read(); exec(src[:src.index('LIB = library()')])
LIB = library()
def met(BT, BR):
    pT = BT.sum(1); l2 = mixing.lam2_cut(BR); pi, pt = mixing.pinf_cut(BR)
    P, R = DV['Al_P'][75].reshape(90, 2).sum(1), DV['Al_R'][75].reshape(90, 2).mean(1)
    e = mixing.full(pT, BR, R, P); er = mixing.reduced(pT, W1, R, P, 1 - l2)
    return dict(p=float(pT @ (W1 / W1.sum())), Nmix=1 / (1 - l2), xi=pi / pt, eta=e, err_pp=(er - e) * 100)
out = {}
for f in sorted(glob.glob('../multipass_reset/review_traces/*.npz')):
    nm = os.path.basename(f)[:-4]; z = np.load(f); out[nm] = met(z['BT'], z['BR'].T)
for ref in ('A_shape0.0', 'A_shape1.0', 'C_jitter0.3', 'W_slope0.2', 'W_slope1.0'):
    out[ref + '_library'] = met(LIB[ref][1], LIB[ref][2])
json.dump(out, open('data/review_traces.json', 'w'), indent=1)
for k in sorted(out): print(f"{k:22s} p {out[k]['p']:.4f} N {out[k]['Nmix']:.2f} xi {out[k]['xi']:.3f} eta {out[k]['eta']:.4f} err {out[k]['err_pp']:+.2f}")

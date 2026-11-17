import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_osc as O
W = _os.path.join(ELLIPS_OUT, '')
from multiprocessing import Pool
O.NG['NPB_22'] = 7; O.NG['NPB_21'] = 7
with Pool(2) as p:
    r = p.starmap(O.fit_material, [('NPB_22', ['1_22']), ('NPB_21', ['1_21'])])
out = {}
wp = np.array([350., 450., 550., 633.])
for name, (q, ds, per) in zip(('NPB_22', 'NPB_21'), r):
    N = O.N_of(q, wp)
    out[name] = dict(q=[float(v) for v in q], d=[float(v) for v in ds], mse=[float(v) for v in per])
    print('%-7s d=%.2f nm  MSE %.2f   %s' % (name, ds[0], per[0], '  '.join('%3.0f: %.3f+%.3fj' % (w, z.real, z.imag) for w, z in zip(wp, N))))
json.dump(out, open(W + 'osc_npb_sep.json', 'w'), indent=1)

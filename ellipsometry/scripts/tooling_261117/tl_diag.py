import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_osc as O
import jnk_ref as R, jnk_model as M
W = _os.path.join(ELLIPS_OUT, '')
OSC = json.load(open(W + 'osc.json'))
for mat, v in OSC.items():
    for sh, d in zip(v['sheets'], v['d']):
        wl, P, D, ang = L.se(sh, 245, 1000, 2); meas = L.ncs(P, D)
        Nf = O.N_of(np.array(v['q']), wl)
        n_, c_, s_ = M.se_ncs(wl, [R.n_air(wl), Nf, R.n_sio2(wl), R.n_si(wl)], [d, 2.0], ang)
        r = np.sqrt(((n_ - meas[0])**2 + (c_ - meas[1])**2 + (s_ - meas[2])**2).mean(1) / 3)
        bands = [(245, 300), (300, 350), (350, 400), (400, 450), (450, 600), (600, 1000)]
        print('%-6s %-5s ' % (mat, sh) + '  '.join('%d-%d:%5.1f' % (a, b, 1e3 * np.sqrt(np.mean(r[(wl >= a) & (wl < b)]**2))) for a, b in bands))

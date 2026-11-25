"""Si-piece thickness in the lab tooling template's convention (NTVE 0 A) next
to the measured-oxide one (33.68 A), so tooling factors can be compared like
with like against earlier runs fitted with that template."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_ce_genosc as C, tl_ce_thickness as T
W = _os.path.join(ELLIPS_OUT, '')

def job(sh, mat):
    return C.fit_thickness(sh, np.array(T.disp(mat)['p']), 0.0, mat in T.ROUGH)

if __name__ == '__main__':
    from multiprocessing import Pool
    TH = json.load(open(W + 'ce_thickness.json'))
    with Pool(4) as p:
        r = p.starmap(job, list(T.PIECES.items()))
    out = {}
    print('%-5s %-6s | NTVE 33.68 A: d  rough  MSE   | NTVE 0 A (template): d  rough  MSE' % ('sheet', 'mat'))
    for (sh, m), (d, mse, rg) in zip(T.PIECES.items(), r):
        t = TH[sh]; out[sh] = dict(mat=m, d_ntve0=d, rough_ntve0=rg, mse_ntve0=mse)
        print('%-5s %-6s |      %6.2f  %4.2f  %5.2f   |               %6.2f  %4.2f  %5.2f'
              % (sh, m, t['d'], t['rough'], t['mse'], d, rg, mse))
    json.dump(out, open(W + 'ce_thickness_ntve0.json', 'w'), indent=1)

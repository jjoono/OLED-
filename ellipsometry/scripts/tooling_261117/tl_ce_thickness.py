"""Every continuous-film Si piece against its material's CE-frame Gen-Osc
dispersion: thickness only, at the lab oxide and at a bare-wafer oxide, so the
one assumption that moves thickness is visible."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_ce_genosc as C

W = _os.path.join(ELLIPS_OUT, '')
PIECES = {'1_21': 'NPB', '1_22': 'NPB', '1_23': 'NPB', '1_24': 'TAPC', '2_21': 'Bphen', '2_22': 'Bphen'}
ROUGH = {'Bphen'}           # the only material whose fit improves with a roughness layer (11 %)
D_OX_BARE = 1.5            # nm, typical untreated native oxide

def disp(mat):
    CE = json.load(open(W + 'ce_genosc.json'))
    return CE[mat]['with_roughness' if mat in ROUGH else 'pick']

def job(sh, mat, d_ox):
    return C.fit_thickness(sh, np.array(disp(mat)['p']), d_ox, mat in ROUGH)

if __name__ == '__main__':
    from multiprocessing import Pool
    jobs = [(sh, m, dx) for sh, m in PIECES.items() for dx in (C.D_OX, D_OX_BARE)]
    with Pool(4) as p:
        r = p.starmap(job, jobs)
    out = {}
    print('%-5s %-6s | at %.2f nm oxide (lab, O2 plasma)            | at %.1f nm oxide' % ('sheet', 'mat', C.D_OX, D_OX_BARE))
    for i, (sh, m) in enumerate(PIECES.items()):
        (d1, m1, r1), (d2, m2, r2) = r[2 * i], r[2 * i + 1]
        out[sh] = dict(mat=m, d=d1, rough=r1, mse=m1, d_bare_oxide=d2, rough_bare_oxide=r2, mse_bare_oxide=m2,
                       dispersion='with_roughness' if m in ROUGH else 'pick')
        print('%-5s %-6s | d %6.2f  rough %4.2f nm  MSE %5.2f        | d %6.2f  rough %4.2f  MSE %5.2f'
              % (sh, m, d1, r1, m1, d2, r2, m2))
    json.dump(out, open(W + 'ce_thickness.json', 'w'), indent=1)

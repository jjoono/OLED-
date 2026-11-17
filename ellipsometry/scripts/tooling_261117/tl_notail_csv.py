import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # tooling.xlsx
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fits, tables, figures
import sys, json, csv, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_osc as O
W = _os.path.join(ELLIPS_OUT, '')
OUT = _os.path.join(ELLIPS_OUT, 'tooling_nk', ''); _os.makedirs(OUT, exist_ok=True)
NT = json.load(open(W + 'osc_notail.json'))
wl = np.arange(250.0, 1001.0, 1.0)
for m in ('TAPC', 'Bphen'):
    N = O.N_of(np.array(NT[m]['q']), wl)
    with open(OUT + 'nk_%s_no_visible_tail.csv' % m, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['wavelength_nm', 'n', 'k'])
        for a, z in zip(wl, N): w.writerow(['%.0f' % a, '%.5f' % z.real, '%.5f' % max(z.imag, 0.0)])
json.dump({k: json.load(open(W + f)) for k, f in [('dispersions', 'osc.json'), ('npb_per_piece', 'osc_npb_sep.json'),
           ('no_visible_tail', 'osc_notail.json'), ('si_pieces', 'si_id.json'), ('glass_global', 'glass_global.json')]},
          open(OUT + 'fit_parameters.json', 'w'), indent=1)
print('ok')

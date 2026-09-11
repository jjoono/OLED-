"""Structural self-check of the generated .mod files."""
import re, os

def check(path):
    t = open(path, 'rb').read().decode('utf-8-sig')
    fb = re.search(r"start_Gen-Osc Fit Parms\r\n(.*?)\r\n\t\tend_Gen-Osc Fit Parms", t, re.S).group(1)
    L = [x.lstrip('\t') for x in fb.split('\r\n')]
    n_osc = int(L[0].split('\t')[0])                       # first field = oscillator count
    einf = float(L[0].split('\t')[1])
    types = [x.strip().strip("'") for x in L if x.startswith("'")]
    n_par = sum(1 for x in L
                if any(("'%s'" % nm) in x for nm in
                       ('Resistivity (Ohm\u00b7cm)', 'Scat. Time (fs)', 'Amp', 'Br', 'Eo', 'Eg')))
    gb = re.search(r"start_Gen-Osc Grade Parms\r\n(.*?)\r\n\t\tend_Gen-Osc Grade Parms", t, re.S).group(1)
    n_grp = sum(1 for x in gb.split('\r\n') if x.strip() == 'F')
    pole = re.search(r"([\d.eE+-]+)\t[TF]\t-1000\.0\t1000\.0\tF\t'UV Pole Amp\.'", t).group(1)
    d = re.search(r"^\t\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'Thickness # 2'", t, re.M).group(1)
    rg = re.search(r"^\t([\d.eE+-]+)\t[TF]\t0\.0\t500\.0\tF\t'Roughness'", t, re.M).group(1)
    ok = (n_osc == len(types)) and (n_par == n_grp) and (float(pole) == 0.0)
    return dict(name=os.path.basename(path), n_osc=n_osc, types='+'.join(types), einf=einf,
                n_par=n_par, n_grp=n_grp, pole=float(pole), d=float(d)/10, rg=float(rg)/10, ok=ok)

for tag, folder in [('v2  Drude+TL', _os.path.join(ELLIPS_OUT, r'mod_files_v2')),
                    ('v3  Drude+TL+TL', _os.path.join(ELLIPS_OUT, r'mod_files_v3'))]:
    print('===== %s =====' % tag)
    print('%-26s %4s %-28s %7s %5s %5s %6s %7s %6s %s'
          % ('file', 'osc', 'types', 'Einf', 'par', 'grp', 'pole', 'd(nm)', 'rg', 'OK'))
    for f in sorted(os.listdir(folder)):
        if not f.endswith('.mod'):
            continue
        r = check(os.path.join(folder, f))
        print('%-26s %4d %-28s %7.3f %5d %5d %6.1f %7.1f %6.1f %s'
              % (r['name'][:26], r['n_osc'], r['types'][:28], r['einf'], r['n_par'],
                 r['n_grp'], r['pole'], r['d'], r['rg'], 'OK' if r['ok'] else 'FAIL'))
    print()

import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

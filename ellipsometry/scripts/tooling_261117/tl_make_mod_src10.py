"""Src10_10nm .mod files: NPB Gen-Osc on Si, as tl_make_mod.py writes the tooling
pieces, at the measured oxide (33.68 A) and in the lab template's NTVE 0 A
convention.  The material is not given; NPB is what tl_src10.py identifies.

SRC10_XLSX=/path/to/Src10_10nm.xlsx (default $ELLIPS_DATA/Src10_10nm.xlsx)."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_make_mod as MM, tl_ce_genosc as C
import numpy as np

W = _os.path.join(ELLIPS_OUT, '')
SRC10 = _os.environ.get('SRC10_XLSX', _os.path.join(ELLIPS_DATA, 'Src10_10nm.xlsx'))
OUTDIR = _os.path.join(ELLIPS_OUT, 'tooling_mod', 'Src10', '')
SHEET = 'Sheet1'


def main():
    L.XL = SRC10                      # every L.se() in this process reads the Src10 export
    _os.makedirs(OUTDIR, exist_ok=True)
    tpl = open(C.TEMPLATE, 'rb').read().decode('utf-8-sig')
    rec = json.load(open(W + 'ce_genosc.json'))['NPB']['pick']
    out = {}
    print('%-22s %8s %8s %8s  %s' % ('oxide', 'd (nm)', 'MSE fit', 'MSE file', 'check'))
    for lab, dox in (('NTVE 33.68 A', C.D_OX), ('NTVE 0 A (template)', 0.0)):
        d, mse, _ = C.fit_thickness(SHEET, np.array(rec['p']), dox, False)
        fn = 'Src10_NPB_%.1fnm_on_NTVE%sA_Si_GenOsc.mod' % (d, '%.1f' % (dox * 10) if dox else '0')
        MM.write_si(tpl, 'NPB', rec, d, OUTDIR + fn, dox, 0.0)
        mf, m = MM.mse_from_file(OUTDIR + fn, SHEET)
        ok = m['structure_ok'] and abs(mf - mse) < 0.05
        out[lab] = dict(d=d, mse=mse, file=fn)
        print('%-22s %8.2f %8.2f %8.2f  %s  %s' % (lab, d, mse, mf, 'OK' if ok else 'MISMATCH', fn))
    json.dump(out, open(W + 'src10_mod.json', 'w'), indent=1)
    print('->', OUTDIR)


if __name__ == '__main__':
    main()

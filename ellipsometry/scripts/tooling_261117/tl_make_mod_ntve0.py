"""Si-piece .mod files in the lab tooling template's convention, NTVE 0 A
(Bphen_S5 in the BPhen_S12_tooling snapshot fits Si with no native oxide).

Same writer and dispersions as tl_make_mod.py; only the oxide (0 A, fit off)
and the thickness/roughness refitted for it (tl_ce_thickness_ntve0.py) differ.
Use these to compare tooling factors with earlier runs of that template; the
33.68 A files fit better (lower MSE) and are the better thickness estimate."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_make_mod as MM, tl_ce_genosc as C, tl_ce_thickness as T

W = _os.path.join(ELLIPS_OUT, '')
OUTDIR = _os.path.join(ELLIPS_OUT, 'tooling_mod', 'Si_NTVE0_template_convention', '')


def main():
    _os.makedirs(OUTDIR, exist_ok=True)
    tpl = open(C.TEMPLATE, 'rb').read().decode('utf-8-sig')
    N0 = json.load(open(W + 'ce_thickness_ntve0.json'))
    print('%-5s %-6s %8s | %-40s | %8s %8s  %s' % ('sheet', 'mat', 'd (nm)', 'file', 'MSE fit', 'MSE file', 'check'))
    for sh, v in N0.items():
        name, d, rg = v['mat'], v['d_ntve0'], v['rough_ntve0']
        fn = '%s_%s_%.1fnm_on_NTVE0A_Si_GenOsc.mod' % (sh, name, d)
        MM.write_si(tpl, name, T.disp(name), d, OUTDIR + fn, 0.0, rg)
        mse_file, m = MM.mse_from_file(OUTDIR + fn, sh)
        ok = m['structure_ok'] and m['dox'] == 0.0 and abs(mse_file - v['mse_ntve0']) < 0.05
        print('%-5s %-6s %8.2f | %-40s | %8.2f %8.2f  %s' % (sh, name, d, fn[:40], v['mse_ntve0'], mse_file, 'OK' if ok else 'MISMATCH'))
    print('->', OUTDIR)


if __name__ == '__main__':
    main()

"""CompleteEASE .mod generator v7 - written from the full-range fit done in
CompleteEASE's own frame (SI_JAW + NTVE_JAW decoded out of the .mod itself),
so Generate alone reproduces the data over the model's valid range.

Parameters come straight from ce_final_result.json; no unit juggling is needed
because ce_osc.py already uses CompleteEASE's own conventions
(Drude rho[Ohm.cm]/tau[fs], TL Amp/Br/Eo/Eg, Gaussian Amp/Br/En).
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import re, os, json

TEMPLATE = _os.path.join(ELLIPS_DATA, r'VendorA_ITO_GenOSC_6functions.mod')
OUTDIR   = _os.path.join(ELLIPS_OUT, r'mod_v8')
RESULT   = _os.path.join(ELLIPS_OUT, r'ce_final_result.json')
NAME = {'#1': 'VendorA_ITO_1', '#2': 'VendorA_ITO_2', '#3': 'VendorB_IZO_1',
        '#4': 'VendorB_IZO_2', '#5': 'VendorA_ITO_2pctO2'}
TAIL = "\tF\tF\t0.0\t0.0\t100000.0\tF\t100.0\t"

def num(x): return repr(float(x))
def pline(v, fit, lo, hi, nm):
    return "\t\t\t%s\t%s\t%s\t%s\tF\t'%s'%s" % (v, 'T' if fit else 'F', lo, hi, nm, TAIL)
def gline(lo, hi, i):
    return "\t\t\t0.0\tF\t%s\t%s\tF\t'Grade %d'%s" % (lo, hi, i, TAIL)

def blocks(p):
    """p = [d, rough, dth, Einf, rho, tau, TL_A, TL_Br, TL_Eo, TL_Eg, G_A, G_Br, G_En]"""
    osc = [
        ("'Drude(RT)'", False, [
            ('Resistivity (Ohm\u00b7cm)', num(p[4]), True, '1.0E-5', '1.0E-1', '1.0E-9', '1000.0'),
            ('Scat. Time (fs)',           num(p[5]), True, '0.5',    '30.0',   '1.0E-6', '1000.0')]),
        ("'Tauc-Lorentz'", True, [
            ('Amp', num(p[6]), True, '10.0',   '600.0',  '0.001',  '1000.0'),
            ('Br',  num(p[7]), True, '0.02',   '6.00',   '1.0E-4', '1000.0'),
            ('Eo',  num(p[8]), True, '3.20',   '5.50',   '1.0E-4', '15.0'),
            ('Eg',  num(p[9]), True, '2.80',   '4.20',   '0.0',    '15.0')]),
        ("'Gaussian'", False, [
            ('Amp',  num(p[10]), True,  '0.0',     '5.0',    '-100.0',  '1000.0'),
            ('iAmp', '0.0',      False, '-1000.0', '1000.0', '-1000.0', '1000.0'),
            ('Br',   num(p[11]), True,  '0.10',    '5.00',   '0.0',     '100.0'),
            ('En',   num(p[12]), True,  '0.80',    '3.50',   '1.0E-8',  '15.0')]),
    ]
    fit = ["\t\t\t%d\t%s\tT\t0.5\t5.0\tF\t'Einf'%s" % (len(osc), num(p[3]), TAIL)]
    grade, npar = [], 0
    for typ, common_eg, pl in osc:
        fit.append('\t\t\t' + typ + '\t')
        for nm, v, f, lo, hi, glo, ghi in pl:
            fit.append(pline(v, f, lo, hi, nm))
            grade.append('\t\t\tF\t'); grade += [gline(glo, ghi, j) for j in (1, 2, 3)]
            npar += 1
        if common_eg: fit.append('\t\t\tF\t')
    fit.append('\t\t\t'); grade.append('\t\t\t')
    return '\r\n'.join(fit), '\r\n'.join(grade), len(osc), npar

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    tpl = open(TEMPLATE, 'rb').read().decode('utf-8-sig')
    R = json.load(open(RESULT))
    print('%-16s %6s %8s %7s %7s' % ('sample', 'MSE', 'd(nm)', 'rough', 'dth'))
    for sh, rec in R.items():
        p = rec['p']
        fb, gb, n_osc, n_par = blocks(p)
        t = tpl
        t = re.sub(r"(start_Gen-Osc Fit Parms\r\n).*?(\r\n\t\tend_Gen-Osc Fit Parms)",
                   lambda m: m.group(1) + fb + m.group(2), t, count=1, flags=re.S)
        t = re.sub(r"(start_Gen-Osc Grade Parms\r\n).*?(\r\n\t\tend_Gen-Osc Grade Parms)",
                   lambda m: m.group(1) + gb + m.group(2), t, count=1, flags=re.S)
        t = re.sub(r"^\t[^\t\n]+(\t[TF]\t)-5\.0\t5\.0(\tF\t'Angle Offset')",
                   lambda m: '\t' + num(p[2]) + m.group(1) + '-2.0\t2.0' + m.group(2), t, count=1, flags=re.M)
        t = re.sub(r"^\t[^\t\n]+(\t[TF]\t)0\.0\t500\.0(\tF\t'Roughness')",
                   lambda m: '\t' + num(p[1]*10) + m.group(1) + '0.0\t150.0' + m.group(2), t, count=1, flags=re.M)
        t = re.sub(r"^\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Thickness # 2')",
                   lambda m: '\t\t' + num(p[0]*10) + '\tT\t250.0\t900.0' + m.group(1), t, count=1, flags=re.M)
        t = t.replace('263.42946105543194\tT\t', '0.0\tF\t', 1)
        t = t.replace('11.015243042816298\tT\t', '11.0\tF\t', 1)
        out = os.path.join(OUTDIR, NAME[sh] + '_GenOsc_v8.mod')
        open(out, 'w', encoding='utf-8-sig', newline='').write(t)
        print('%-16s %6.2f %8.2f %7.2f %+7.3f' % (NAME[sh], rec['mse'], p[0], p[1], p[2]))
    print('saved ->', OUTDIR)

if __name__ == '__main__':
    main()

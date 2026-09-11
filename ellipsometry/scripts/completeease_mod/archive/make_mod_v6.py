"""CompleteEASE .mod generator v6.

The .SE data actually spans ~190-1800 nm; the xlsx extract I fitted from was
truncated at 1078 nm.  So the Drude IS determined by the measurement, and my
starting Drude (screened plasma energy 0.84 eV = 1484 nm) is far too weak --
CompleteEASE's own fit of #1 gives 1.29 eV = 958 nm, driven by the real NIR.

v6 therefore:
  * frees Drude Resistivity AND Scat. Time (v5 wrongly froze tau)
  * frees Tauc-Lorentz Eo so the deep UV (down to 190 nm) can adjust
  * starts every sample's Drude amplitude scaled by 2.40, the factor
    CompleteEASE needed for #1 -- all five share the same NIR blind spot
  * KEEPS Gaussian Amp >= 0, the one constraint that makes k<0 impossible
"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')   # measurement exports (.xlsx), CompleteEASE .mod
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')    # fitted results, figures, intermediates

import re, os, math, json

HBAR_EVS = 6.582119569e-16
HBAR_EVFS = 0.6582119569
EPS0 = 8.8541878128e-12
GD = 0.10
DRUDE_BOOST = 1.674 / 0.69841     # from CompleteEASE's own NIR-informed fit of #1

TEMPLATE = _os.path.join(ELLIPS_DATA, r'VendorA_ITO_GenOSC_6functions.mod')
OUTDIR   = _os.path.join(ELLIPS_OUT, r'mod_v6')
PARAMS   = _os.path.join(ELLIPS_OUT, r'genosc_params.json')
GEO = {'1': (51.0, 3.0, 0.00), '2': (51.0, 3.0, -0.20), '3': (42.0, 2.0, 0.18),
       '4': (42.0, 2.0, -0.14), '5': (42.0, 2.0, 0.00)}
NAME = {'1': 'VendorA_ITO_1', '2': 'VendorA_ITO_2', '3': 'VendorB_IZO_1',
        '4': 'VendorB_IZO_2', '5': 'VendorA_ITO_2pctO2'}

TAIL = "\tF\tF\t0.0\t0.0\t100000.0\tF\t100.0\t"

def pline(val, fit, lo, hi, name):
    return "\t\t\t%s\t%s\t%s\t%s\tF\t'%s'%s" % (val, 'T' if fit else 'F', lo, hi, name, TAIL)

def gline(lo, hi, i):
    return "\t\t\t0.0\tF\t%s\t%s\tF\t'Grade %d'%s" % (lo, hi, i, TAIL)

def num(x):
    return repr(float(x))

def build(p):
    """p = [einf, A_Drude, A_TL, E0_TL, C_TL, Eg_TL, A_G, Ec_G, sigma_G]"""
    tau = HBAR_EVFS / GD
    A_d = p[1] * DRUDE_BOOST           # NIR was missing from the data I fitted
    rho = HBAR_EVS**2 / (A_d * EPS0 * (tau * 1e-15)) * 100.0
    brG = p[8] * 2.0 * math.sqrt(math.log(2.0))
    eg  = p[5]

    # (name, value, fit?, lo, hi)   grade groups keep the CompleteEASE defaults
    osc = [
        ("'Drude(RT)'", False, [
            ('Resistivity (Ohm\u00b7cm)', num(rho), True, '1.0E-4', '5.0E-2', '1.0E-9', '1000.0'),
            ('Scat. Time (fs)',           num(tau), True,  '1.0',    '20.0',   '1.0E-6', '1000.0'),
        ]),
        ("'Tauc-Lorentz'", True, [
            ('Amp', num(p[2]), True,  '20.0',  '400.0', '0.001',  '1000.0'),
            ('Br',  num(p[4]), True, '0.10', '3.00', '1.0E-4', '1000.0'),
            ('Eo',  num(p[3]), True, '3.50', '5.50', '1.0E-4', '15.0'),
            ('Eg',  num(eg),   True, num(round(eg - 0.25, 4)), num(round(eg + 0.25, 4)), '0.0', '15.0'),
        ]),
        ("'Gaussian'", False, [
            ('Amp',  num(p[6]), True,  '0.0',     '1.0',    '-100.0',  '1000.0'),
            ('iAmp', '0.0',     False, '-1000.0', '1000.0', '-1000.0', '1000.0'),
            ('Br',   num(brG),  True,  '0.20',    '3.00',   '0.0',     '100.0'),
            ('En',   num(p[7]), True,  '1.00',    '3.40',   '1.0E-8',  '15.0'),
        ]),
    ]

    fit = ["\t\t\t%d\t%s\tT\t1.0\t5.0\tF\t'Einf'%s" % (len(osc), num(p[0]), TAIL)]
    grade, npar = [], 0
    for typ, common_eg, plist in osc:
        fit.append('\t\t\t' + typ + '\t')
        for nm, v, f, lo, hi, glo, ghi in plist:
            fit.append(pline(v, f, lo, hi, nm))
            grade.append('\t\t\tF\t')
            grade += [gline(glo, ghi, j) for j in (1, 2, 3)]
            npar += 1
        if common_eg:
            fit.append('\t\t\tF\t')          # Common Eg flag (Tauc-Lorentz only)
    fit.append('\t\t\t')
    grade.append('\t\t\t')
    return '\r\n'.join(fit), '\r\n'.join(grade), len(osc), npar, rho, tau

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    tpl = open(TEMPLATE, 'rb').read().decode('utf-8-sig')
    R = json.load(open(PARAMS))
    print('%-16s %3s %5s %10s %8s  %s' % ('sample', 'osc', 'par', 'rho', 'tau', 'geometry'))
    for s in '12345':
        p = R[s]['p']
        fitblk, gradeblk, n_osc, n_par, rho, tau = build(p)
        d_nm, rg_nm, dth = GEO[s]
        t = tpl
        t = re.sub(r"(start_Gen-Osc Fit Parms\r\n).*?(\r\n\t\tend_Gen-Osc Fit Parms)",
                   lambda m: m.group(1) + fitblk + m.group(2), t, count=1, flags=re.S)
        t = re.sub(r"(start_Gen-Osc Grade Parms\r\n).*?(\r\n\t\tend_Gen-Osc Grade Parms)",
                   lambda m: m.group(1) + gradeblk + m.group(2), t, count=1, flags=re.S)
        t = re.sub(r"^\t[^\t\n]+(\t[TF]\t)-5\.0\t5\.0(\tF\t'Angle Offset')",
                   lambda m: '\t' + num(dth) + m.group(1) + '-1.0\t1.0' + m.group(2),
                   t, count=1, flags=re.M)
        t = re.sub(r"^\t[^\t\n]+(\t[TF]\t)0\.0\t500\.0(\tF\t'Roughness')",
                   lambda m: '\t' + num(rg_nm * 10) + m.group(1) + '0.0\t100.0' + m.group(2),
                   t, count=1, flags=re.M)
        t = re.sub(r"^\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Thickness # 2')",
                   lambda m: '\t\t' + num(d_nm * 10) + '\tT\t300.0\t800.0' + m.group(1),
                   t, count=1, flags=re.M)
        t = t.replace('263.42946105543194\tT\t', '0.0\tF\t', 1)
        t = t.replace('11.015243042816298\tT\t', '11.0\tF\t', 1)
        out = os.path.join(OUTDIR, NAME[s] + '_GenOsc_v6.mod')
        open(out, 'w', encoding='utf-8-sig', newline='').write(t)
        print('%-16s %3d %5d %10.4g %8.3f  d=%.0fA rough=%.0fA dth=%+.2f'
              % (NAME[s], n_osc, n_par, rho, tau, d_nm * 10, rg_nm * 10, dth))
    print('saved ->', OUTDIR)

if __name__ == '__main__':
    main()

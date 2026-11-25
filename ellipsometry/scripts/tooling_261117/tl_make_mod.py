"""Write CompleteEASE .mod files for the tooling set, then read them back.

Structure comes from the lab's own Gen-Osc file (ADNC_ITO_GenOSC_v2_claude.mod,
saved by CompleteEASE, model version 2: SI_JAW / NTVE_JAW / Gen-Osc).  Only the
Gen-Osc layer is rewritten, with the syntax the 260819 generators established
by trial (README_260820_session.md section 2):

  * Fit Parms line 0 = <oscillator count> TAB <Einf> ...
  * Gaussian = four parameter lines Amp, iAmp, Br, En (iAmp hidden in the GUI),
    no trailing flag line (that is Tauc-Lorentz only)
  * Grade Parms holds one 4-line group per oscillator PARAMETER
  * the UV pole lives in 'Gen-Osc Permanent Poles', not in Fit Parms
  * thicknesses in Angstrom; fit flags on parameter lines are honoured

Einf is written fixed at 1.  A parameter whose 90 % CI exceeds 50 % of its value,
or that sits on a bound, is written with its fit flag OFF so a refit in
CompleteEASE cannot slide it.  Every file is parsed back, the stack rebuilt from
what the FILE says, and its MSE against the data must equal the fit's.

Glass pieces need a different substrate and the back surface, neither of which
the template carries; see make_glass()."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, re, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_genosc as G, tl_ce_genosc as C
import jnk_model as M

W = _os.path.join(ELLIPS_OUT, '')
OUTDIR = _os.path.join(ELLIPS_OUT, 'tooling_mod', '')
TAIL = "\tF\tF\t0.0\t0.0\t100000.0\tF\t100.0\t"
num = lambda x: repr(float(x))
pl = lambda v, f, lo, hi, nm: "\t\t\t%s\t%s\t%s\t%s\tF\t'%s'%s" % (v, 'T' if f else 'F', lo, hi, nm, TAIL)
gl = lambda lo, hi, i: "\t\t\t0.0\tF\t%s\t%s\tF\t'Grade %d'%s" % (lo, hi, i, TAIL)


def lock_flags(rec):
    """True = leave free in CompleteEASE.  Off when the CI is wide or a bound is hit."""
    p, ci = np.array(rec['p']), np.array(rec['ci'])
    free = [(ci[i] <= 0.5 * abs(p[i])) and (i not in rec['at_bound']) for i in range(len(p))]
    return free


def genosc_blocks(name, rec):
    p = np.array(rec['p']); free = lock_flags(rec)
    lo, hi = G.bounds(name, (len(p) - 2) // 3)
    ng = (len(p) - 2) // 3
    fit = ["\t\t\t%d\t%s\tF\t0.0\t10.0\tF\t'Einf'%s" % (ng, num(1.0), TAIL)]
    grade = []
    for g in range(ng):
        a, b, e = p[2 + 3 * g: 5 + 3 * g]; fa, fb, fe = free[2 + 3 * g: 5 + 3 * g]
        fit.append("\t\t\t'Gaussian'\t")
        fit.append(pl(num(a), fa, num(lo[2 + 3 * g]), num(hi[2 + 3 * g]), 'Amp'))
        fit.append(pl('0.0', False, '-1000.0', '1000.0', 'iAmp'))
        fit.append(pl(num(b), fb, num(lo[3 + 3 * g]), num(hi[3 + 3 * g]), 'Br'))
        fit.append(pl(num(e), fe, num(lo[4 + 3 * g]), num(hi[4 + 3 * g]), 'En'))
        for glo, ghi in (('-100.0', '1000.0'), ('-1000.0', '1000.0'), ('0.0', '100.0'), ('1.0E-8', '15.0')):
            grade.append('\t\t\tF\t'); grade += [gl(glo, ghi, j) for j in (1, 2, 3)]
    fit.append('\t\t\t'); grade.append('\t\t\t')
    return '\r\n'.join(fit), '\r\n'.join(grade), (p[0], free[0]), (p[1], free[1])


def write_si(tpl, name, rec, d_nm, out_path, d_ox_nm, rough_nm=0.0):
    fb, gb, (pa, fpa), (pe, fpe) = genosc_blocks(name, rec)
    t = tpl
    t = re.sub(r"(start_Gen-Osc Fit Parms\r\n).*?(\r\n\t\tend_Gen-Osc Fit Parms)",
               lambda m: m.group(1) + fb + m.group(2), t, count=1, flags=re.S)
    t = re.sub(r"(start_Gen-Osc Grade Parms\r\n).*?(\r\n\t\tend_Gen-Osc Grade Parms)",
               lambda m: m.group(1) + gb + m.group(2), t, count=1, flags=re.S)
    t = re.sub(r"^\t\t\t[^\t\n]+\t[TF](\t-1000\.0\t1000\.0\tF\t'UV Pole Amp\.')",
               lambda m: '\t\t\t%s\t%s%s' % (num(pa), 'T' if fpa else 'F', m.group(1)), t, count=1, flags=re.M)
    t = re.sub(r"^\t\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'UV Pole En\.')",
               lambda m: '\t\t\t%s\t%s\t6.0\t15.0%s' % (num(pe), 'T' if fpe else 'F', m.group(1)), t, count=1, flags=re.M)
    t = re.sub(r"^\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Thickness # 2')",
               lambda m: '\t\t%s\tT\t%s\t%s%s' % (num(d_nm * 10), num(max(d_nm - 10, 0.5) * 10), num((d_nm + 10) * 10), m.group(1)),
               t, count=1, flags=re.M)
    t = re.sub(r"^\t\t[^\t\n]+\t[TF](\t[^\t\n]+\t[^\t\n]+\tT\t'Native Oxide')",
               lambda m: '\t\t%s\tF%s' % (num(d_ox_nm * 10), m.group(1)), t, count=1, flags=re.M)
    t = re.sub(r"^\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Angle Offset')",
               lambda m: '\t0.0\tF\t-2.0\t2.0' + m.group(1), t, count=1, flags=re.M)
    t = re.sub(r"^\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Roughness')",
               lambda m: '\t%s\t%s\t0.0\t100.0%s' % (num(rough_nm * 10), 'T' if rough_nm > 0.05 else 'F', m.group(1)),
               t, count=1, flags=re.M)
    open(out_path, 'w', encoding='utf-8-sig', newline='').write(t)
    return t


# ------------------------------------------------------------------ read back
def parse(path):
    t = open(path, 'rb').read().decode('utf-8-sig')
    fb = re.search(r"start_Gen-Osc Fit Parms\r\n(.*?)\r\n\t\tend_Gen-Osc Fit Parms", t, re.S).group(1)
    Ls = [x.strip('\t') for x in fb.split('\r\n')]
    n_osc = int(Ls[0].split('\t')[0]); einf = float(Ls[0].split('\t')[1])
    oscs = []; cur = None
    for ln in Ls[1:]:
        if ln.startswith("'") and ln.endswith("'"):
            cur = [ln.strip("'"), {}]; oscs.append(cur)
        elif "\t'" in ln and cur is not None:
            f = ln.split('\t'); cur[1][f[5].strip("'")] = float(f[0])
    gb = re.search(r"start_Gen-Osc Grade Parms\r\n(.*?)\r\n\t\tend_Gen-Osc Grade Parms", t, re.S).group(1)
    n_grp = sum(1 for x in gb.split('\r\n') if x.strip() == 'F')
    pa = float(re.search(r"\n\t\t\t([-\d.eE+]+)\t[TF]\t-1000\.0\t1000\.0\tF\t'UV Pole Amp\.'", t).group(1))
    pe = float(re.search(r"\n\t\t\t([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'UV Pole En\.'", t).group(1))
    d = float(re.search(r"\n\t\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'Thickness # 2'", t).group(1)) / 10
    dox = float(re.search(r"\n\t\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tT\t'Native Oxide'", t).group(1)) / 10
    dth = float(re.search(r"\n\t([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'Angle Offset'", t).group(1))
    rg = float(re.search(r"\n\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'Roughness'", t).group(1)) / 10
    ok = (n_osc == len(oscs)) and (n_grp == 4 * len(oscs)) and all(o[0] == 'Gaussian' and len(o[1]) == 4 for o in oscs)
    return dict(einf=einf, oscs=oscs, pole=(pa, pe), d=d, dox=dox, dth=dth, rough=rg, structure_ok=ok, n_grp=n_grp)


def mse_from_file(path, sh):
    m = parse(path)
    E_osc = []
    for typ, ps in m['oscs']:
        E_osc += [ps['Amp'], ps['Br'], ps['En']]
    p = np.array([m['pole'][0], m['pole'][1]] + E_osc)
    assert abs(m['einf'] - G.EINF) < 1e-12
    data = C.data_for([sh])
    rough = m['rough'] > 0
    x = np.concatenate([p, [m['d']] + ([m['rough']] if rough else [])])
    r = C.residual(x, data, len(p), None, rough, m['dox'])
    return M.mse_ncs(r, 2 if rough else 1), m


def main():
    _os.makedirs(OUTDIR, exist_ok=True)
    tpl = open(C.TEMPLATE, 'rb').read().decode('utf-8-sig')
    CE = json.load(open(W + 'ce_genosc.json'))
    TH = json.load(open(W + 'ce_thickness.json'))
    print('%-5s %-6s %8s | %-44s | %8s %8s  %s' % ('sheet', 'mat', 'd (nm)', 'file', 'MSE fit', 'MSE file', 'check'))
    for sh, v in TH.items():
        name, d, mse_fit, rg = v['mat'], v['d'], v['mse'], v['rough']
        fn = '%s_%s_%.1fnm_on_NTVE%.1fA_Si_GenOsc.mod' % (sh, name, d, C.D_OX * 10)
        path = OUTDIR + fn
        write_si(tpl, name, CE[name][v['dispersion']], d, path, C.D_OX, rg)
        mse_file, m = mse_from_file(path, sh)
        ok = m['structure_ok'] and abs(mse_file - mse_fit) < 0.05
        print('%-5s %-6s %8.2f | %-44s | %8.2f %8.2f  %s' % (sh, name, d, fn[:44], mse_fit, mse_file, 'OK' if ok else 'MISMATCH'))
    print('->', OUTDIR)


if __name__ == '__main__':
    main()

"""Glass-piece .mod files on the lab's own tooling template (replaces the
inferred-format tl_make_mod_glass.py output).

The structure is the CompleteEASE-saved Bphen_S5 model pulled out of the
BPhen_S12_tooling snapshot (tl_s5_extract.py):
  Layer0  EMA(Si_JAW, glass Gen-Osc), EMA % (Mat 2) -> 100 = glass
  Layer1  NTVE_JAW, 0 A, fit off
  Layer2  the organic film, Gen-Osc
  back reflections ON (the flag before '# Back Reflections', T in the template),
  '% 1st Reflection' = the piece's fitted collected fraction, fit ON
Only values change, plus: the template's two-sample MultiSamp block is swapped
for the single-sample one of the ADNC file (both CompleteEASE-saved), and the
glass material is renamed consistently from 'no5_glass_sub_genosc' because its
values are this set's glass, not no5.

The glass material carries this set's glass as Einf 1 + UV pole.  The
template's one Tauc-Lorentz is kept (same block structure) but parked: Amp at
its lower bound 0.001 and Eg 6.0 / Eo 7.0 eV, above the measured range
(<= 5.04 eV), so the glass has k = 0 everywhere in the data.  Left at the no5
values (Eg 3.27 eV), even Amp 0.001 gives k ~ 1e-6 below 379 nm, which over
the 1 mm substrate takes ~16 % off the 300 nm back-surface beam."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, re, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_make_mod as MM, tl_ce_genosc as C, tl_genosc as GO, tl_glass as GL, tl_glass_global as GG, tl_load as L
import jnk_model as M, ce_osc as osc

W = _os.path.join(ELLIPS_OUT, '')
S5 = _os.path.join(ELLIPS_OUT, 'S5_template.mod')
OUTDIR = _os.path.join(ELLIPS_OUT, 'tooling_mod', 'glass', '')
GLASS = 'tooling2611_glass'
num = MM.num
HC = 1239.841984
TL_EO, TL_EG = 7.0, 6.0          # parked above the data (max 5.04 eV); see docstring


def sub1(rx, rep, t, flags=re.M):
    t2, n = re.subn(rx, rep, t, count=1, flags=flags)
    assert n == 1, rx
    return t2


def write(sh, mat, rec, d_nm, f, glass, path):
    t = open(S5, 'rb').read().decode('utf-8-sig')
    adnc = open(C.TEMPLATE, 'rb').read().decode('utf-8-sig')
    single = re.search(r"\tstart_MultiSamp\r\n.*?\tend_MultiSamp", adnc, re.S).group(0)
    t = sub1(r"\tstart_MultiSamp\r\n.*?\tend_MultiSamp", lambda m: single, t, re.S)
    t = t.replace('no5_glass_sub_genosc', GLASS)
    # ---- model parameters
    t = sub1(r"^\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Angle Offset')", lambda m: '\t0.0\tF\t-2.0\t2.0' + m.group(1), t)
    t = sub1(r"^\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Roughness')", lambda m: '\t0.0\tF\t0.0\t100.0' + m.group(1), t)
    t = sub1(r"\r\n\t[TF]\t(\r\n\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+\tF\t'# Back Reflections')", lambda m: '\r\n\tT\t' + m.group(1), t, 0)
    t = sub1(r"^\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'# Back Reflections')", lambda m: '\t5.0\tF\t0.0\t20.0' + m.group(1), t)
    t = sub1(r"^\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'% 1st Reflection')", lambda m: '\t%s\tT\t0.0\t100.0%s' % (num(100 * f), m.group(1)), t)
    # ---- substrate: 100 % glass, glass = Einf 1 + UV pole (+ nil TL)
    L0 = re.search(r"start_Layer0.*?end_Layer0", t, re.S).group(0); L0n = L0
    L0n = sub1(r"(start_%s Fit Parms\r\n\t+)1\t[^\t]+\t" % GLASS, lambda m: m.group(1) + '1\t1.0\t', L0n, 0)
    L0n = sub1(r"^(\t+)[^\t\n]+(\t[TF]\t0\.001\t1000\.0\tF\t'Amp')", lambda m: m.group(1) + '0.001' + m.group(2), L0n)
    L0n = sub1(r"^(\t+)[^\t\n]+(\t[TF]\t1\.0E-4\t15\.0\tF\t'Eo')", lambda m: m.group(1) + num(TL_EO) + m.group(2), L0n)
    L0n = sub1(r"^(\t+)[^\t\n]+(\t[TF]\t0\.0\t15\.0\tF\t'Eg')", lambda m: m.group(1) + num(TL_EG) + m.group(2), L0n)
    L0n = sub1(r"^(\t+)[^\t\n]+(\t[TF]\t-1000\.0\t1000\.0\tF\t'UV Pole Amp\.')", lambda m: m.group(1) + num(glass['pole_amp']) + m.group(2), L0n)
    L0n = sub1(r"^(\t+)[^\t\n]+(\t[TF]\t[^\t\n]+\t[^\t\n]+\tF\t'UV Pole En\.')", lambda m: m.group(1) + num(glass['pole_en']) + m.group(2), L0n)
    L0n = sub1(r"^(\t+)[^\t\n]+(\t[TF]\t-1000\.0\t1000\.0\tF\t'IR Pole Amp\.')", lambda m: m.group(1) + '0.0' + m.group(2), L0n)
    L0n = sub1(r"^(\t+)[^\t\n]+\t[TF](\t-10\.0\t100\.0\tF\t'EMA % \(Mat 2\)')", lambda m: m.group(1) + '100.0\tF' + m.group(2), L0n)
    t = t.replace(L0, L0n)
    # ---- NTVE off
    t = sub1(r"^\t\t[^\t\n]+\t[TF](\t[^\t\n]+\t[^\t\n]+\tF\t'Thickness # 1')", lambda m: '\t\t0.0\tF' + m.group(1), t)
    # ---- film
    L2 = re.search(r"start_Layer2.*?end_Layer2", t, re.S).group(0); L2n = L2
    fb, gb, (pa, fpa), (pe, fpe) = MM.genosc_blocks(mat, rec)
    L2n = sub1(r"(start_Gen-Osc Fit Parms\r\n).*?(\r\n\t\tend_Gen-Osc Fit Parms)", lambda m: m.group(1) + fb + m.group(2), L2n, re.S)
    L2n = sub1(r"(start_Gen-Osc Grade Parms\r\n).*?(\r\n\t\tend_Gen-Osc Grade Parms)", lambda m: m.group(1) + gb + m.group(2), L2n, re.S)
    L2n = sub1(r"^\t\t\t[^\t\n]+\t[TF](\t-1000\.0\t1000\.0\tF\t'UV Pole Amp\.')", lambda m: '\t\t\t%s\t%s%s' % (num(pa), 'T' if fpa else 'F', m.group(1)), L2n)
    L2n = sub1(r"^\t\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'UV Pole En\.')", lambda m: '\t\t\t%s\t%s\t6.0\t15.0%s' % (num(pe), 'T' if fpe else 'F', m.group(1)), L2n)
    L2n = sub1(r"^\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Thickness # 2')",
               lambda m: '\t\t%s\tT\t%s\t%s%s' % (num(d_nm * 10), num(max(d_nm - 10, 0.0) * 10), num((d_nm + 10) * 10), m.group(1)), L2n)
    t = t.replace(L2, L2n)
    open(path, 'w', encoding='utf-8-sig', newline='').write(t)


def parse(path):
    t = open(path, 'rb').read().decode('utf-8-sig')
    g = lambda rx, s: float(re.search(rx, s).group(1))
    L0 = re.search(r"start_Layer0.*?end_Layer0", t, re.S).group(0)
    L2 = re.search(r"start_Layer2.*?end_Layer2", t, re.S).group(0)
    gl = dict(einf=g(r"start_%s Fit Parms\r\n\t+1\t([-\d.eE+]+)" % GLASS, L0),
              tl=[g(r"\n\t+([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'%s'" % k, L0) for k in ('Amp', 'Br', 'Eo', 'Eg')],
              pa=g(r"\n\t+([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'UV Pole Amp\.'", L0),
              pe=g(r"\n\t+([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'UV Pole En\.'", L0),
              ir=g(r"\n\t+([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'IR Pole Amp\.'", L0),
              ema=g(r"\n\t+([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'EMA % \(Mat 2\)'", L0))
    fb = re.search(r"start_Gen-Osc Fit Parms\r\n(.*?)\r\n\t\tend_Gen-Osc Fit Parms", L2, re.S).group(1)
    Ls = [x.strip('\t') for x in fb.split('\r\n')]; oscs = []; cur = None
    for ln in Ls[1:]:
        if ln.startswith("'") and ln.endswith("'"): cur = {}; oscs.append(cur)
        elif "\t'" in ln and cur is not None:
            fl = ln.split('\t'); cur[fl[5].strip("'")] = float(fl[0])
    p = [g(r"\n\t\t\t([-\d.eE+]+)\t[TF]\t-1000\.0\t1000\.0\tF\t'UV Pole Amp\.'", L2),
         g(r"\n\t\t\t([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'UV Pole En\.'", L2)] + sum([[o['Amp'], o['Br'], o['En']] for o in oscs], [])
    return dict(glass=gl, p=np.array(p), einf_film=float(Ls[0].split('\t')[1]),
                d=g(r"\n\t\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'Thickness # 2'", t) / 10,
                ntve=g(r"\n\t\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'Thickness # 1'", t),
                f=g(r"\n\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'% 1st Reflection'", t) / 100,
                back_on=re.search(r"\r\n\t([TF])\t\r\n\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+\tF\t'# Back Reflections'", t).group(1) == 'T',
                n_samples=int(re.search(r"start_MultiSamp\r\n\t\t(\d+)\t", t).group(1)))


def mse_file(path, sh):
    m = parse(path); gl = m['glass']
    wl, P, D, ang = L.se(sh, 300.0, 1000.0, 3); meas = L.ncs(P, D)
    E = HC / wl
    a, b = osc.tauc_lorentz(E, *gl['tl'])
    assert np.all(b == 0.0), 'glass absorbs inside the data range - the fit assumed it does not'
    ng = np.sqrt(gl['einf'] + gl['pa'] / (gl['pe']**2 - E**2) - gl['ir'] / E**2 + a) + 0j
    assert abs(gl['ema'] - 100) < 1e-9 and m['ntve'] == 0.0 and m['back_on'] and m['n_samples'] == 0
    N, C_, S_ = GG.norm(*GL.ncs_glass(wl, GO.N_of(m['p'], wl), m['d'], ng, ang, m['f']))
    r = np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    return M.mse_ncs(r, 2 + 2 / 8), m


def main():
    _os.makedirs(OUTDIR, exist_ok=True)
    CE = json.load(open(W + 'ce_genosc.json')); GJ = json.load(open(W + 'ce_glass2.json'))
    print('%-5s %-6s %7s %6s | %-34s | %7s %7s  %s' % ('sheet', 'mat', 'd (nm)', '%1st', 'file', 'MSE fit', 'file', 'check'))
    for sh, mat in GJ['mat'].items():
        d, f = GJ['d'][sh], GJ['f'][sh]
        tag = ' (bare)' if d < 1.0 else ''
        fn = '%s_%s_%.1fnm_on_glass_GenOsc.mod' % (sh, mat, d)
        write(sh, mat, CE[mat]['pick'], d, f, GJ['glass_genosc'], OUTDIR + fn)
        mf, m = mse_file(OUTDIR + fn, sh)
        ok = abs(mf - GJ['mse'][sh]) < 0.3
        print('%-5s %-6s %7.2f %6.1f | %-34s | %7.2f %7.2f  %s%s' % (sh, mat, d, 100 * f, fn[:34], GJ['mse'][sh], mf, 'OK' if ok else 'MISMATCH', tag))
    print('->', OUTDIR)


if __name__ == '__main__':
    main()

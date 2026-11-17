"""Glass-piece .mod files - UNVERIFIED FORMAT, test-load one before using the rest.

None of the supplied files has a glass substrate, so two things are inferred
from the Si template rather than copied from a CompleteEASE-saved glass model:

  * the substrate: the SI_JAW table layer is replaced by a Gen-Osc layer with no
    oscillators (Einf + UV pole - a one-term Sellmeier fitted to the glass
    Cauchy), and the NTVE_JAW layer is dropped; layers renumbered, count bumped
    down, the same way make_mod_ag added layers;
  * the back surface: the T/F flag on the line before '# Back Reflections' is
    taken to be the switch, and '% 1st Reflection' is set to the fitted
    collected fraction and left free.

If CompleteEASE refuses the file, or Generate does not land near the MSE
printed here, a glass-substrate .mod saved from CompleteEASE settles both."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, re, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_make_mod as MM, tl_ce_genosc as C, tl_genosc as GO, tl_glass as GL, tl_glass_global as GG, tl_load as L
import jnk_model as M

W = _os.path.join(ELLIPS_OUT, '')
OUTDIR = _os.path.join(ELLIPS_OUT, 'tooling_mod', 'glass_UNVERIFIED_FORMAT', '')
num, TAIL = MM.num, MM.TAIL


def glass_substrate_block(layer2, einf, pa, pe):
    """Turn the template's Gen-Osc film block into a 0-oscillator substrate."""
    b = layer2
    b = re.sub(r"^\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+\tF\t'Thickness # 2'[^\n]*",
               "\t\t5000000.0\tF\t-10.0\t100000.0\tF\t'Substrate Thickness'" + TAIL, b, count=1, flags=re.M)
    b = re.sub(r"(start_Gen-Osc Fit Parms\r\n).*?(\r\n\t\tend_Gen-Osc Fit Parms)",
               lambda m: m.group(1) + "\t\t\t0\t%s\tF\t0.0\t10.0\tF\t'Einf'%s\r\n\t\t\t" % (num(einf), TAIL) + m.group(2),
               b, count=1, flags=re.S)
    b = re.sub(r"(start_Gen-Osc Grade Parms\r\n).*?(\r\n\t\tend_Gen-Osc Grade Parms)",
               lambda m: m.group(1) + "\t\t\t" + m.group(2), b, count=1, flags=re.S)
    b = re.sub(r"^\t\t\t[^\t\n]+\t[TF](\t-1000\.0\t1000\.0\tF\t'UV Pole Amp\.')",
               lambda m: '\t\t\t%s\tF%s' % (num(pa), m.group(1)), b, count=1, flags=re.M)
    b = re.sub(r"^\t\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'UV Pole En\.')",
               lambda m: '\t\t\t%s\tF\t6.0\t20.0%s' % (num(pe), m.group(1)), b, count=1, flags=re.M)
    return b.replace('start_Layer2', 'start_Layer0').replace('end_Layer2', 'end_Layer0')


def write_glass(si_text, gg, f, out_path):
    Ls = si_text.split('\r\n')
    i0 = Ls.index('\tstart_Layer0'); i1 = Ls.index('\tend_Layer1'); j0 = Ls.index('\tstart_Layer2'); j1 = Ls.index('\tend_Layer2')
    film = '\r\n'.join(Ls[j0:j1 + 1])
    sub = glass_substrate_block(film, gg['einf'], gg['pole_amp'], gg['pole_en'])
    film = film.replace('start_Layer2', 'start_Layer1').replace('end_Layer2', 'end_Layer1').replace("'Thickness # 2'", "'Thickness # 1'")
    t = '\r\n'.join(Ls[:i0]) + '\r\n' + sub + '\r\n\t\r\n' + film + '\r\n' + '\r\n'.join(Ls[j1 + 1:])
    t = re.sub(r"(start_Model Parms\r\n\t)\d+(\t0\t)", lambda m: m.group(1) + '2' + m.group(2), t, count=1)
    t = re.sub(r"\r\n\tF\t(\r\n\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+\tF\t'# Back Reflections')",
               lambda m: '\r\n\tT\t' + m.group(1), t, count=1)
    t = re.sub(r"^\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'% 1st Reflection')",
               lambda m: '\t%s\tT\t0.0\t100.0%s' % (num(100 * f), m.group(1)), t, count=1, flags=re.M)
    open(out_path, 'w', encoding='utf-8-sig', newline='').write(t)


def parse_glass(path):
    t = open(path, 'rb').read().decode('utf-8-sig')
    sub = re.search(r"start_Layer0(.*?)end_Layer0", t, re.S).group(1)
    film = re.search(r"start_Layer1(.*?)end_Layer1", t, re.S).group(1)
    g_einf = float(re.search(r"start_Gen-Osc Fit Parms\r\n\t+0\t([-\d.eE+]+)", sub).group(1))
    g_pa = float(re.search(r"\n\t\t\t([-\d.eE+]+)\t[TF]\t-1000\.0\t1000\.0\tF\t'UV Pole Amp\.'", sub).group(1))
    g_pe = float(re.search(r"\n\t\t\t([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'UV Pole En\.'", sub).group(1))
    fb = re.search(r"start_Gen-Osc Fit Parms\r\n(.*?)\r\n\t\tend_Gen-Osc Fit Parms", film, re.S).group(1)
    Ls = [x.strip('\t') for x in fb.split('\r\n')]; oscs = []; cur = None
    for ln in Ls[1:]:
        if ln.startswith("'") and ln.endswith("'"): cur = [ln.strip("'"), {}]; oscs.append(cur)
        elif "\t'" in ln and cur is not None:
            fl = ln.split('\t'); cur[1][fl[5].strip("'")] = float(fl[0])
    pa = float(re.search(r"\n\t\t\t([-\d.eE+]+)\t[TF]\t-1000\.0\t1000\.0\tF\t'UV Pole Amp\.'", film).group(1))
    pe = float(re.search(r"\n\t\t\t([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'UV Pole En\.'", film).group(1))
    d = float(re.search(r"\n\t\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'Thickness # 1'", film).group(1)) / 10
    f = float(re.search(r"\n\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'% 1st Reflection'", t).group(1)) / 100
    on = re.search(r"\r\n\t([TF])\t\r\n\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+\tF\t'# Back Reflections'", t).group(1)
    nlay = int(re.search(r"start_Model Parms\r\n\t(\d+)\t", t).group(1))
    p = [pa, pe] + sum([[o[1]['Amp'], o[1]['Br'], o[1]['En']] for o in oscs], [])
    return dict(glass=(g_einf, g_pa, g_pe), p=np.array(p), d=d, f=f, backside_on=on == 'T', n_layers=nlay)


def mse_glass(path, sh):
    m = parse_glass(path)
    wl, P, D, ang = L.se(sh, 300.0, 1000.0, 3); meas = L.ncs(P, D)
    E = 1239.841984 / wl; e, a, b = m['glass']
    ng = np.sqrt(e + a / (b**2 - E**2)) + 0j
    N, C_, S_ = GG.norm(*GL.ncs_glass(wl, GO.N_of(m['p'], wl), m['d'], ng, ang, m['f']))
    r = np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    return M.mse_ncs(r, 1 + 3 / 6), m


def main():
    _os.makedirs(OUTDIR, exist_ok=True)
    CE = json.load(open(W + 'ce_genosc.json')); GJ = json.load(open(W + 'ce_glass.json'))
    tpl = open(C.TEMPLATE, 'rb').read().decode('utf-8-sig')
    print('%-5s %-6s %7s | %-40s | %7s %7s  %s' % ('sheet', 'mat', 'd (nm)', 'file', 'MSE fit', 'MSE file', 'check'))
    for sh, mat in GJ['mat'].items():
        rec = CE[mat]['pick']        # the dispersion tl_ce_glass fitted with (no roughness layer on glass)
        d = GJ['d'][sh]
        si_text = MM.write_si(tpl, mat, rec, d, OUTDIR + '_tmp.mod', 0.0, 0.0)
        fn = '%s_%s_%.1fnm_on_glass_GenOsc_UNVERIFIED.mod' % (sh, mat, d)
        write_glass(si_text, GJ['glass_genosc'], GJ['f'], OUTDIR + fn)
        mse_file, m = mse_glass(OUTDIR + fn, sh)
        ok = m['backside_on'] and m['n_layers'] == 2 and abs(mse_file - GJ['mse'][sh]) < 0.3
        print('%-5s %-6s %7.2f | %-40s | %7.2f %7.2f  %s' % (sh, mat, d, fn[:40], GJ['mse'][sh], mse_file, 'OK' if ok else 'MISMATCH'))
    _os.remove(OUTDIR + '_tmp.mod')
    print('->', OUTDIR)


if __name__ == '__main__':
    main()

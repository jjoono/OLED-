"""CompleteEASE .mod files for the GraHIL and the perovskite-NC pieces, read back and re-scored.

  GraHIL_<d>nm_on_NTVE15A_Si_GenOsc.mod
      Layer0 SI_JAW / Layer1 NTVE_JAW 15 A / Layer2 'GraHIL' (Gen-Osc)
  PerovNC_on_GraHIL_on_NTVE15A_Si_GenOsc.mod
      Layer0 SI_JAW / Layer1 NTVE_JAW 15 A / Layer2 'GraHIL' (Gen-Osc, as above, fit off)
      / Layer3 'Perov_lower' (tabulated) / Layer4 'Perov_dense' (Gen-Osc) + Roughness

The dispersions are the compact Gen-Osc refits of gp_genosc.py.  CompleteEASE has
no layer for "effective medium + scattering loss", so the lower part of the
perovskite is written as a table: Bruggeman(Perov_dense, void, f_b) plus the
extinction k_sc (633 / lambda)^p (held constant below 450 nm), computed from the
same fit.  It does not follow if Perov_dense is refitted in CompleteEASE.

Syntax, all taken from files CompleteEASE itself saved:
  * frame, Gen-Osc layer, tabulated layer (the SI_JAW / NTVE_JAW blocks): the lab's
    ADNC_ITO_GenOSC_v2 model; Gen-Osc rules as in tooling_261117/tl_make_mod.py
  * a named material uses its name in the section tags (start_<name> Fit Parms ...):
    'no5_glass_sub_genosc' and 'Si_JAW' in the lab's tooling snapshot
  * numbers as CompleteEASE prints them (num)
  * Mat Table arrays: base64 of gzip (10-byte zero header, deflate level 6) of
    big-endian float32; wavelength in Angstrom, ascending; e1, e2.  Re-encoding the
    template's own tables with enc() gives back its bytes exactly (checked in main).
  * layers are added by repeating start_LayerN ... end_LayerN and raising the count
    on the first Model Parms line (completeease_mod/make_mod_ag.py)

A Gaussian whose 90 % interval exceeds half its value, or that sits on a bound,
is written with its fit flag off.  Every file is parsed back, the stack rebuilt
from what the FILE says (tables decoded, oscillators re-evaluated), and its MSE
against the data must equal the fit's."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, re, json, base64, gzip, zlib, struct, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gp_load as G, gp_model as P, gp_grahil as GH, gp_perov as PV, gp_genosc as GO

W = _os.path.join(ELLIPS_OUT, ''); OUTDIR = _os.path.join(ELLIPS_OUT, 'grahil_perov_mod', '')
TEMPLATE = P.C.TEMPLATE
NG_GRAHIL, NG_PEROV = GO.NG_GRAHIL, 7       # which refits of gp_genosc.py to write
TAIL = "\tF\tF\t0.0\t0.0\t100000.0\tF\t100.0\t"
def num(x):
    """a float the way CompleteEASE writes one (1.0E-8, not 1e-08)"""
    r = repr(float(x))
    if 'e' not in r: return r
    m, e = r.split('e'); return (m if '.' in m else m + '.0') + 'E' + str(int(e))
pl = lambda v, f, lo, hi, nm: "\t\t\t%s\t%s\t%s\t%s\tF\t'%s'%s" % (v, 'T' if f else 'F', lo, hi, nm, TAIL)
gl = lambda lo, hi, i: "\t\t\t0.0\tF\t%s\t%s\tF\t'Grade %d'%s" % (lo, hi, i, TAIL)
E_TAB = np.arange(0.70, 6.6001, 0.005)[::-1]        # table grid: uniform in eV like the JAW tables, wavelength ascending


# ------------------------------------------------------------------ encoding
def enc(a):
    d = np.asarray(a, '>f4').tobytes(); co = zlib.compressobj(6, zlib.DEFLATED, -15); z = co.compress(d) + co.flush()
    return base64.b64encode(bytes.fromhex('1f8b0800000000000000') + z + struct.pack('<II', zlib.crc32(d) & 0xffffffff, len(d) & 0xffffffff)).decode()


def dec(b64):
    return np.frombuffer(gzip.decompress(base64.b64decode(b64)), '>f4').astype(float)


# ------------------------------------------------------------------ blocks
def free_flag(p):
    return (p['ci'] <= 0.5 * abs(p['v'])) and not p['at_bound']


def genosc_layer(L2, idx, name, d_nm, fit_d, einf, pole_amp, pole_en, gauss, lock_all=False, ir_amp=None):
    """einf, pole_amp, pole_en and the entries Amp / Br / En of gauss are parameter records of gp_genosc
    (v, ci, lo, hi, at_bound); pole_en may be a plain number (fixed)."""
    on = lambda p: 'T' if (free_flag(p) and not lock_all) else 'F'
    pe = pole_en if isinstance(pole_en, dict) else dict(v=pole_en, ci=np.inf, at_bound=True)
    fit = ["\t\t\t%d\t%s\t%s\t0.0\t10.0\tF\t'Einf'%s" % (len(gauss), num(einf['v']), on(einf), TAIL)]
    grade = []
    for g in gauss:
        fit.append("\t\t\t'Gaussian'\t")
        for key in ('Amp', 'iAmp', 'Br', 'En'):
            if key == 'iAmp':
                fit.append(pl('0.0', False, '-1000.0', '1000.0', 'iAmp'))
            else:
                p = g[key]; fit.append(pl(num(p['v']), on(p) == 'T', num(p['lo']), num(p['hi']), key))
        for glo, ghi in (('-100.0', '1000.0'), ('-1000.0', '1000.0'), ('0.0', '100.0'), ('1.0E-8', '15.0')):
            grade.append('\t\t\tF\t'); grade += [gl(glo, ghi, j) for j in (1, 2, 3)]
    fit.append('\t\t\t'); grade.append('\t\t\t')
    b = L2
    b = re.sub(r"(start_Gen-Osc Fit Parms\r\n).*?(\r\n\t\tend_Gen-Osc Fit Parms)", lambda m: m.group(1) + '\r\n'.join(fit) + m.group(2), b, count=1, flags=re.S)
    b = re.sub(r"(start_Gen-Osc Grade Parms\r\n).*?(\r\n\t\tend_Gen-Osc Grade Parms)", lambda m: m.group(1) + '\r\n'.join(grade) + m.group(2), b, count=1, flags=re.S)
    b = re.sub(r"^\t\t\t[^\t\n]+\t[TF](\t-1000\.0\t1000\.0\tF\t'UV Pole Amp\.')",
               lambda m: '\t\t\t%s\t%s%s' % (num(pole_amp['v']), on(pole_amp), m.group(1)), b, count=1, flags=re.M)
    b = re.sub(r"^\t\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'UV Pole En\.')",
               lambda m: '\t\t\t%s\t%s\t%s\t%s%s' % (num(pe['v']), on(pe), num(pe.get('lo', 1.0e-8)), num(pe.get('hi', 15.0)), m.group(1)), b, count=1, flags=re.M)
    b = re.sub(r"^\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+\tF\t'Thickness # 2'",
               lambda m: "\t\t%s\t%s\t%s\t%s\tF\t'Thickness # %d'" % (num(d_nm * 10), 'T' if fit_d else 'F', num(0.8 * d_nm * 10), num(1.2 * d_nm * 10), idx), b, count=1, flags=re.M)
    if ir_amp is not None:
        b, k = re.subn(r"^\t\t\t[^\t\n]+\t[TF](\t-1000\.0\t1000\.0\tF\t'IR Pole Amp\.')", lambda m: '\t\t\t%s\t%s%s' % (num(ir_amp['v']), on(ir_amp), m.group(1)), b, count=1, flags=re.M)
        assert k == 1
    b = b.replace("'Layer'\t'Gen-Osc'\t'GENOSC'\t''\t'Gen-Osc'", "'Layer'\t'%s'\t'GENOSC'\t''\t'%s'" % (name, name), 1)
    for sec in ('Misc Parms', 'Fit Parms', 'Grade Parms', 'Permanent Poles'):
        b = b.replace('start_Gen-Osc ' + sec, 'start_%s %s' % (name, sec)).replace('end_Gen-Osc ' + sec, 'end_%s %s' % (name, sec))
    return b.replace('start_Layer2', 'start_Layer%d' % idx).replace('end_Layer2', 'end_Layer%d' % idx)


def table_layer(L1, idx, name, comment, d_nm, fit_d, eps):
    """a tabulated film on E_TAB, built on the template's NTVE_JAW block"""
    assert "'" not in comment
    wl_A = P.HC / E_TAB * 10.0
    b = L1
    b = re.sub(r"^\t\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+\tT\t'Native Oxide'\tF\tT",
               lambda m: "\t\t%s\t%s\t%s\t%s\tF\t'Thickness # %d'\tF\tF" % (num(d_nm * 10), 'T' if fit_d else 'F', num(0.8 * d_nm * 10), num(1.2 * d_nm * 10), idx), b, count=1, flags=re.M)
    b = b.replace("'Layer'\t'NTVE_JAW'\t'EV'\t''\t'NTVE_JAW'", "'Layer'\t'%s'\t'EV'\t''\t'%s'" % (name, name), 1)
    b = re.sub(r"(start_File Info\r\n\t\t\t)'[^']*'", lambda m: m.group(1) + "'%s'" % comment, b, count=1)
    b = re.sub(r"(start_Mat Table\r\n\t\t\t)\d+", lambda m: m.group(1) + str(len(E_TAB)), b, count=1)
    for tag, arr in (('Wvl', wl_A), ('e1', eps.real), ('e2', eps.imag)):
        b = re.sub(r"(start_%s Array\r\n\t\t\t\t)'[^']+'" % tag, lambda m: m.group(1) + "'%s'" % enc(arr), b, count=1)
    b = b.replace('start_NTVE_JAW Fit Parms', 'start_%s Fit Parms' % name).replace('end_NTVE_JAW Fit Parms', 'end_%s Fit Parms' % name)
    return b.replace('start_Layer1', 'start_Layer%d' % idx).replace('end_Layer1', 'end_Layer%d' % idx)


def assemble(tpl, layers, d_ox_nm, rough_nm, fit_rough):
    """layers: finished blocks for Layer2.. (bottom to top), placed over the template's substrate and oxide"""
    Ls = tpl.split('\r\n'); i0 = Ls.index('\tstart_Layer2'); i1 = Ls.index('\tend_Layer2')
    t = '\r\n'.join(Ls[:i0]) + '\r\n' + '\r\n\t\r\n'.join(layers) + '\r\n' + '\r\n'.join(Ls[i1 + 1:])
    t = re.sub(r"(start_Model Parms\r\n\t)\d+(\t0\t)", lambda m: m.group(1) + str(2 + len(layers)) + m.group(2), t, count=1)
    t = re.sub(r"^\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Angle Offset')", lambda m: '\t0.0\tF\t-2.0\t2.0' + m.group(1), t, count=1, flags=re.M)
    t = re.sub(r"^\t[^\t\n]+\t[TF]\t[^\t\n]+\t[^\t\n]+(\tF\t'Roughness')",
               lambda m: '\t%s\t%s\t0.0\t300.0%s' % (num(rough_nm * 10), 'T' if fit_rough else 'F', m.group(1)), t, count=1, flags=re.M)
    t = re.sub(r"^\t\t[^\t\n]+\t[TF](\t[^\t\n]+\t[^\t\n]+\tT\t'Native Oxide')", lambda m: '\t\t%s\tF%s' % (num(d_ox_nm * 10), m.group(1)), t, count=1, flags=re.M)
    return t


def split_template(tpl):
    Ls = tpl.split('\r\n')
    blk = lambda k: '\r\n'.join(Ls[Ls.index('\tstart_Layer%d' % k): Ls.index('\tend_Layer%d' % k) + 1])
    return blk(1), blk(2)


# ------------------------------------------------------------------ read back
def parse(path):
    """-> dict(n_layers, rough, angle_offset, layers=[dict(idx, name, type, d, eps(wl_nm))] bottom to top)"""
    t = open(path, 'rb').read().decode('utf-8-sig'); assert t.count('\r\n') == t.count('\n')
    out = dict(n_layers=int(re.search(r"start_Model Parms\r\n\t(\d+)\t0\t", t).group(1)),
               rough=float(re.search(r"\n\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'Roughness'", t).group(1)) / 10,
               angle_offset=float(re.search(r"\n\t([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'Angle Offset'", t).group(1)), layers=[])
    for m in re.finditer(r"\n\tstart_Layer(\d+)\r\n(.*?)\r\n\tend_Layer\1\r\n", t, re.S):
        body = m.group(2); Ls = body.split('\r\n'); f = Ls[1].strip('\t').split('\t')
        assert f[0] == "'Layer'" and f[1] == f[4]
        lay = dict(idx=int(m.group(1)), name=f[1].strip("'"), type=f[2].strip("'"), d=float(Ls[0].strip('\t').split('\t')[0]) / 10,
                   d_name=Ls[0].strip('\t').split('\t')[5].strip("'"), d_fit=Ls[0].strip('\t').split('\t')[1] == 'T')
        nm = re.escape(lay['name'])
        if lay['type'] == 'EV':
            npt = int(re.search(r"start_Mat Table\r\n\t+(\d+)\t", body).group(1))
            arr = {tag: dec(re.search(r"start_%s Array\r\n\t+'([^']+)'" % tag, body).group(1)) for tag in ('Wvl', 'e1', 'e2')}
            assert all(len(a) == npt for a in arr.values()) and np.all(np.diff(arr['Wvl']) > 0)
            assert ('start_%s Fit Parms' % lay['name']) in body
            lay['eps'] = (lambda a: (lambda wl: np.interp(wl, a['Wvl'] / 10, a['e1']) + 1j * np.interp(wl, a['Wvl'] / 10, a['e2'])))(arr)
            lay['table'] = (float(arr['Wvl'][0] / 10), float(arr['Wvl'][-1] / 10), npt)
        elif lay['type'] == 'GENOSC':
            fb = re.search(r"start_%s Fit Parms\r\n(.*?)\r\n\t\tend_%s Fit Parms" % (nm, nm), body, re.S).group(1)
            rows = [x.strip('\t') for x in fb.split('\r\n')]; n_osc = int(rows[0].split('\t')[0]); einf = float(rows[0].split('\t')[1])
            oscs = []; flags = []
            for ln in rows[1:]:
                if ln == "'Gaussian'": oscs.append({})
                elif "\t'" in ln:
                    q = ln.split('\t'); oscs[-1][q[5].strip("'")] = float(q[0]); flags.append(q[1] == 'T')
            gb = re.search(r"start_%s Grade Parms\r\n(.*?)\r\n\t\tend_%s Grade Parms" % (nm, nm), body, re.S).group(1)
            assert n_osc == len(oscs) and all(set(o) == {'Amp', 'iAmp', 'Br', 'En'} and o['iAmp'] == 0 for o in oscs)
            assert sum(1 for x in gb.split('\r\n') if x.strip() == 'F') == 4 * n_osc
            assert ('start_%s Permanent Poles' % lay['name']) in body and ('start_%s Misc Parms' % lay['name']) in body
            pa = float(re.search(r"\n\t\t\t([-\d.eE+]+)\t[TF]\t-1000\.0\t1000\.0\tF\t'UV Pole Amp\.'", body).group(1))
            pe = float(re.search(r"\n\t\t\t([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'UV Pole En\.'", body).group(1))
            gs = [(o['Amp'], o['Br'], o['En']) for o in oscs]
            ir = float(re.search(r"\n\t\t\t([-\d.eE+]+)\t[TF]\t-1000\.0\t1000\.0\tF\t'IR Pole Amp\.'", body).group(1))
            lay['eps'] = (lambda e, a, b, g, i: (lambda wl: GO.eps_genosc(P.HC / wl, e, a, b, g, i)))(einf, pa, pe, gs, ir)
            lay.update(einf=einf, pole=(pa, pe), ir=ir, n_osc=n_osc, n_free=int(sum(flags)))
        else:
            raise ValueError(lay['type'])
        out['layers'].append(lay)
    assert out['n_layers'] == len(out['layers']) and [l['idx'] for l in out['layers']] == list(range(out['n_layers']))
    assert [l['name'] for l in out['layers'][:2]] == ['SI_JAW', 'NTVE_JAW']
    return out


def mse_from_file(path, sheet, lo, hi, step, n_parm):
    m = parse(path); d = G.load(sheet, lo, hi, step); wl = d['wl']
    sq = lambda e: (lambda n: np.where(n.imag < 0, -n, n))(np.sqrt(e.astype(complex)))
    film = [sq(l['eps'](wl)) for l in m['layers'][2:]][::-1]                    # top first
    L = [P.n_air(wl)] + ([P.bruggeman(film[0], P.n_air(wl), 0.5)] if m['rough'] > 0 else []) + film + [P.n_ox(wl), P.n_si(wl)]
    D = ([m['rough']] if m['rough'] > 0 else []) + [l['d'] for l in m['layers'][2:]][::-1] + [m['layers'][1]['d']]
    r = P.resid(P.ncs(wl, d['ang'], L, D), P.data_ncs(d))
    return float(P.mse_ncs(r, n_parm)), m


# ------------------------------------------------------------------ main
def gauss_of(parms, ng):
    return [{k: parms['%s%d' % (k, i + 1)] for k in ('Amp', 'Br', 'En')} for i in range(ng)]


def main():
    _os.makedirs(OUTDIR, exist_ok=True)
    tpl = open(TEMPLATE, 'rb').read().decode('utf-8-sig'); L1, L2 = split_template(tpl)
    for tag in ('Wvl', 'e1', 'e2'):                                             # the encoder reproduces CompleteEASE's own bytes
        b64 = re.search(r"start_%s Array\r\n\t+'([^']+)'" % tag, L1).group(1); assert enc(dec(b64)) == b64
    J = json.load(open(W + 'genosc.json')); g = J['grahil'][str(NG_GRAHIL)]; p = J['perov'][str(NG_PEROV)]
    gp, pp = g['parms'], p['parms']; rows = []

    # ---- GraHIL on Si
    def grahil_block(lock):
        return genosc_layer(L2, 2, 'GraHIL', gp['d']['v'], not lock, gp['einf'], gp['pole_amp'], g['pole_en'], gauss_of(gp, NG_GRAHIL), lock, gp['ir_amp'])
    f1 = 'GraHIL_%.1fnm_on_NTVE%.0fA_Si_GenOsc.mod' % (gp['d']['v'], PV.D_OX * 10)
    open(OUTDIR + f1, 'w', encoding='utf-8-sig', newline='').write(assemble(tpl, [grahil_block(False)], PV.D_OX, 0.0, False))
    mse, m = mse_from_file(OUTDIR + f1, 'layer1', GH.LO, GH.HI, GH.STEP, len(gp))
    rows.append((f1, 'layer1', '%g-%g' % (GH.LO, GH.HI), g['mse'], mse, m))

    # ---- perovskite NC film on GraHIL on Si
    v = {k: pp[k]['v'] for k in PV.STRUCT}; ng = NG_PEROV
    gs = [tuple(pp['%s%d' % (k, i + 1)]['v'] for k in ('Amp', 'Br', 'En')) for i in range(ng)]
    wl_t = P.HC / E_TAB; Nd = GO.n_of(E_TAB, pp['einf']['v'], pp['pole_amp']['v'], pp['pole_en']['v'], gs)
    Nlow = P.bruggeman(Nd, np.ones(len(wl_t), complex), v['f_b']) + 1j * v['k_sc'] * (633.0 / np.maximum(wl_t, p['sc_hold']))**v['p']
    com = 'Perovskite NC film, lower part: effective medium of Perov_dense and void with scattering loss, fitted 205-1000 nm (table extrapolated outside)'
    lay = [grahil_block(True),
           table_layer(L1, 3, 'Perov_lower', com, v['d_b'], True, Nlow**2),
           genosc_layer(L2, 4, 'Perov_dense', v['d_t'], True, pp['einf'], pp['pole_amp'], pp['pole_en'], gauss_of(pp, ng))]
    lay[0] = re.sub(r"^\t\t[^\t\n]+\t[TF](\t[^\t\n]+\t[^\t\n]+\tF\t'Thickness # 2')", lambda q: '\t\t%s\tF%s' % (num(p['d_g'] * 10), q.group(1)), lay[0], count=1, flags=re.M)
    f2 = 'PerovNC_on_GraHIL_on_NTVE%.0fA_Si_GenOsc.mod' % (PV.D_OX * 10)
    open(OUTDIR + f2, 'w', encoding='utf-8-sig', newline='').write(assemble(tpl, lay, PV.D_OX, v['d_r'], v['d_r'] > 0.3))
    mse, m = mse_from_file(OUTDIR + f2, 'layer2', PV.LO, PV.HI, PV.STEP, len(pp))
    rows.append((f2, 'layer2', '%g-%g' % (PV.LO, PV.HI), p['mse'], mse, m))

    summ = []
    print('%-48s %-7s %-9s %8s %8s  %s' % ('file', 'sheet', 'nm', 'MSE fit', 'MSE file', 'check'))
    for fn, sh, rng, mfit, mfile, m in rows:
        ok = abs(mfile - mfit) < 0.05
        print('%-48s %-7s %-9s %8.2f %8.2f  %s' % (fn, sh, rng, mfit, mfile, 'OK' if ok else 'MISMATCH'))
        for l in m['layers']:
            extra = ('Einf %.3f, UV pole %.2f @ %.2f eV, IR pole %.4f, %d Gaussians (%d of their parameters free)' % (l['einf'], *l['pole'], l['ir'], l['n_osc'], l['n_free'])) if l['type'] == 'GENOSC' else \
                    ('table %.0f-%.0f nm, %d points' % l['table'])
            print('     Layer%d %-12s %-7s %9.2f nm (%s, fit %s)  %s' % (l['idx'], l['name'], l['type'], l['d'], l['d_name'], 'on' if l['d_fit'] else 'off', extra))
        print('     Roughness %.2f nm' % m['rough'])
        summ.append(dict(file=fn, sheet=sh, range_nm=rng, mse_fit=mfit, mse_file=mfile, ok=bool(ok), roughness_nm=m['rough'],
                         layers=[dict(idx=l['idx'], name=l['name'], type=l['type'], d_nm=l['d']) for l in m['layers']]))
    json.dump(summ, open(OUTDIR + 'mod_check.json', 'w'), indent=1)
    print('->', OUTDIR)


if __name__ == '__main__':
    main()

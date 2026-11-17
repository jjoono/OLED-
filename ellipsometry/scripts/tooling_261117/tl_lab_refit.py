"""The lab's evaporator tooling models with this set's data put in, refitted
from scratch.  Models are read straight from the snapshots in
Tooling_Evaporator.zip (EVAP_ZIP, default $ELLIPS_DATA/Tooling_Evaporator.zip):

  NPB_tooling_S8_S13  uniaxial Gen-Osc (o: Einf + 6 Gaussians, e: Einf + 5),
                      the lab fits the whole dispersion; NTVE 14.5 A on Si
  TAPC_S11_tooling    isotropic Gen-Osc, fixed; thickness and the glass
                      backside fraction fitted; NTVE 14.5 A on Si
  BPhen_S12_tooling   isotropic Gen-Osc, fixed; thickness fitted; NTVE 0 on Si

Data sets are swapped one for one (same substrate, closest nominal thickness):
  NPB    #1 S8 30 nm Si -> 1_21, #2 S13 10 nm Si -> 1_23, #3 S13 50 nm Si -> 1_22,
         #4 S13 10 nm glass -> 1_33, #5 S8 30 nm glass -> 1_31, #6 S13 50 nm glass -> 1_32
  TAPC   #1 Si -> 1_24, #2 glass -> 1_34 (#3, a second glass, has no counterpart)
  Bphen  2_21 2_22 2_23 2_24 on Si, 2_31 2_32 2_33 2_34 on glass (the model is
         thickness-only, so the pieces do not influence one another)
Everything else is the lab's: layer stack (EMA substrate with Si_JAW or the
model's own glass at EMA % 0 / 100, NTVE_JAW, the film), per-data-set values,
fit flags, parameter bounds, 245.8-1688 nm.

The NPB snapshot's glass absorbs (k 0.0025-0.005: a 1 mm substrate would be
opaque) and its back reflections are off; the lab left its glass sets out of
the fit.  So NPB runs twice:
  NPB as saved     glass sets with the NPB model's glass, no backside
  NPB + TAPC glass glass and backside handling of the lab's TAPC model
                   (transparent no5 glass, backside on, % 1st Reflection fitted)
each with the dispersion fixed (old n,k) and refitted with the lab's flags.
Backside beams: Re(n_glass) at the interfaces, Beer-Lambert through the 1 mm
substrate, film axis-averaged on the back path; 'normalised' N,C,S.

Thickness prior: the deposition is expected within 5 % of nominal for 30 and
50 nm and within 10 % for the thin 10 nm films, so every thickness is bounded
to nominal x (1 +- TOL[nominal]).  The thickness-
only fits are also run unbounded, and every bounded result that ends on a
bound is flagged: there the number is the constraint, not the data."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, re, json, zipfile, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_ce_genosc as C, tl_glass_global as GG, tl_bphen10 as B
import jnk_model as M, ce_osc as osc
from tmm_fix import tmm
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
EVAP_ZIP = _os.environ.get('EVAP_ZIP', _os.path.join(ELLIPS_DATA, 'Tooling_Evaporator.zip'))
HC = 1239.841984
WL = (245.0, 1690.0, 2)
T_SUB = 1.0e6                                   # nm: 'Substrate Thickness' 1.0E7 A in every model
TOL = {10: 0.10, 30: 0.05, 50: 0.05}             # relative, by nominal thickness
NOMINAL = {'1_21': 30, '1_31': 30, '1_22': 50, '1_32': 50, '1_23': 10, '1_33': 10, '1_24': 50, '1_34': 50,
           '2_21': 50, '2_22': 50, '2_31': 50, '2_32': 50, '2_23': 10, '2_24': 10, '2_33': 10, '2_34': 10}
def dbounds(sh, bounded):
    n = NOMINAL[sh]
    return (n * (1 - TOL[n]), n * (1 + TOL[n])) if bounded else (0.0, 300.0)
PL = re.compile(r"^\t*([-\d.eE+]+|NaN)\t([TF])\t([-\d.eE+]+)\t([-\d.eE+]+)\t[TF]\t'([^']*)'", re.M)


def fitlog(name):
    with zipfile.ZipFile(EVAP_ZIP) as z:
        with zipfile.ZipFile(z.open(name + '.SEsnap')) as s:
            return s.read('_FitLog').decode('latin-1').replace('\r', '')


def params(seg):
    return [dict(v=float(m.group(1)), fit=m.group(2) == 'T', lo=float(m.group(3)), hi=float(m.group(4)), name=m.group(5))
            for m in PL.finditer(seg)]


def genosc(text, start):
    """Gen-Osc block whose 'Fit Parms' starts at `start`: Einf, oscillators, three poles."""
    j = text.index('Permanent Poles', start); e = text.index('end_', start)     # Fit Parms has no nested blocks
    m = re.search(r"^\t*\d+\t([-\d.eE+]+)\t([TF])\t([-\d.eE+]+)\t([-\d.eE+]+)\t[TF]\t'Einf'", text[start:e], re.M)
    einf = dict(v=float(m.group(1)), fit=m.group(2) == 'T', lo=float(m.group(3)), hi=float(m.group(4)), name='Einf')
    body = params(text[start:e]); poles = params(text[j:j + 2000])[:3]
    assert [p['name'] for p in poles] == ['UV Pole Amp.', 'UV Pole En.', 'IR Pole Amp.']
    types = re.findall(r"^\t*'(Gaussian|Tauc-Lorentz)'\t*$", text[start:e], re.M)
    g = dict(einf=einf, osc=[], uv_amp=poles[0], uv_en=poles[1], ir_amp=poles[2])
    k = 0
    for t in types:
        n = 4 if t == 'Gaussian' else 4                    # Gaussian: Amp iAmp Br En; TL: Amp Br Eo Eg
        g['osc'].append((t, {p['name']: p for p in body[k:k + n]})); k += n
    assert k == len(body), 'unparsed Gen-Osc parameters'
    return g


def parse(name):
    t = fitlog(name)
    ms = t[t.index('start_Model Structure'):t.index('end_Model Structure')]
    L0 = ms[ms.index('start_Layer0'):ms.index('end_Layer0')]
    L1 = ms[ms.index('start_Layer1'):ms.index('end_Layer1')]
    L2 = ms[ms.index('start_Layer2'):]
    gname = re.search(r"'([^']+)'\t'GENOSC'", L0).group(1)
    glass = genosc(L0, L0.index('start_%s Fit Parms' % gname))
    starts = [m.start() for m in re.finditer(r'start_Gen-Osc Fit Parms', L2)]
    if "'Biaxial'\t'BIAXIAL'" in L2:
        mode = re.search(r"start_Biaxial Fit Parms\n\t*'([^']+)'", L2).group(1)
        assert mode == 'Uniaxial' and len(starts) == 3
        film = dict(kind='uniaxial', o=genosc(L2, starts[0]), e=genosc(L2, starts[2]))
    else:
        film = dict(kind='iso', o=genosc(L2, starts[0]))
    mp = t[t.index('start_MultiSamp'):t.index('end_MultiSamp')]
    npar, nset = map(int, re.search(r'start_MultiSamp\n\t*(\d+)\t(\d+)', mp).groups())
    msa = {}
    for m in re.finditer(r"^\t*'([^']+)'\t*\n((?:.*\n){%d})" % nset, mp, re.M):
        msa[m.group(1)] = params(m.group(2))
    parms = t[t.index('start_Model Parms'):t.index('start_MultiSamp')]
    back = re.search(r"\n\t*([TF])\t*\n[^\n]*'# Back Reflections'", parms).group(1) == 'T'
    sel = re.search(r'start_Dataset Selection\n\t*\d+\t*\n((?:\t*[TF]\t*\n)+)', t).group(1).split()
    return dict(name=name, glass=glass, glass_name=gname, film=film, msa=msa, nset=nset, back=back, selected=sel,
                ntve=params(L1)[0])


# ------------------------------------------------------------------ optics
def gauss(E, A, Br, En):
    if Br < 1e-9:
        return 0.0, 0.0
    return osc.gaussian(E, A, Br, En)


def eps(g, E, vals=None):
    """vals: optional {id(param dict): value} overriding the stored values."""
    v = (lambda p: vals.get(id(p), p['v'])) if vals else (lambda p: p['v'])
    e = v(g['einf']) + v(g['uv_amp']) / (v(g['uv_en'])**2 - E**2) - v(g['ir_amp']) / E**2 + 0j
    for t, q in g['osc']:
        if t == 'Gaussian':
            a, b = gauss(E, v(q['Amp']), v(q['Br']), v(q['En']))
        else:
            a, b = osc.tauc_lorentz(E, v(q['Amp']), v(q['Br']), v(q['Eo']), v(q['Eg']))
        e = e + a + 1j * b
    return e


def ncs_si(wl, ang, eo, ee, d, dox):
    one = np.ones(len(wl), complex); ox = C.N_ce('NTVE_JAW', wl)**2; si = C.N_ce('SI_JAW', wl)**2
    lo, le, ds = [one, eo], [one, ee], [d]
    if dox > 0:
        lo.append(ox); le.append(ox); ds.append(dox)
    return B.ncs_regions(wl, ang, [(1.0, lo + [si], le + [si], ds)], GG.CONV['normalised'])


def ncs_glass(wl, ang, eo, ee, d, eg, f, back, nb=20):
    one = np.ones(len(wl), complex); ngc = B.sqrtc(eg); ng = ngc.real + 0j; kg = ngc.imag
    out = []
    for th in ang:
        rp = B.tmm_uni(wl, [one, eo, eg], [one, ee, eg], [d], th)[0]
        rs = B.tmm_uni(wl, [one, eo, eg], [one, eo, eg], [d], th)[1]
        Ip, Is, X = np.abs(rp)**2, np.abs(rs)**2, np.conj(rp) * rs
        if back and f > 0:
            Ni = B.sqrtc((2 * eo + ee) / 3)
            thg = np.rad2deg(np.arcsin(np.sin(np.deg2rad(th)) / ng.real))
            att = np.exp(-4 * np.pi * kg * T_SUB / (wl * np.cos(np.deg2rad(thg))))
            tdp, tds = M.t_amp(wl, [one, Ni, ng], [d], th)
            tup, tus = M.t_amp(wl, [ng, Ni, one], [d], thg)
            rip, ris = tmm(wl, [ng, Ni, one], [d], thg)
            rbp, rbs = tmm(wl, [ng, one], [], thg); rbp, rbs = rbp * att, rbs * att
            ap = tdp * rbp * tup; as_ = tds * rbs * tus
            for _ in range(nb):
                Ip = Ip + f * np.abs(ap)**2; Is = Is + f * np.abs(as_)**2; X = X + f * np.conj(ap) * as_
                ap = ap * rip * rbp; as_ = as_ * ris * rbs
        I = Ip + Is
        out.append(GG.CONV['normalised']((Is - Ip) / I, 2 * X.real / I, 2 * X.imag / I))
    return tuple(np.stack([o[k] for o in out], 1) for k in range(3))


# ------------------------------------------------------------------ runs
RUNS = {
    'NPB as saved':          dict(snap='NPB_tooling_S8_S13', sets=['1_21', '1_23', '1_22', '1_33', '1_31', '1_32'], glass_from=None),
    'NPB + TAPC glass':      dict(snap='NPB_tooling_S8_S13', sets=['1_21', '1_23', '1_22', '1_33', '1_31', '1_32'], glass_from='TAPC_S11_tooling'),
    'TAPC':                  dict(snap='TAPC_S11_tooling', sets=['1_24', '1_34'], template_sets=[0, 1], glass_from=None),
    'Bphen':                 dict(snap='BPhen_S12_tooling', sets=['2_21', '2_22', '2_23', '2_24', '2_31', '2_32', '2_33', '2_34'], glass_from=None),
    # the Bphen model fits Si with NTVE 0; the NPB and TAPC models use 14.5 A - same footing:
    'Bphen, NTVE 14.5 A':    dict(snap='BPhen_S12_tooling', sets=['2_21', '2_22', '2_23', '2_24', '2_31', '2_32', '2_33', '2_34'], glass_from=None),
}
is_glass = lambda sh: sh.split('_')[1].startswith('3')


def setup(run):
    R = RUNS[run]; mdl = parse(R['snap'])
    glass, back = mdl['glass'], mdl['back']
    if R['glass_from']:
        glass = parse(R['glass_from'])['glass']; back = True
    if run.startswith('Bphen'):
        back = True                         # glass sets: the TAPC model's handling (the S12 sets were all Si)
    msa = mdl['msa']; tidx = R.get('template_sets', list(range(len(R['sets']))))
    sets = []
    for j, sh in enumerate(R['sets']):
        g = is_glass(sh)
        if run.startswith('Bphen'):
            dox = 0.0 if (g or run == 'Bphen') else 1.45
        else:
            dox = msa['Thickness # 1'][tidx[j]]['v'] / 10 if 'Thickness # 1' in msa else mdl['ntve']['v'] / 10
            if g: dox = 0.0
        ema = 100.0 if g else 0.0
        if not run.startswith('Bphen') and 'EMA % (Mat 2)' in msa:
            assert msa['EMA % (Mat 2)'][tidx[j]]['v'] == ema, (run, sh)
        wl, P, D, ang = L.se(sh, *WL); E = HC / wl
        sets.append(dict(sh=sh, glass=g, dox=dox, wl=wl, E=E, ang=ang, meas=L.ncs(P, D),
                         eg=eps(glass, E) if g else None, back=back and g))
    return mdl, sets


def model_ncs(film, s, d, f, vals=None):
    eo = eps(film['o'], s['E'], vals); ee = eps(film['e'], s['E'], vals) if film['kind'] == 'uniaxial' else eo
    if s['glass']:
        return ncs_glass(s['wl'], s['ang'], eo, ee, d, s['eg'], f, s['back'])
    return ncs_si(s['wl'], s['ang'], eo, ee, d, s['dox'])


def resid(m, s):
    return np.concatenate([(m[k] - s['meas'][k]).ravel() for k in range(3)])


def fit_thickness(film, s, bounded=False):
    """One set, dispersion fixed: d (and f on glass with backside), grid start in d."""
    nf = 1 if s['back'] else 0; best = None; dlo, dhi = dbounds(s['sh'], bounded)
    grid = np.linspace(dlo + 1e-3, dhi - 1e-3, 5) if bounded else np.arange(1.0, 121.0, 4.0)
    for d0 in grid:
        x0 = [d0] + [0.7] * nf
        r = least_squares(lambda x: resid(model_ncs(film, s, x[0], x[1] if nf else 0.0), s), x0,
                          bounds=([dlo] + [0.0] * nf, [dhi] + [1.0] * nf))
        if best is None or r.cost < best.cost: best = r
    return best


def flagged(film):
    out = []
    for key in ('o', 'e') if film['kind'] == 'uniaxial' else ('o',):
        g = film[key]
        for p in [g['einf'], g['uv_amp'], g['uv_en'], g['ir_amp']] + [q for _, o in g['osc'] for q in o.values()]:
            if p['fit'] and p['name'] != 'iAmp':        # iAmp is flagged but never in the lab's fit results
                out.append((key, p))
    return out


def on_bound(sh, d, bounded):
    lo, hi = dbounds(sh, bounded)
    return bool(bounded and min(d - lo, hi - d) < 1e-3)


def run(name, refit, bounded=True, physical=False):
    """physical: the refit keeps the lab's structure and flags but bounds Gaussian
    amplitudes >= 0, Einf >= 1 and Br >= 0.05 eV (the lab bounds allow amplitudes
    to -100, and a free refit used them: negative eps2 = gain)."""
    mdl, sets = setup(name); film = mdl['film']
    res = [fit_thickness(film, s, bounded) for s in sets]
    d = [float(r.x[0]) for r in res]; f = [float(r.x[1]) if s['back'] else None for r, s in zip(res, sets)]
    out = dict(run=name, refit=False, bounded=bounded, sets=[s['sh'] for s in sets], d=d, f=f,
               d_on_bound=[on_bound(s['sh'], v, bounded) for s, v in zip(sets, d)],
               mse={s['sh']: float(M.mse_ncs(r.fun, len(r.x))) for s, r in zip(sets, res)})
    if not refit:
        allr = np.concatenate([r.fun for r in res]); out['mse_all'] = float(M.mse_ncs(allr, sum(len(r.x) for r in res)))
        return out
    fp = flagged(film); npf = len(fp)
    gi = [j for j, s in enumerate(sets) if s['back']]
    x0 = np.r_[[p['v'] for _, p in fp], d, [f[j] for j in gi]]
    lo = np.r_[[p['lo'] for _, p in fp], [dbounds(s['sh'], bounded)[0] for s in sets], [0.0] * len(gi)]
    hi = np.r_[[p['hi'] for _, p in fp], [dbounds(s['sh'], bounded)[1] for s in sets], [1.0] * len(gi)]
    if physical:
        floor = {'Amp': 0.0, 'Einf': 1.0, 'Br': 0.05}
        lo[:npf] = [max(l, floor.get(p['name'], l)) for l, (_, p) in zip(lo[:npf], fp)]
    span = np.maximum(hi - lo, 1e-6); x0 = np.clip(x0, lo + 1e-6 * span, hi - 1e-6 * span)
    def unpack(x):
        return {id(p): v for (_, p), v in zip(fp, x[:npf])}
    def per(x, j):
        s = sets[j]; ff = x[npf + len(sets) + gi.index(j)] if j in gi else 0.0
        return resid(model_ncs(film, s, x[npf + j], ff, unpack(x)), s)
    allres = lambda x: np.concatenate([per(x, j) for j in range(len(sets))])
    r = least_squares(allres, x0, bounds=(lo, hi), x_scale='jac', max_nfev=400)
    x = r.x; vals = unpack(x)
    J = r.jac; s2 = np.sum(r.fun**2) / max(r.fun.size - x.size, 1)
    ci = 1.645 * np.sqrt(np.clip(np.diag(np.linalg.pinv(J.T @ J) * s2), 0, None))
    wl = np.array([350.0, 450.0, 550.0, 633.0]); E = HC / wl
    nk = {k: [[float(v.real), float(v.imag)] for v in B.sqrtc(eps(film[k], E, vals))] for k in (('o', 'e') if film['kind'] == 'uniaxial' else ('o',))}
    nk0 = {k: [[float(v.real), float(v.imag)] for v in B.sqrtc(eps(film[k], E))] for k in nk}
    dd = [float(v) for v in x[npf:npf + len(sets)]]
    return dict(run=name, refit=True, bounded=bounded, physical=physical, sets=[s['sh'] for s in sets], d=dd,
                d_on_bound=[on_bound(s['sh'], v, bounded) for s, v in zip(sets, dd)],
                d_ci90=[float(v) for v in ci[npf:npf + len(sets)]],
                f=[float(x[npf + len(sets) + gi.index(j)]) if j in gi else None for j in range(len(sets))],
                mse={s['sh']: float(M.mse_ncs(per(x, j), npf / len(sets) + 1 + (1 if j in gi else 0))) for j, s in enumerate(sets)},
                mse_all=float(M.mse_ncs(r.fun, x.size)), n_fit_film=npf, nfev=int(r.nfev),
                film_params=[dict(axis=k, name=p['name'], old=p['v'], new=float(v), lo=p['lo'], hi=p['hi'],
                                  at_bound=bool(min(abs(v - p['lo']), abs(p['hi'] - v)) < 1e-4 * max(1.0, abs(p['hi']))), ci90=float(c))
                             for (k, p), v, c in zip(fp, x[:npf], ci[:npf])],
                nk_new=nk, nk_old=nk0, nk_wl=wl.tolist())


def main():
    from multiprocessing import Pool
    extra = sys.argv[1:] == ['extra']
    if extra:     # second pass: physically bounded refit, Bphen on the NPB/TAPC oxide convention
        jobs = [('NPB + TAPC glass', True, True, True), ('NPB as saved', True, True, True),
                ('Bphen, NTVE 14.5 A', False, True), ('Bphen, NTVE 14.5 A', False, False)]
    else:
        jobs = ([('NPB as saved', True, True), ('NPB + TAPC glass', True, True)]
                + [(n, False, b) for n in ('NPB as saved', 'NPB + TAPC glass', 'TAPC', 'Bphen') for b in (True, False)])
    with Pool(4) as pool:
        res = pool.starmap(run, jobs)
    json.dump(res, open(W + ('lab_refit_extra.json' if extra else 'lab_refit.json'), 'w'), indent=1)
    for r in res:
        print('\n== %s, %s, %s   overall MSE %.2f' % (r['run'], ('dispersion REFITTED with the lab flags (%d film parameters)%s' % (r['n_fit_film'], ', physical bounds' if r.get('physical') else '')) if r['refit'] else 'dispersion fixed (old n,k)',
              'd within nominal +-5 %% (10 nm: +-10 %%)' if r['bounded'] else 'd free', r['mse_all']))
        for j, sh in enumerate(r['sets']):
            print('   %s %-5s nominal %2d  d = %6.2f nm %s %s MSE %7.2f' % (sh, 'glass' if is_glass(sh) else 'Si', NOMINAL[sh], r['d'][j],
                  'ON BOUND' if r['d_on_bound'][j] else '        ',
                  ('f = %.2f' % r['f'][j]) if r['f'][j] is not None else '        ', r['mse'][sh]))
        if r['refit']:
            for k in r['nk_new']:
                print('   n,k %s old -> new: ' % k + '  '.join('%.0f: %.3f+%.3fj -> %.3f+%.3fj' % (w, *a, *b) for w, a, b in zip(r['nk_wl'], r['nk_old'][k], r['nk_new'][k])))
            nb = [p for p in r['film_params'] if p['at_bound']]
            print('   film parameters on a bound after the refit: %d of %d' % (len(nb), r['n_fit_film']) + ('  (' + ', '.join('%s %s=%.3g' % (p['axis'], p['name'], p['new']) for p in nb) + ')' if nb else ''))


if __name__ == '__main__':
    main()

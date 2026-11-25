"""The lab's own Bphen tooling model (snapshot BPhen_S12_tooling, model 'Bphen_S5'):

  substrate  EMA(Si_JAW, 'no5_glass_sub_genosc'), EMA % (Mat 2) = 0 -> Si, 100 -> glass
  layer 1    NTVE_JAW, 0 A
  layer 2    Bphen Gen-Osc: Einf 0.6616 + 3 Gaussians + UV pole 186.07 eV^2 @ 11 eV
  two samples fitted together: 54.82 nm and 13.80 nm

Its Bphen and glass are read out of the snapshot's _FitLog and used as they are,
to see whether the lab's own model reads the Bphen 10 nm pieces differently."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, re, json, zlib, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_ce_genosc as C, tl_genosc as GO, tl_glass as GL, tl_glass_global as GG, tl_bphen10 as B
import jnk_model as M, ce_osc as osc
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
HC = 1239.841984
SNAP_MOD = _os.path.join(ELLIPS_OUT, 'S5_template.mod')     # written by tl_s5_extract


def _num(rx, t):
    return float(re.search(rx, t).group(1))


def read_s5(path=SNAP_MOD):
    t = open(path, 'rb').read().decode('utf-8-sig')
    film = re.search(r"start_Layer2(.*?)end_Layer2", t, re.S).group(1)
    fb = re.search(r"start_Gen-Osc Fit Parms\r\n(.*?)end_Gen-Osc Fit Parms", film, re.S).group(1)
    einf = float(fb.strip().split('\t')[1])
    g = []
    for blk in fb.split("'Gaussian'")[1:]:
        v = {m.group(2): float(m.group(1)) for m in re.finditer(r"([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'(\w+)'", blk)}
        g.append((v['Amp'], v['Br'], v['En']))
    pole = (_num(r"([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'UV Pole Amp\.'", film),
            _num(r"([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'UV Pole En\.'", film))
    sub = re.search(r"start_no5_glass_sub_genosc Fit Parms(.*?)end_no5_glass_sub_genosc Permanent Poles", t, re.S).group(1)
    ge = float(sub.strip().split('\n')[0].split('\t')[1])
    tl = {k: _num(r"([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'%s'" % k, sub) for k in ('Amp', 'Br', 'Eo', 'Eg')}
    gpole = (_num(r"([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'UV Pole Amp\.'", sub),
             _num(r"([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'UV Pole En\.'", sub),
             _num(r"([-\d.eE+]+)\t[TF]\t[^\t]+\t[^\t]+\tF\t'IR Pole Amp\.'", sub))
    return dict(bphen=dict(einf=einf, gauss=g, pole=pole), glass=dict(einf=ge, tl=tl, pole=gpole))


def eps_s5_bphen(s5, wl):
    E = HC / np.asarray(wl, float); b = s5['bphen']
    e = b['einf'] + b['pole'][0] / (b['pole'][1]**2 - E**2) + 0j
    for a, br, en in b['gauss']:
        x, y = osc.gaussian(E, a, br, en); e = e + x + 1j * y
    return e


def n_s5_glass(s5, wl):
    """CompleteEASE pole convention: A / (En^2 - E^2); the IR pole sits at E = 0."""
    E = HC / np.asarray(wl, float); g = s5['glass']
    a, b = osc.tauc_lorentz(E, g['tl']['Amp'], g['tl']['Br'], g['tl']['Eo'], g['tl']['Eg'])
    e = g['einf'] + g['pole'][0] / (g['pole'][1]**2 - E**2) - g['pole'][2] / E**2 + a + 1j * b
    return B.sqrtc(e)


def fit_si(sheet, disp, d_ox, s5, with_rough=False, n_start=12):
    wl, P, D, ang = L.se(sheet, 245.0, 1000.0, 2); meas = L.ncs(P, D)
    air = np.ones(len(wl), complex); ox = C.N_ce('NTVE_JAW', wl)**2; si = C.N_ce('SI_JAW', wl)**2
    if disp == 'S5':
        ef = eps_s5_bphen(s5, wl)
    else:
        CE = json.load(open(W + 'ce_genosc.json')); ef = B.eps_bphen(np.array(CE['Bphen']['with_roughness']['p']), wl)[0]
    def stack(x):
        lay, ds = [air], []
        if with_rough and x[1] > 1e-4:
            lay.append(B.bruggeman_eps(ef, air, 0.5)); ds.append(x[1])
        lay += [ef]; ds += [x[0]]
        if d_ox > 0:
            lay.append(ox); ds.append(d_ox)
        return lay + [si], ds
    def res(x):
        lay, ds = stack(x)
        N, C_, S_ = B.ncs_regions(wl, ang, [(1.0, lay, lay, ds)], GG.CONV['normalised'])
        return np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    lo = [0.0] + ([0.0] if with_rough else []); hi = [120.0] + ([20.0] if with_rough else [])
    best = None
    for d0 in np.linspace(1, 70, n_start):
        x0 = [d0] + ([1.0] if with_rough else [])
        r = least_squares(res, x0, bounds=(lo, hi))
        if best is None or r.cost < best.cost: best = r
    return [float(v) for v in best.x], float(M.mse_ncs(best.fun, len(best.x)))


def fit_glass(sheet, disp, glass, s5, n_start=10):
    """d and the collected backside fraction f free; glass either the lab's no5 or this set's own."""
    wl, P, D, ang = L.se(sheet, 300.0, 1000.0, 3); meas = L.ncs(P, D)
    if disp == 'S5':
        Nf = B.sqrtc(eps_s5_bphen(s5, wl))
    else:
        CE = json.load(open(W + 'ce_genosc.json')); Nf = GO.N_of(np.array(CE['Bphen']['pick']['p']), wl)
    if glass == 'no5':
        ng = n_s5_glass(s5, wl)
    elif glass == 'no5 k=0':                     # its index only, UV absorption dropped
        ng = n_s5_glass(s5, wl).real + 0j
    else:
        GJ = json.load(open(W + 'ce_glass.json')); A, Bc = GJ['cauchy']; ng = A + Bc / (wl / 1000)**2 + 0j
    def res(x):
        N, C_, S_ = GG.norm(*GL.ncs_glass(wl, Nf, x[0], ng, ang, x[1]))
        return np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    best = None
    for d0 in np.linspace(0.5, 70, n_start):
        r = least_squares(res, [d0, 0.5], bounds=([0.0, 0.0], [150.0, 1.0]))
        if best is None or r.cost < best.cost: best = r
    return [float(v) for v in best.x], float(M.mse_ncs(best.fun, 2))


def main():
    from multiprocessing import Pool
    s5 = read_s5()
    print('lab Bphen (S5): Einf %.4f, Gaussians %s, UV pole %.2f eV^2 @ %.2f eV'
          % (s5['bphen']['einf'], ', '.join('%.0f nm' % (HC / g[2]) for g in s5['bphen']['gauss']), *s5['bphen']['pole']))
    wl = np.array([350., 450., 550., 633.])
    CE = json.load(open(W + 'ce_genosc.json')); GJ = json.load(open(W + 'ce_glass.json'))
    nb5 = B.sqrtc(eps_s5_bphen(s5, wl)); nb = GO.N_of(np.array(CE['Bphen']['with_roughness']['p']), wl)
    print('Bphen n,k  lab S5 vs this set: ' + '  '.join('%3.0f: %.3f+%.3fj / %.3f+%.3fj' % (w, a.real, a.imag, b.real, b.imag) for w, a, b in zip(wl, nb5, nb)))
    ng5 = n_s5_glass(s5, wl); A, Bc = GJ['cauchy']
    print('glass n  lab no5 vs this set: ' + '  '.join('%3.0f: %.4f / %.4f' % (w, a.real, A + Bc / (w / 1000)**2) for w, a in zip(wl, ng5)))
    si_jobs = [(sh, disp, dox, s5, rough) for sh in ('2_21', '2_22', '2_23', '2_24')
               for disp, dox, rough in (('S5', 0.0, False), ('S5', C.D_OX, False), ('S5', 0.0, True), ('this set', C.D_OX, False))]
    gl_jobs = [(sh, disp, g, s5) for sh in ('2_31', '2_32', '2_33', '2_34') for disp, g in (('S5', 'no5'), ('S5', 'no5 k=0'), ('S5', 'own'), ('this set', 'own'))]
    with Pool(4) as pool:
        rs = pool.starmap(fit_si, si_jobs); rg = pool.starmap(fit_glass, gl_jobs)
    print('\nSi pieces')
    for (sh, disp, dox, _, rough), (x, mse) in zip(si_jobs, rs):
        print('  %s  Bphen %-8s NTVE %4.1f A %-6s  d = %6.2f nm%s   MSE %6.2f'
              % (sh, disp, dox * 10, '+rough' if rough else '', x[0], ('  rough %.2f' % x[1]) if rough else '', mse))
    print('\nglass pieces (backside fraction free)')
    for (sh, disp, g, _), (x, mse) in zip(gl_jobs, rg):
        print('  %s  Bphen %-8s glass %-7s  d = %6.2f nm  f = %.2f   MSE %6.2f' % (sh, disp, g, x[0], x[1], mse))
    json.dump(dict(si=[dict(sheet=j[0], disp=j[1], d_ox=j[2], rough=j[4], x=r[0], mse=r[1]) for j, r in zip(si_jobs, rs)],
                   glass=[dict(sheet=j[0], disp=j[1], glass=j[2], x=r[0], mse=r[1]) for j, r in zip(gl_jobs, rg)]),
              open(W + 's5_fits.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

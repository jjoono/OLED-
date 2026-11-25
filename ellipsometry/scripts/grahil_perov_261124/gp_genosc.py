"""The two dispersions as compact CompleteEASE Gen-Osc materials, for .mod files.

gp_grahil / gp_perov describe GraHIL and the dense perovskite with a fixed
Gaussian basis (32 / 44 Gaussians): model-free, but nothing one would edit in a
Gen-Osc dialog.  Here each is refitted AGAINST THE DATA as

    eps(E) = Einf + A_uv / (E_uv^2 - E^2) + sum of a few free Gaussians   (CompleteEASE definitions, lib/ce_osc)

starting from a least-squares match to the basis result.

  GraHIL        layer1, 200-1000 nm; Einf free (with Einf = 1 the pole alone is too dispersive for this film:
                MSE 12-14), UV pole at 12 eV, IR pole (the free-carrier dispersion of PEDOT: a Gaussian
                put there slides below the data range), 3-4 Gaussians, thickness free
  perovskite    layer2, 205-1000 nm, the structure of gp_perov 'base' (roughness / dense / lower / GraHIL),
                all structure parameters free again; Einf and a UV pole of free energy (6.2-12 eV), 7-9 Gaussians.
                Broad Gaussians leave an absorption tail below the gap (k = 0.007 at 633 nm without a guard),
                which the transparent film does not have: eps2 above PEN_FROM nm is penalised.
                GraHIL underneath is the Gen-Osc refit above (NG_GRAHIL), at the 108.3 nm of gp_perov.
                The scattering extinction of the lower layer is held constant below SC_HOLD nm, as in the
                delivered n,k tables, so that the tabulated lower layer of the .mod is what was fitted.

-> genosc.json (parameters, 90 % intervals, bound hits, MSE next to the basis MSE, n,k differences)."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gp_load as G, gp_model as P, gp_grahil as GH, gp_perov as PV
import ce_osc as osc
from tmm_fix import tmm
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
SC_HOLD = 450.0
PEN_W, PEN_FROM = 5.0, 580.0          # weight on the dense perovskite's eps2 below the gap (an N,C,S residual is ~0.01); from this wavelength up
E_POLE_G = 12.0                      # GraHIL: pole energy fixed, as in gp_grahil
# Gaussians: start (Amp, Br, En), En window, Br window
G_GRAHIL = [((0.025, 2.5, 1.5), (0.6, 3.2), (0.3, 4.0)), ((0.015, 2.0, 4.9), (3.2, 5.25), (0.2, 4.0)),
            ((0.06, 0.45, 5.5), (5.25, 5.9), (0.1, 2.0)), ((0.2, 0.85, 6.5), (5.9, 7.5), (0.1, 2.5))]
G_PICK = {3: (0, 2, 3), 4: (0, 1, 2, 3)}
NG_GRAHIL = 4                        # the GraHIL refit that goes into the perovskite stack and into the .mod files
G_PEROV = [((0.40, 0.075, 2.333), (2.28, 2.37), (0.03, 0.20)), ((0.50, 0.14, 2.42), (2.375, 2.48), (0.05, 0.40)),
           ((0.45, 0.24, 2.55), (2.49, 2.70), (0.10, 0.60)), ((1.1, 0.45, 2.87), (2.72, 3.05), (0.15, 1.00)),
           ((1.2, 0.50, 3.26), (3.08, 3.50), (0.15, 0.90)), ((2.2, 0.60, 3.75), (3.52, 3.98), (0.20, 1.00)),
           ((2.8, 0.75, 4.30), (4.00, 4.65), (0.25, 1.30)), ((1.5, 0.80, 5.05), (4.70, 5.45), (0.30, 2.00)),
           ((2.1, 1.1, 5.92), (5.50, 6.60), (0.30, 3.00))]
P_PICK = {7: (0, 1, 2, 3, 4, 6, 8), 8: (0, 1, 2, 3, 4, 5, 6, 8), 9: tuple(range(9))}
HEAD_G = dict(einf=(1.3, 1.0, 2.2), pole_amp=(80.0, 0.0, 400.0), ir_amp=(0.01, 0.0, 1.0))    # start, lo, hi
HEAD_P = dict(einf=(2.0, 1.0, 4.0), pole_amp=(5.0, 0.0, 400.0), pole_en=(7.0, 6.2, 12.0))


def eps_genosc(E, einf, pa, pe, gs, ir=0.0):
    """CompleteEASE poles: A / (En^2 - E^2); the IR pole sits at E = 0"""
    e = einf + pa / (pe**2 - E**2) - ir / E**2 + 0j
    for a, br, en in gs:
        e1, e2 = osc.gaussian(E, a, br, en); e = e + e1 + 1j * e2
    return e


def n_of(E, einf, pa, pe, gs, ir=0.0):
    n = np.sqrt(eps_genosc(E, einf, pa, pe, gs, ir)); return np.where(n.imag < 0, -n, n)


def gbounds(spec, amax):
    lo, hi, x0 = [], [], []
    for (a, br, en), (e_lo, e_hi), (b_lo, b_hi) in spec:
        lo += [0.0, b_lo, e_lo]; hi += [amax, b_hi, e_hi]; x0 += [a, br, en]
    return np.array(lo), np.array(hi), np.array(x0)


def head_bounds(head):
    return np.array([v[1] for v in head.values()]), np.array([v[2] for v in head.values()]), np.array([v[0] for v in head.values()])


def disp(E, xh, xg, head, pe_fixed):
    """xh: the free head parameters in the order of head (einf, pole_amp[, pole_en][, ir_amp])"""
    h = dict(zip(head, xh)); return n_of(E, h['einf'], h['pole_amp'], h.get('pole_en', pe_fixed), np.asarray(xg).reshape(-1, 3), h.get('ir_amp', 0.0))


def n_from_record(rec, wl, head, pe_fixed=None):
    """index of a fitted record (fit_grahil / fit_perov) on any wavelength grid"""
    p = rec['parms']; ng = rec['ng']
    return disp(P.HC / np.asarray(wl, float), [p[k]['v'] for k in head], [p['%s%d' % (k, i + 1)]['v'] for i in range(ng) for k in ('Amp', 'Br', 'En')], head, pe_fixed)


def to_basis(E, eps_target, spec, amax, head, pe_fixed):
    """Head + Gaussians that best reproduce a tabulated eps: the start for the data fit."""
    gl_, gh_, g0 = gbounds(spec, amax); hl, hh, h0 = head_bounds(head); nh = len(head)
    lo = np.r_[hl, gl_]; hi = np.r_[hh, gh_]; x0 = np.r_[h0, g0]
    def res(x):
        e = disp(E, x[:nh], x[nh:], head, pe_fixed)**2 - eps_target; return np.r_[e.real, e.imag]
    return least_squares(res, np.clip(x0, lo + 1e-9, hi - 1e-9), bounds=(lo, hi), x_scale='jac', max_nfev=600).x, lo, hi


def summary(r, lo, hi, names, n_meas):
    J = r.jac; s2 = np.sum(r.fun[:n_meas]**2) / max(n_meas - r.x.size, 1)
    ci = 1.645 * np.sqrt(np.clip(np.diag(np.linalg.pinv(J.T @ J) * s2), 0, None))
    at = [bool(min(x - l, h - x) < 1e-6 * max(h - l, 1e-12)) for x, l, h in zip(r.x, lo, hi)]
    return {k: dict(v=float(x), ci=float(c), lo=float(l), hi=float(h), at_bound=b) for k, x, c, l, h, b in zip(names, r.x, ci, lo, hi, at)}


# ------------------------------------------------------------------ GraHIL
def fit_grahil(ng):
    S = GH.setup(); wl = S['wl']; E = P.HC / wl; meas = P.data_ncs(S['d'])
    Nb, d_b = PV.grahil_N(wl); spec = [G_GRAHIL[i] for i in G_PICK[ng]]; nh = len(HEAD_G)
    xs, lo, hi = to_basis(E, Nb**2, spec, 5.0, HEAD_G, E_POLE_G)
    lo = np.r_[20.0, lo]; hi = np.r_[300.0, hi]; x0 = np.r_[d_b, xs]
    def model(x):
        n = disp(E, x[1:1 + nh], x[1 + nh:], HEAD_G, E_POLE_G)
        return P.ncs(wl, S['ang'], [S['air'], n, S['ox'], S['si']], [x[0], PV.D_OX]), n
    r = least_squares(lambda x: P.resid(model(x)[0], meas), np.clip(x0, lo + 1e-9, hi - 1e-9), bounds=(lo, hi), x_scale='jac', max_nfev=400)
    names = ['d'] + list(HEAD_G) + ['%s%d' % (k, i + 1) for i in range(ng) for k in ('Amp', 'Br', 'En')]
    n = model(r.x)[1]; q = (wl >= 250)
    return dict(ng=ng, mse=float(P.mse_ncs(r.fun, len(r.x))), parms=summary(r, lo, hi, names, r.fun.size), pole_en=E_POLE_G, d_ox=PV.D_OX,
                dn_max=float(np.abs(n.real - Nb.real)[q].max()), dk_max=float(np.abs(n.imag - Nb.imag)[q].max()),
                nk={str(w): [float(np.interp(w, wl, n.real)), float(np.interp(w, wl, n.imag))] for w in (300, 450, 550, 633, 800, 1000)})


# ------------------------------------------------------------------ perovskite
def fit_perov(ng, g_rec):
    """g_rec: the GraHIL record whose Gen-Osc form is the GraHIL layer here (so the .mod is what was fitted)"""
    S = PV.setup(); wl = S['wl']; E = S['E']; S['Ng'] = n_from_record(g_rec, wl, HEAD_G, E_POLE_G)
    base = [r for r in json.load(open(W + 'perov_fits.json'))['variants'] if r['variant'] == 'base'][0]
    Nb = PV.nk_of(base, wl); spec = [G_PEROV[i] for i in P_PICK[ng]]; nh = len(HEAD_P)
    xs, lo, hi = to_basis(E, Nb**2, spec, 60.0, HEAD_P, None)
    names = PV.STRUCT + list(HEAD_P) + ['%s%d' % (k, i + 1) for i in range(ng) for k in ('Amp', 'Br', 'En')]; ns = len(PV.STRUCT)
    lo = np.r_[[PV.S_LO[k] for k in PV.STRUCT], lo]; hi = np.r_[[PV.S_HI[k] for k in PV.STRUCT], hi]
    x0 = np.clip(np.r_[[base['struct'][k] for k in PV.STRUCT], xs], lo + 1e-9, hi - 1e-9)
    def split(x):
        v = dict(zip(PV.STRUCT, x[:ns])); v['sc_hold'] = SC_HOLD
        return v, disp(E, x[ns:ns + nh], x[ns + nh:], HEAD_P, None)
    sub = wl >= PEN_FROM
    def res(x):
        v, n = split(x); a, b = PV.residuals(S, v, n); return np.r_[a, b, PEN_W * (n**2).imag[sub]]
    r = least_squares(res, x0, bounds=(lo, hi), x_scale='jac', max_nfev=400)
    v, n = split(r.x); a, b = PV.residuals(S, v, n); rr = a.reshape(3, len(wl), 3)
    band = lambda p, q: float(1e3 * np.sqrt(np.mean(rr[:, (wl >= p) & (wl < q), :]**2)))
    q = (wl >= 450)
    return dict(ng=ng, mse=float(P.mse_ncs(a, len(r.x))), mse_basis=base['mse_ncs'], rms_R=float(1e3 * np.sqrt(np.mean((b / PV.WR)**2))),
                bands={'205-300': band(205, 300), '300-450': band(300, 450), '450-570': band(450, 570), '570-1000': band(570, 1001)},
                parms=summary(r, lo, hi, names, a.size + b.size), sc_hold=SC_HOLD, k_subgap_max=float(n.imag[sub].max()), d_g=PV.D_G, d_ox=PV.D_OX,
                dn_max_450_1000=float(np.abs(n.real - Nb.real)[q].max()), dk_max_450_1000=float(np.abs(n.imag - Nb.imag)[q].max()),
                nk={str(w): [float(np.interp(w, wl, n.real)), float(np.interp(w, wl, n.imag))] for w in (300, 400, 450, 500, 520, 530, 540, 550, 600, 633, 800, 1000)})


def job(a):
    return (a[0], a[1], fit_grahil(a[1]) if a[0] == 'grahil' else fit_perov(a[1], a[2]))


def main():
    from multiprocessing import Pool
    with Pool(4) as pool:
        res = pool.map(job, [('grahil', n) for n in G_PICK], chunksize=1)
        g_rec = [r for k, n, r in res if n == NG_GRAHIL][0]
        res += pool.map(job, [('perov', n, g_rec) for n in P_PICK], chunksize=1)
    gb = [v for v in json.load(open(W + 'grahil_fits.json'))['variants'] if v['variant'] == 'iso'][0]
    out = dict(grahil={str(n): r for k, n, r in res if k == 'grahil'}, perov={str(n): r for k, n, r in res if k == 'perov'}, grahil_basis_mse=gb['mse'], grahil_basis_d=gb['d'])
    json.dump(out, open(W + 'genosc.json', 'w'), indent=1)
    tri = lambda p, n: '  '.join('(%.4g, %.3f, %.3f)' % tuple(p['%s%d' % (k, i + 1)]['v'] for k in ('Amp', 'Br', 'En')) for i in range(int(n)))
    print('GraHIL (basis: MSE %.2f, d %.2f nm)' % (gb['mse'], gb['d']))
    for n, r in out['grahil'].items():
        p = r['parms']; print('  %s Gaussians: MSE %.2f  d %.2f +- %.2f nm  n,k(550) %.4f %.4f  max |dn| %.4f |dk| %.4f vs basis (250-1000 nm)  at bound: %s'
                              % (n, r['mse'], p['d']['v'], p['d']['ci'], *r['nk']['550'], r['dn_max'], r['dk_max'], [k for k in p if p[k]['at_bound']]))
        print('     Einf %.4f  UV pole %.2f @ %.1f eV  IR pole %.5f +- %.5f   ' % (p['einf']['v'], p['pole_amp']['v'], r['pole_en'], p['ir_amp']['v'], p['ir_amp']['ci']) + tri(p, n))
    print('dense perovskite in the gp_perov structure')
    for n, r in out['perov'].items():
        p = r['parms']; print('  %s Gaussians: MSE %.2f (basis %.2f)  bands %s  rms R %.1f | ' % (n, r['mse'], r['mse_basis'], ' '.join('%.1f' % v for v in r['bands'].values()), r['rms_R'])
                              + ' '.join('%s=%.3g' % (k, p[k]['v']) for k in PV.STRUCT))
        print('     n,k 540: %.3f %.3f  633: %.3f %.4f  max k above %g nm %.4f  max |dn| %.3f |dk| %.3f vs basis (450-1000 nm)  at bound: %s' % (*r['nk']['540'], *r['nk']['633'], PEN_FROM, r['k_subgap_max'], r['dn_max_450_1000'], r['dk_max_450_1000'], [k for k in p if p[k]['at_bound']]))
        print('     Einf %.4f  pole %.2f @ %.2f eV   ' % (p['einf']['v'], p['pole_amp']['v'], p['pole_en']['v']) + tri(p, n))


if __name__ == '__main__':
    main()

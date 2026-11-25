"""Gen-Osc fits in CompleteEASE's own frame, for writing .mod files.

Substrate and native oxide are the SI_JAW / NTVE_JAW tables decoded out of the
lab's own .mod files (all four supplied files carry identical tables), and the
oxide thickness is the lab's measured value for O2-plasma-treated wafers,
read from Nativeoxide_180sO2plasma_Si_wafer.mod.  A .mod written from these
parameters therefore reproduces the data on Generate alone.

Dispersion rules (see tl_genosc.py): Einf = 1 fixed, UV pole 6-12 eV, Gaussians
with amplitude >= 0 and FWHM 0.10-1.20 eV centred above the absorption onset -
plus a penalty on eps2 below the onset, because a broad band centred above the
gap can still leak a visible tail, and these films are transparent there."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, re, json, base64, gzip, zlib, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_genosc as G
import jnk_model as M
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
TEMPLATE = _os.path.join(ELLIPS_DATA, 'ADNC_ITO_GenOSC_v2_claude.mod')
OXMOD = _os.path.join(ELLIPS_DATA, 'Nativeoxide_180sO2plasma_Si_wafer.mod')
HC = 1239.841984
PEN_W = 5.0                # weight on sub-gap eps2 (an NCS residual is ~0.002)
NGS = {'NPB': (4, 5, 6), 'TAPC': (3, 4, 5), 'Bphen': (3, 4, 5)}


# ------------------------------------------------------------------ the CE frame
def _decode(b64):
    raw = base64.b64decode(b64)
    try:
        d = gzip.decompress(raw)
    except Exception:
        d = zlib.decompress(raw)
    return np.frombuffer(d, '>f4').astype(float)


def ce_tables(path=TEMPLATE):
    t = open(path, 'rb').read().decode('utf-8-sig')
    out = {}
    for name in ('SI_JAW', 'NTVE_JAW'):
        m = re.search(r"'Layer'\t'%s'.*?start_Mat Table(.*?)end_Mat Table" % name, t, re.S)
        arr = {tag: _decode(re.search(r"start_%s Array\r\n\t+'([^']+)'" % tag, m.group(1)).group(1))
               for tag in ('Wvl', 'e1', 'e2')}
        out[name] = (arr['Wvl'] / 10.0, arr['e1'], arr['e2'])       # Angstrom -> nm
    return out


_T = None
def N_ce(name, wl):
    global _T
    if _T is None:
        _T = ce_tables()
    w, e1, e2 = _T[name]; o = np.argsort(w)
    e = np.interp(wl, w[o], e1[o]) + 1j * np.interp(wl, w[o], e2[o])
    N = np.sqrt(e.astype(complex)); return np.where(N.imag < 0, -N, N)


def lab_oxide_nm(path=OXMOD):
    t = open(path, 'rb').read().decode('utf-8-sig')
    return float(re.search(r"\n\t\t([\d.eE+-]+)\t[TF]\t[^\t]+\t[^\t]+\t[TF]\t'Native Oxide'", t).group(1)) / 10.0


D_OX = lab_oxide_nm() if _os.path.exists(OXMOD) else 3.368


# ------------------------------------------------------------------ fits
def data_for(shs, step=2):
    out = []
    for sh in shs:
        wl, P, D, ang = L.se(sh, 245.0, 1000.0, step)
        out.append((wl, L.ncs(P, D), ang, np.ones(len(wl), complex), N_ce('NTVE_JAW', wl), N_ce('SI_JAW', wl)))
    return out


def stack(p, d, rough, wl, air, ox, si, d_ox):
    Nf = G.N_of(p, wl)
    if rough > 1e-6:
        return [air, M.bruggeman(Nf, air, 0.5), Nf, ox, si], [rough, d, d_ox]
    return [air, Nf, ox, si], [d, d_ox]


def residual(x, data, npar, name, with_rough, d_ox):
    p = x[:npar]; nsh = len(data); out = []
    for j, (wl, meas, ang, air, ox, si) in enumerate(data):
        rough = x[npar + nsh + j] if with_rough else 0.0
        ls, ds = stack(p, x[npar + j], rough, wl, air, ox, si, d_ox)
        n_, c_, s_ = M.se_ncs(wl, ls, ds, ang)
        out += [(n_ - meas[0]).ravel(), (c_ - meas[1]).ravel(), (s_ - meas[2]).ravel()]
    if name:                                     # transparent below the onset
        eon = HC / G.ONSET_NM[name]
        Eg = np.linspace(1.24, eon - 0.10, 40)
        out.append(PEN_W * G.eps_genosc(p, HC / Eg).imag)
    return np.concatenate(out)


def fit_dispersion(name, ng, d0, with_rough=False, d_ox=None, n_start=10):
    d_ox = D_OX if d_ox is None else d_ox
    shs = G.MAT[name]; data = data_for(shs)
    lo, hi = G.bounds(name, ng); npar = len(lo); nsh = len(shs)
    LO = np.concatenate([lo, [1.0] * nsh] + ([[0.0] * nsh] if with_rough else []))
    HI = np.concatenate([hi, [200.0] * nsh] + ([[15.0] * nsh] if with_rough else []))
    rng = np.random.default_rng(zlib.crc32(repr((name, ng, with_rough)).encode())); best = None
    for p0 in G.starts(name, ng, rng, n_start):
        x0 = np.concatenate([p0, d0] + ([[0.5] * nsh] if with_rough else []))
        x0 = np.clip(x0, LO + 1e-6, HI - 1e-6)
        try:
            r = least_squares(residual, x0, args=(data, npar, name, with_rough, d_ox), bounds=(LO, HI),
                              x_scale='jac', max_nfev=2500)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    x = best.x
    J = best.jac; res = best.fun; s2 = np.sum(res**2) / max(res.size - x.size, 1)
    ci = 1.645 * np.sqrt(np.clip(np.diag(np.linalg.pinv(J.T @ J) * s2), 0, None))
    per = []
    for j in range(nsh):
        xx = np.concatenate([x[:npar], [x[npar + j]]] + ([[x[npar + nsh + j]]] if with_rough else []))
        per.append(float(M.mse_ncs(residual(xx, [data[j]], npar, None, with_rough, d_ox), npar / nsh + 1)))
    at = [i for i in range(npar) if min(abs(x[i] - LO[i]), abs(HI[i] - x[i])) < 1e-4 * max(1.0, abs(HI[i]))]
    return dict(name=name, ng=ng, rough=with_rough, d_ox=d_ox, p=[float(v) for v in x[:npar]],
                d=[float(v) for v in x[npar:npar + nsh]],
                r=[float(v) for v in x[npar + nsh:]] if with_rough else [0.0] * nsh,
                ci=[float(v) for v in ci], mse_per=per, mse=float(np.mean(per)), at_bound=at, sheets=shs)


def fit_thickness(sh, p, d_ox=None, with_rough=False):
    """A piece against a fixed dispersion: thickness (and roughness) only."""
    d_ox = D_OX if d_ox is None else d_ox
    data = data_for([sh]); npar = len(p)
    def res(x):
        return residual(np.concatenate([p, x]), data, npar, None, with_rough, d_ox)
    best = None
    for d0 in np.arange(1.0, 120.0, 3.0):
        x0, lo, hi = ([d0, 1.0], [0.2, 0.0], [300.0, 15.0]) if with_rough else ([d0], [0.2], [300.0])
        r = least_squares(res, x0, bounds=(lo, hi))
        if best is None or r.cost < best.cost: best = r
    d = float(best.x[0]); rg = float(best.x[1]) if with_rough else 0.0
    return d, float(M.mse_ncs(best.fun, len(best.x))), rg


def main():
    from multiprocessing import Pool
    CS = json.load(open(W + 'cauchy_si.json'))
    print('CompleteEASE frame: SI_JAW + NTVE_JAW %.2f nm (lab O2-plasma value)' % D_OX)
    jobs = [(n, ng, np.array([CS[s]['d'] for s in G.MAT[n]]), False) for n in G.MAT for ng in NGS[n]]
    jobs += [(n, None, np.array([CS[s]['d'] for s in G.MAT[n]]), True) for n in G.MAT]   # roughness test later
    with Pool(4) as pool:
        res = pool.starmap(fit_dispersion, [j for j in jobs if j[1] is not None])
    out = {}
    for name in G.MAT:
        rs = [r for r in res if r['name'] == name]
        best = min(r['mse'] for r in rs)
        pick = min((r for r in rs if r['mse'] <= 1.10 * best), key=lambda r: r['ng'])
        out[name] = dict(pick=pick, all=rs)
        print('%-6s ' % name + '  '.join('ng=%d MSE %.2f%s' % (r['ng'], r['mse'], ' <' if r is pick else '') for r in rs))
    # roughness: does the data want any, with the picked oscillator count?
    with Pool(3) as pool:
        rr = pool.starmap(fit_dispersion, [(n, out[n]['pick']['ng'], np.array(out[n]['pick']['d']), True, None, 4)
                                           for n in G.MAT])
    for r in rr:
        out[r['name']]['with_roughness'] = r
        print('%-6s with roughness: MSE %.2f (vs %.2f)  roughness %s nm  d %s'
              % (r['name'], r['mse'], out[r['name']]['pick']['mse'], ', '.join('%.2f' % v for v in r['r']),
                 ', '.join('%.2f' % v for v in r['d'])))
    json.dump(out, open(W + 'ce_genosc.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

"""Bphen 10 nm (2_23, 2_24 on Si; 2_33, 2_34 on glass): which physical picture fits?

Every homogeneous single-layer variant tried earlier (any dispersion, free oxide,
void fraction, angle offset) left MSE 34-42.  Three pictures were not tried:

  rough      a continuous Bphen film under a Bruggeman 50 % roughness layer, and
             the same with the film itself porous (Bruggeman void fraction)
  uniaxial   crystallized Bphen is oriented: Gaussians scaled by an order
             parameter S, eps2_o = eps2 (1 - S), eps2_e = eps2 (1 + 2 S), optic
             axis along the normal (S = 0 isotropic, -0.5 lying, 1 standing)
  islands    crystallites larger than the beam's coherence length next to bare
             (or thinly wetted) substrate: the two regions' reflections add
             incoherently, area-weighted - CompleteEASE's MultiModel.  This is the
             one picture a single coherent stack cannot imitate at all.

Si pieces in CompleteEASE's frame (SI_JAW, NTVE_JAW 33.68 A); glass pieces with
the glass and backside fraction shared by the other glass pieces (tl_ce_glass).
The Bphen dispersion is the one fitted to the 50 nm pieces."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')

import sys, json, zlib, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import tl_load as L, tl_genosc as GO, tl_ce_genosc as C, tl_glass as GL, tl_glass_global as GG
import jnk_model as M, ce_osc as osc
from tmm_fix import tmm
from scipy.optimize import least_squares

W = _os.path.join(ELLIPS_OUT, '')
HC = 1239.841984


# ------------------------------------------------------------------ optics
def eps_bphen(p, wl, S=0.0):
    """-> eps_o, eps_e with the Gaussians weighted by the order parameter."""
    E = HC / np.asarray(wl, float)
    base = GO.EINF + p[0] / (p[1]**2 - E**2) + 0j
    go = np.zeros_like(E, complex)
    for i in range((len(p) - 2) // 3):
        a, b = osc.gaussian(E, p[2 + 3 * i], p[3 + 3 * i], p[4 + 3 * i]); go += a + 1j * b
    return base + (1 - S) * go, base + (1 + 2 * S) * go


def sqrtc(e):
    n = np.sqrt(np.asarray(e, complex)); return np.where(n.imag < 0, -n, n)


def tmm_uni(wl, eo, ee, ds, theta):
    """Stack of uniaxial layers, optic axis along z (isotropic: eo = ee).
    Same propagation convention as tmm_fix.tmm; interface scalars cancel in r."""
    th = np.deg2rad(theta); k0 = 2 * np.pi / wl
    s2 = (sqrtc(eo[0]) * np.sin(th))**2
    qs = [sqrtc(o - s2) for o in eo]
    qp = [sqrtc(o * (1 - s2 / e)) for o, e in zip(eo, ee)]
    eta_s = qs; eta_p = [o / q for o, q in zip(eo, qp)]
    n = len(wl)
    def run(eta, q):
        Mx = np.zeros((n, 2, 2), complex); Mx[:, 0, 0] = 1; Mx[:, 1, 1] = 1
        mm = lambda A, B: np.einsum('nij,njk->nik', A, B)
        for i in range(len(eo) - 1):
            r = (eta[i] - eta[i + 1]) / (eta[i] + eta[i + 1])
            I = np.zeros((n, 2, 2), complex); I[:, 0, 0] = 1; I[:, 0, 1] = r; I[:, 1, 0] = r; I[:, 1, 1] = 1
            Mx = mm(Mx, I)
            if i < len(ds):
                ph = k0 * q[i + 1] * ds[i]
                P = np.zeros((n, 2, 2), complex); P[:, 0, 0] = np.exp(-1j * ph); P[:, 1, 1] = np.exp(1j * ph)
                Mx = mm(Mx, P)
        return Mx[:, 1, 0] / Mx[:, 0, 0]
    # p-pol: eta = eps_o / q_p gives the H-field r with the same sign as tmm_fix
    return -run(eta_p, qp), run(eta_s, qs)


def ncs_regions(wl, ang, regions, conv):
    """regions = [(weight, eo_layers, ee_layers, ds), ...] -> N, C, S (incoherent sum)."""
    out = []
    for th in ang:
        Ip = Is = 0.0; X = 0.0
        for w, eo, ee, ds in regions:
            rp, rs = tmm_uni(wl, eo, ee, ds, th)
            Ip = Ip + w * np.abs(rp)**2; Is = Is + w * np.abs(rs)**2; X = X + w * np.conj(rp) * rs
        I = Ip + Is
        out.append(conv((Is - Ip) / I, 2 * X.real / I, 2 * X.imag / I))
    return tuple(np.stack([o[k] for o in out], 1) for k in range(3))


# ------------------------------------------------------------------ models
def si_stack_layers(wl, P, d_ox, film):
    """film = list of (eo, ee, d) from the top down; returns eo, ee, ds for air..Si."""
    air = np.ones(len(wl), complex)
    ox = C.N_ce('NTVE_JAW', wl)**2; si = C.N_ce('SI_JAW', wl)**2
    eo = [air] + [f[0] for f in film] + [ox, si]
    ee = [air] + [f[1] for f in film] + [ox, si]
    return eo, ee, [f[2] for f in film] + [d_ox]


def bruggeman_eps(ea, eb, fa):
    return M.bruggeman(sqrtc(ea), sqrtc(eb), fa)**2


MODELS = {
    #            names                      x0                     lo                   hi
    'film':           (['d'],               [8.0],               [0.1],              [80.0]),
    'film+rough':     (['d', 'r'],          [6.0, 3.0],          [0.0, 0.0],         [80.0, 40.0]),
    'porous+rough':   (['d', 'fb', 'r'],    [10.0, 0.6, 3.0],    [0.0, 0.05, 0.0],   [80.0, 1.0, 40.0]),
    'uniaxial':       (['d', 'S'],          [8.0, 0.2],          [0.1, -0.5],        [80.0, 1.0]),
    'uniaxial+rough': (['d', 'S', 'r'],     [6.0, 0.2, 3.0],     [0.0, -0.5, 0.0],   [80.0, 1.0, 40.0]),
    'islands':        (['a', 'h', 'w'],     [0.3, 30.0, 1.0],    [0.0, 1.0, 0.0],    [1.0, 300.0, 20.0]),
    'islands+rough':  (['a', 'h', 'w', 'r'],[0.3, 30.0, 1.0, 5.0],[0.0, 1.0, 0.0, 0.0],[1.0, 300.0, 20.0, 60.0]),
}


def regions_for(model, x, wl, p, sub):
    """sub(wl, film) -> (eo, ee, ds) for one region on this substrate."""
    v = dict(zip(MODELS[model][0], x)); air = np.ones(len(wl), complex)
    eb_o, eb_e = eps_bphen(p, wl, v.get('S', 0.0))
    def rough_layer(eo, ee, r):
        return (bruggeman_eps(eo, air, 0.5), bruggeman_eps(ee, air, 0.5), r)
    if model.startswith('islands'):
        isl = [(eb_o, eb_e, v['h'])]
        if v.get('r', 0) > 1e-3:
            isl = [rough_layer(eb_o, eb_e, v['r'])] + isl
        bare = [(eb_o, eb_e, v['w'])] if v['w'] > 1e-3 else []
        return [(v['a'], *sub(wl, isl)), (1 - v['a'], *sub(wl, bare))]
    if model.startswith('porous'):
        eb_o = eb_e = bruggeman_eps(eb_o, air, v['fb'])
    film = [(eb_o, eb_e, v['d'])]
    if v.get('r', 0) > 1e-3:
        film = [rough_layer(eb_o, eb_e, v['r'])] + film
    return [(1.0, *sub(wl, film))]


def fit(sheet, model, conv_name='normalised', n_start=24):
    on_glass = sheet.split('_')[1].startswith('3')
    CE = json.load(open(W + 'ce_genosc.json')); p = np.array(CE['Bphen']['with_roughness']['p'])
    wl, P, D, ang = L.se(sheet, 300.0 if on_glass else 245.0, 1000.0, 3 if on_glass else 2)
    meas = L.ncs(P, D)
    conv = GG.CONV[conv_name]
    if on_glass:
        GJ = json.load(open(W + 'ce_glass.json')); A, B = GJ['cauchy']; fback = GJ['f']
        eg = (A + B / (wl / 1000)**2)**2 + 0j
        def sub(wl_, film):
            air = np.ones(len(wl_), complex)
            return [air] + [f[0] for f in film] + [eg], [air] + [f[1] for f in film] + [eg], [f[2] for f in film]
    else:
        def sub(wl_, film):
            return si_stack_layers(wl_, None, C.D_OX, film)
    def model_ncs(x):
        regs = regions_for(model, x, wl, p, sub)
        if on_glass:     # glass: incoherent back surface on every region (isotropic film approx. for the back path)
            out_ = []
            for th in ang:
                Ip = Is = 0.0; X = 0.0
                for w, eo, ee, ds in regs:
                    rp, rs = tmm_uni(wl, eo, ee, ds, th)
                    Ip = Ip + w * np.abs(rp)**2; Is = Is + w * np.abs(rs)**2; X = X + w * np.conj(rp) * rs
                    # back surface, through an isotropic version of the same stack
                    Nf = [sqrtc(e) for e in eo]
                    ng = Nf[-1]
                    thg = np.rad2deg(np.arcsin(np.sin(np.deg2rad(th)) / ng.real))
                    tdp, tds = M.t_amp(wl, Nf, ds, th)
                    tup, tus = M.t_amp(wl, Nf[::-1], ds[::-1], thg)
                    rip, ris = tmm(wl, Nf[::-1], ds[::-1], thg)
                    rbp, rbs = tmm(wl, [ng, np.ones(len(wl), complex)], [], thg)
                    ap = tdp * rbp * tup; as_ = tds * rbs * tus
                    for _ in range(12):
                        Ip = Ip + w * fback * np.abs(ap)**2; Is = Is + w * fback * np.abs(as_)**2
                        X = X + w * fback * np.conj(ap) * as_
                        ap = ap * rip * rbp; as_ = as_ * ris * rbs
                I = Ip + Is
                out_.append(conv((Is - Ip) / I, 2 * X.real / I, 2 * X.imag / I))
            return tuple(np.stack([o[k] for o in out_], 1) for k in range(3))
        return ncs_regions(wl, ang, regs, conv)
    def res(x):
        N, C_, S_ = model_ncs(x)
        return np.concatenate([(N - meas[0]).ravel(), (C_ - meas[1]).ravel(), (S_ - meas[2]).ravel()])
    names, x0, lo, hi = MODELS[model]
    lo, hi = np.array(lo, float), np.array(hi, float)
    rng = np.random.default_rng(zlib.crc32(repr((sheet, model, conv_name)).encode())); best = None
    for k in range(n_start):
        s = np.array(x0, float) if k == 0 else lo + rng.random(len(lo)) * (hi - lo)
        try:
            r = least_squares(res, np.clip(s, lo + 1e-6, hi - 1e-6), bounds=(lo, hi), max_nfev=600)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    rr = best.fun.reshape(3, len(wl), len(ang))
    band = lambda a, b: float(1e3 * np.sqrt(np.mean(rr[:, (wl >= a) & (wl < b), :]**2)))
    return dict(sheet=sheet, model=model, conv=conv_name, x=dict(zip(names, [float(v) for v in best.x])),
                mse=float(M.mse_ncs(best.fun, len(names))),
                bands={'UV<350': band(200, 350), '350-500': band(350, 500), '500-1000': band(500, 1001)})


def main():
    from multiprocessing import Pool
    sheets = ['2_23', '2_24', '2_33', '2_34']
    jobs = [(sh, m) for sh in sheets for m in MODELS]
    with Pool(4) as pool:
        res = pool.starmap(fit, jobs)
    json.dump(res, open(W + 'bphen10.json', 'w'), indent=1)
    for sh in sheets:
        print('== %s (%s)' % (sh, 'glass' if sh.split('_')[1].startswith('3') else 'Si'))
        for r in sorted([r for r in res if r['sheet'] == sh], key=lambda r: r['mse']):
            print('   %-15s MSE %7.2f | UV %6.1f  350-500 %6.1f  500-1000 %6.1f | %s'
                  % (r['model'], r['mse'], r['bands']['UV<350'], r['bands']['350-500'], r['bands']['500-1000'],
                     '  '.join('%s=%.3g' % kv for kv in r['x'].items())))


if __name__ == '__main__':
    main()

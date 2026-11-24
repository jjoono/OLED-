"""Green PeLED (Nature 611, 688): where the photons go, for several descriptions of the perovskite layer.

    glass (incoherent) / FTO 350 / HIL / perovskite 270 / ZADN 45 / LiF 1.2 / Al

Every emission event is followed to its end with the dipole model of gp_dipole:

  air      leaves through the glass into air (with the incoherent reflections between the glass/air face and the stack)
  sub      trapped in the glass (in-plane index between 1 and n_glass)
  Ap       absorbed by the perovskite itself -> can be re-emitted (photon recycling)
  Apar     absorbed elsewhere (FTO, HIL, Al): lost
  S        scattered in the lower, porous part of the perovskite (its extinction k_sc is not absorption)

Scattered light is not lost: it restarts as an isotropic, incoherent source inside the lower part, at
the same wavelength, and is followed again (a geometric series, summed exactly).  Light trapped in the
glass either counts as lost ('sub lost': it leaves the pixel sideways) or keeps bouncing between the
glass/air face and the stack until it is absorbed or scattered there ('sub kept': a laterally infinite
device).  Re-absorbed light is re-emitted with the internal radiative efficiency eta_pl and the
internal spectrum:  EQE = eta_pl * air / (1 - eta_pl * Ap)   (charge balance 1).

The internal spectrum is not assumed: it is the one that makes the model's normal-direction output the
measured EL spectrum (Extended Data Fig. 7c).  The test is the angular distribution (Extended Data
Fig. 7g) and whether the measured EQE (28.9 %) is reachable with eta_pl <= 1.

Emitters: isotropic dipoles, uniform in depth within the dense part and within the lower part;
x_d = share of the first-generation photons born in the dense part ('uniform' = in proportion to the
perovskite volume).  Fates are per photon (each position normalised by its own total rate), which is
right for a high radiative efficiency.  Lateral thickness variation: fates and patterns are averaged
over a Gaussian spread of the perovskite thickness."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gp_model as P, gp_device as Dv, gp_device_check as DC, gp_perov as PV, gp_dipole as DP

W = _os.path.join(ELLIPS_OUT, '')
WL = np.arange(505.0, 585.1, 2.5)
TH = np.arange(0.0, 81.0, 10.0)
N_GLASS, D_FTO, D_ZADN, D_LIF, D_PER = 1.5, 350.0, 45.0, 1.2, 270.0
OFFS = np.arange(-60.0, 60.1, 10.0)             # perovskite-thickness offsets (nm) from which the spreads are built
NZ = 6
EQE_MEAS, EQE_MEAN, PLQE_MEAS = 0.289, 0.252, 0.46        # champion and mean (40 devices) EQE; external PL yield of the film on quartz (Extended Data Fig. 6)
KEYS = ('air', 'sub', 'Ap', 'Apar', 'S')


# ------------------------------------------------------------------ structure
def dense_index(fit, table='basis', n_scale=1.0, k_edge_scale=1.0):
    """dense perovskite on WL.  table 'stitched': the luminescence-stitched band edge of gp_stitch.py.
    k_edge_scale multiplies k (a sensitivity test of the absorption just above the edge)."""
    if table == 'stitched':
        a = np.loadtxt(_os.path.join(ELLIPS_OUT, 'grahil_perov', 'Perov_dense_nk_stitched.csv'), delimiter=',', skiprows=1)
        Np = np.interp(WL, a[:, 0], a[:, 1]) + 1j * np.interp(WL, a[:, 0], a[:, 2])
    else:
        Np = PV.nk_of(fit, WL)
    return Np.real * n_scale + 1j * Np.imag * k_edge_scale


def build(kind, fit=None, d_hil=75.0, d_per=D_PER, d_zadn=D_ZADN, d_fto=D_FTO, n_zadn=2.0, n_fto=1.9 + 0.005j, undercoat=False, f_b=None, ksc_scale=1.0,
          dense_share=None, n_scale=1.0, k_edge_scale=1.0, table='basis', film=False, extra=0.0):
    """-> n [media, wl], d [layers], roles [layers] ('par' | 'perov' | 'lower' | 'out'), share_sc [wl], src = {'L': layer or None, 'D': layer},
    vol = {'L', 'D'} perovskite volume (nm) of each part.
    film=True: the bare film on a transparent substrate (substrate / perovskite / air), for the photoluminescence yield."""
    one = np.ones(len(WL), complex); al = Dv.n_al(WL)
    hil = 1.5 * one if kind == 'paper' else PV.grahil_N(WL)[0]
    if kind == 'paper':
        per = [(2.3 + 1j * np.interp(WL, *DC.K_PAPER), d_per + extra, 'perov')]; share = np.zeros(len(WL)); vol = dict(D=d_per)
    else:
        Np = dense_index(fit, table, n_scale, k_edge_scale)
        if kind == 'dense':
            per = [(Np, d_per + extra, 'perov')]; share = np.zeros(len(WL)); vol = dict(D=d_per)
        else:
            v = fit['struct']; fb = v['f_b'] if f_b is None else f_b; host = v.get('n_m', 1.0) * one
            sh = v['d_t'] / (v['d_t'] + v['d_b']) if dense_share is None else dense_share
            dd, dl = d_per * sh, d_per * (1 - sh)
            n0 = P.bruggeman(Np, host, fb); ksc = ksc_scale * v['k_sc'] * (633.0 / WL)**v['p']
            per = [(n0 + 1j * ksc, dl + extra, 'lower'), (Np, dd, 'perov')]             # substrate side first
            share = ksc / np.maximum(ksc + n0.imag, 1e-30); vol = dict(L=fb * dl, D=dd)
    if film:
        lay = per; last = one                                                            # substrate (taken as the glass) / film / air
    else:
        under = [(1.9 * one, 25.0, 'par'), (1.46 * one, 25.0, 'par')] if undercoat else []   # SnO2 / SiO2 colour-suppression coat of FTO glass
        lay = under + [(n_fto * one, d_fto, 'par'), (hil, d_hil, 'par')] + per + [(n_zadn * one, d_zadn, 'par'), (1.4 * one, D_LIF, 'par')]; last = al
    n = np.array([N_GLASS * one] + [x[0] for x in lay] + [last]); d = [x[1] for x in lay]; roles = [x[2] for x in lay]
    src = dict(L=(1 + roles.index('lower')) if 'lower' in roles else None, D=1 + roles.index('perov'))
    return n, d, roles, share, src, vol, ('out' if film else 'par')


# ------------------------------------------------------------------ one emission event followed through the stack
def _fresnel_ga(s):
    ca = np.sqrt(np.clip(1 - s**2, 0, None)); cg = np.sqrt(np.clip(1 - (s / N_GLASS)**2, 0, None))
    with np.errstate(divide='ignore', invalid='ignore'):
        rs = ((N_GLASS * cg - ca) / (N_GLASS * cg + ca))**2; rp = ((cg - N_GLASS * ca) / (cg + N_GLASS * ca))**2
    return dict(TE=np.where(s < 1, rs, 1.0), TM=np.where(s < 1, rp, 1.0))


def book(r, pw, s, roles, share, sub_lost, last='par'):
    """flows per in-plane index -> integrated fates (not yet normalised) and the flow into air per s.
    last = 'out': the last medium is air (a bare film): what enters it has left the film and is counted in 'air'."""
    prop = s < N_GLASS; Rga = _fresnel_ga(s); out = dict(K=0.0, air=0.0, sub=0.0, Ap=0.0, Apar=0.0, S=0.0); air_s = 0.0
    for pol in ('TE', 'TM'):
        x = r[pol]; Ff = np.where(prop, x['F_first'], 0.0); Rst = np.where(prop, pw[pol]['R'], 0.0)
        air = Ff * (1 - Rga[pol]) / np.maximum(1 - Rga[pol] * Rst, 1e-12); back = Ff - air
        stuck = (1 - Rst) < 1e-6                                      # a stack that absorbs nothing cannot take the trapped light back
        lost = np.where(s >= 1, back, 0.0) if sub_lost else np.where(stuck & (s >= 1), back, 0.0); back = back - lost
        den = np.maximum(np.where(prop, 1 - Rst, 1.0), 1e-6)
        ab = [x['A'][i] + back * np.where(prop, pw[pol]['A'][i], 0.0) / den for i in range(len(roles))]
        lastf = x['F_last'] + back * np.where(prop, pw[pol]['T'], 0.0) / den
        par = (lastf if last == 'par' else 0.0) + sum(a for a, ro in zip(ab, roles) if ro == 'par')
        if last == 'out': air = air + lastf
        low = sum(a for a, ro in zip(ab, roles) if ro == 'lower') if 'lower' in roles else 0.0
        per = sum(a for a, ro in zip(ab, roles) if ro == 'perov')
        for k, y in (('K', x['K']), ('air', air), ('sub', lost), ('Ap', per + low * (1 - share)), ('Apar', par), ('S', low * share)):
            out[k] = out[k] + DP.integrate(y + 0.0 * s, s)
        air_s = air_s + air
    return out, air_s


def components(kind, fit=None, horiz=2 / 3, **kw):
    """Per source part ('L', 'D'), per substrate mode, per offset and wavelength: the fate fractions and the
    angular pattern (per photon, per steradian in air).  -> dict of arrays."""
    sa = np.sin(np.deg2rad(np.where(TH == 0, 0.5, TH))); res = {}
    for io, off in enumerate(OFFS):
        n, d, roles, share, src, vol, last = build(kind, fit, extra=off, **kw)
        for iw, lam in enumerate(WL):
            nn = n[:, iw]; smax = float(np.max(nn[1:-1].real)); crit = [1.0, N_GLASS] + [float(x) for x in nn[1:-1].real]
            s = DP.s_grid(smax, 900, crit); pw = DP.plane_wave(lam, nn, d, s); pwa = DP.plane_wave(lam, nn, d, sa)
            for part, j in src.items():
                if j is None: continue
                dj = d[j - 1]
                for zz in (np.arange(NZ) + 0.5) / NZ * dj:
                    r = DP.isotropic(DP.dipole(lam, nn, d, j, zz, s), horiz); ra = DP.isotropic(DP.dipole(lam, nn, d, j, zz, sa), horiz)
                    for mode, sl in (('lost', True), ('kept', False)):
                        f, _ = book(r, pw, s, roles, share[iw], sl, last); _, a = book(ra, pwa, sa, roles, share[iw], sl, 'par')
                        tot = sum(f[k] for k in KEYS); assert abs(tot - f['K']) < 1e-6 * abs(f['K']) + 1e-9, (tot, f['K'])
                        key = (part, mode)
                        if key not in res: res[key] = dict(f=np.zeros((len(OFFS), len(WL), len(KEYS))), I=np.zeros((len(OFFS), len(WL), len(TH))), F=np.zeros((len(OFFS), len(WL))))
                        res[key]['f'][io, iw] += np.array([f[k] for k in KEYS]) / f['K'] / NZ
                        res[key]['I'][io, iw] += a * np.sqrt(1 - sa**2) / (2 * np.pi * sa) / f['K'] / NZ
                        res[key]['F'][io, iw] += f['K'] / NZ                    # mean rate enhancement (Purcell), s <= largest real index
    return dict(res=res, vol=vol, has_L=src['L'] is not None)


# ------------------------------------------------------------------ evaluation
def spread(arr, sig):
    w = np.exp(-0.5 * (OFFS / sig)**2) if sig > 0 else (OFFS == 0).astype(float); w = w / w.sum()
    return np.tensordot(w, arr, axes=(0, 0))


def mix(C, x_d, sig, sub, scatter):
    """fates [wl, 5] and pattern [wl, theta] per photon, after the scattered light has been followed"""
    fD = spread(C['res'][('D', sub)]['f'], sig); ID = spread(C['res'][('D', sub)]['I'], sig)
    if not C['has_L']:
        return fD, ID, 1.0, fD[:, 4] * 0
    fL = spread(C['res'][('L', sub)]['f'], sig); IL = spread(C['res'][('L', sub)]['I'], sig)
    xd = C['vol']['D'] / (C['vol']['D'] + C['vol']['L']) if x_d is None else x_d
    f1 = xd * fD + (1 - xd) * fL; I1 = xd * ID + (1 - xd) * IL
    if scatter != 'follow':
        return f1, I1, xd, f1[:, 4]
    g = 1.0 / (1 - fL[:, 4]); f = f1 + (f1[:, 4] * g)[:, None] * fL; f[:, 4] = 0.0
    return f, I1 + (f1[:, 4] * g)[:, None] * IL, xd, f1[:, 4]


def evaluate(C, x_d=None, sig=15.0, sub='lost', scatter='follow', meas=None):
    """x_d: share of first-generation photons born in the dense part (None = by perovskite volume).
    scatter: 'follow' (scattered light re-emitted in the lower part) or 'lost'.  -> dict of results."""
    th, Im, el = meas; use = th <= 70; u60 = th <= 60
    f, I, xd, s_first = mix(C, x_d, sig, sub, scatter)
    spec = el / I[:, 0]; spec = spec / np.trapezoid(spec, WL)                      # internal spectrum (per photon)
    pat = np.trapezoid(spec[:, None] * I, WL, axis=0); patn = pat / pat[0]
    g = np.trapezoid(spec[:, None] * f, WL, axis=0)                                # spectrum-averaged fates
    air, Ap = g[0], g[2]
    ext = spec[:, None] * I; pk = [float(WL[int(np.argmax(ext[:, i]))]) for i in (0, 4, 7)]
    return dict(x_d=float(xd), sig=sig, sub=sub, scatter=scatter, I=[float(v) for v in patn], rms=float(np.sqrt(np.mean((patn - Im)[use]**2))),
                rms60=float(np.sqrt(np.mean((patn - Im)[u60]**2))), fates={k: float(v) for k, v in zip(KEYS, g)}, scattered_first=float(np.trapezoid(spec * s_first, WL)),
                eqe_pl1=float(air / (1 - Ap)), eta_for_max=float(EQE_MEAS / (air + EQE_MEAS * Ap)), eta_for_mean=float(EQE_MEAN / (air + EQE_MEAN * Ap)),
                internal_peak=float(WL[int(np.argmax(spec))]), el_peak_0_40_70=pk, internal_spectrum=[float(v) for v in spec / spec.max()])


def film_yield(C, spec, sig=15.0):
    """Bare film in an integrating sphere: escape probability per emission and the internal radiative efficiency
    that the measured external yield (PLQE_MEAS) then implies.  Whatever is not re-absorbed leaves the film;
    light trapped in the substrate leaves through its edges ('edge counted') or is fed back into the film ('kept')."""
    out = {}
    for sub in ('lost', 'kept'):
        f, _, xd, _ = mix(C, None, sig, sub, 'follow'); g = np.trapezoid(spec[:, None] * f, WL, axis=0) / np.trapezoid(spec, WL)
        p = g[0] + g[1]; out['edge counted' if sub == 'lost' else 'edge re-absorbed'] = dict(escape=float(p), reabsorbed=float(g[2]), eta=float(PLQE_MEAS / (p + PLQE_MEAS * (1 - p))))
    return out


SEMG = dict(d_hil=50.0, d_per=254.0, d_zadn=42.0, d_fto=367.0)        # thicknesses read off the device cross-section (SI Fig. 5)
STRUCTS = [                                   # label, kind, fit variant, options
    ('one layer, paper optics',              'paper', None, dict(n_fto=1.9 + 0j)),
    ('one layer, dense n,k',                 'dense', 'base', {}),
    ('two layers, no scattering',            'SE', 'base', dict(ksc_scale=0.0)),
    ('two layers',                           'SE', 'base', {}),
    ('SEM geometry',                         'SE', 'base', SEMG),
    ('SEM geometry, FTO undercoat',          'SE', 'base', dict(SEMG, undercoat=True)),
    ('SEM geometry, no scattering',          'SE', 'base', dict(SEMG, ksc_scale=0.0)),
    ('SEM geometry, scattering x 0.5',       'SE', 'base', dict(SEMG, ksc_scale=0.5)),
    ('SEM geometry, lower fraction 0.35',    'SE', 'base', dict(SEMG, f_b=0.35)),
    ('SEM geometry, lower fraction 0.70',    'SE', 'base', dict(SEMG, f_b=0.70)),
    ('SEM geometry, dense share 0.20',       'SE', 'base', dict(SEMG, dense_share=0.20)),
    ('SEM geometry, dense share 0.32',       'SE', 'base', dict(SEMG, dense_share=0.32)),
    ('SEM geometry, dense share 0.53',       'SE', 'base', dict(SEMG, dense_share=0.53)),
    ('SEM geometry, HIL 35 nm',              'SE', 'base', dict(SEMG, d_hil=35.0)),
    ('SEM geometry, HIL 65 nm',              'SE', 'base', dict(SEMG, d_hil=65.0)),
    ('SEM geometry, HIL 75 nm',              'SE', 'base', dict(SEMG, d_hil=75.0)),
    ('SEM geometry, perovskite 270 nm',      'SE', 'base', dict(SEMG, d_per=270.0)),
    ('SEM geometry, fit "organic matrix"',   'SE', 'organic matrix', SEMG),
    ('SEM geometry, stitched band edge',     'SE', 'base', dict(SEMG, table='stitched')),
    ('SEM geometry, k x 1.5',                'SE', 'base', dict(SEMG, k_edge_scale=1.5)),
    ('SEM geometry, horizontal dipoles 0.75', 'SE', 'base', dict(SEMG, horiz=0.75)),
    ('SEM geometry, FTO k 0.02',             'SE', 'base', dict(SEMG, n_fto=1.9 + 0.02j)),
    ('SEM geometry, ZADN n 1.85',            'SE', 'base', dict(SEMG, n_zadn=1.85)),
    ('SEM geometry, dense n x 1.05',         'SE', 'base', dict(SEMG, n_scale=1.05)),
]
FILMS = [                                     # bare films for the photoluminescence yield: label, the device structure it belongs to, kind, fit, options
    ('film, paper optics',                   'one layer, paper optics', 'paper', None, dict(film=True)),
    ('film, two layers, no scattering',      'two layers, no scattering', 'SE', 'base', dict(film=True, ksc_scale=0.0)),
    ('film, two layers',                     'SEM geometry', 'SE', 'base', dict(film=True, d_per=254.0)),
    ('film, two layers, scattering x 0.5',   'SEM geometry, scattering x 0.5', 'SE', 'base', dict(film=True, d_per=254.0, ksc_scale=0.5)),
    ('film, fit "organic matrix"',           'SEM geometry, fit "organic matrix"', 'SE', 'organic matrix', dict(film=True, d_per=254.0)),
]


def job(a):
    label, kind, var, kw = a
    return label, components(kind, DC.fit_of(var) if var else None, **kw)


def main():
    from multiprocessing import Pool
    meas = DC.measured(); th, Im, el = meas; el = np.interp(WL, DC.WL, el); meas = (th, Im, el)
    only = sys.argv[1:]
    todo = [s for s in STRUCTS if not only or s[0] in only]; films = [f for f in FILMS if not only]
    with Pool(8) as pool: comps = dict(pool.map(job, todo + [(f[0], f[2], f[3], f[4]) for f in films], chunksize=1))
    out = {}
    print('measured: ' + ' '.join('%.3f' % v for v in Im))
    print('%-38s %-5s %-7s %5s | %5s %5s | %5s %5s %5s %5s | %5s %6s %6s | %s' % ('structure', 'sub', 'scatter', 'x_d', 'rms70', 'rms60', 'air', 'sub', 'Ap', 'Apar', 'EQE1', 'eta25', 'eta29', 'I(20) I(40) I(60) I(70)'))
    for label, kind, var, kw in todo:
        C = comps[label]; rows = []
        modes = [('lost', 'follow'), ('kept', 'follow')] + ([('lost', 'lost')] if C['has_L'] else [])
        for sub, sc in modes:
            for xd in [None] + ([0.0, 0.2, 0.4, 0.6, 0.8, 1.0] if C['has_L'] else []):
                for sig in (0.0, 15.0, 25.0):
                    rows.append(evaluate(C, xd, sig, sub, sc, meas) | dict(x_d_mode='volume' if xd is None else 'fixed'))
        out[label] = rows
        for sub, sc in modes:
            sel = [r for r in rows if r['sub'] == sub and r['scatter'] == sc and r['sig'] == 15.0]
            show = [r for r in sel if r['x_d_mode'] == 'volume'] + ([min((r for r in sel if r['x_d_mode'] == 'fixed'), key=lambda r: r['rms'])] if (C['has_L'] and sub == 'lost') else [])
            for r in show:
                f = r['fates']
                print('%-38s %-5s %-7s %5.2f%s| %5.3f %5.3f | %5.3f %5.3f %5.3f %5.3f | %5.3f %6.3f %6.3f | %s' % (label[:38], sub, sc, r['x_d'], 'v' if r['x_d_mode'] == 'volume' else '*', r['rms'], r['rms60'],
                      f['air'], f['sub'], f['Ap'], f['Apar'], r['eqe_pl1'], r['eta_for_mean'], r['eta_for_max'], ' '.join('%.3f' % r['I'][i] for i in (2, 4, 6, 7))))
    # ---- photoluminescence yield of the bare film -> internal radiative efficiency -> EQE the device model then predicts
    fy = {}
    if films: print('\n%-36s %-17s %7s %7s | device structure and the EQE it then predicts (sub lost / kept)' % ('bare film (external PLQE %.0f %%)' % (100 * PLQE_MEAS), 'substrate light', 'escape', 'eta'))
    for label, dev, kind, var, kw in films:
        if dev not in out: continue
        ref = [r for r in out[dev] if r['sub'] == 'lost' and r['scatter'] == 'follow' and r['sig'] == 15.0 and r['x_d_mode'] == 'volume'][0]
        spec = np.array(ref['internal_spectrum']); y = film_yield(comps[label], spec); fy[label] = dict(device=dev, **y)
        for k, v in y.items():
            pred = []
            for sub in ('lost', 'kept'):
                r = [r for r in out[dev] if r['sub'] == sub and r['scatter'] == 'follow' and r['sig'] == 15.0 and r['x_d_mode'] == 'volume'][0]; f = r['fates']
                pred.append(v['eta'] * f['air'] / (1 - v['eta'] * f['Ap']))
            v['eqe_predicted'] = dict(sub_lost=float(pred[0]), sub_kept=float(pred[1]))
            print('%-36s %-17s %7.3f %7.3f | %-36s %.3f / %.3f' % (label, k, v['escape'], v['eta'], dev, *pred))
    prev = json.load(open(W + 'device_model.json')) if (only and _os.path.exists(W + 'device_model.json')) else {}
    res = prev.get('results', {}); res.update(out)
    json.dump(dict(theta=TH.tolist(), measured=[float(v) for v in Im], eqe_measured=EQE_MEAS, eqe_mean=EQE_MEAN, plqe_measured=PLQE_MEAS, results=res, film_yield=fy or prev.get('film_yield', {})),
              open(W + 'device_model.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

"""Does the ellipsometric n,k reproduce the angular emission of the green PeLED?

Device (Nature 611, 688, Methods): glass / FTO 350 / BufHIL 75 / perovskite 270 /
ZADN 45 / LiF 1.2 / Al 100.  Measured: Extended Data Fig. 7g (intensity against
angle, slightly narrower than Lambertian) and 7c (EL spectrum, 540 nm, 22 nm
FWHM), both digitised.

Scenarios for the perovskite layer
  paper       n = 2.3 constant, k from the absorbance (SI Fig. 6a), one 270 nm layer, BufHIL n = 1.5
              - the optics the paper's own simulation assumed
  dense       the ellipsometric dense-perovskite n,k as one 270 nm layer
  SE          the ellipsometric structure (gp_perov): the lower, low-index, scattering part on the
              HIL side and the dense part on the ETL side, scaled to 270 nm in total
BufHIL: the ellipsometric GraHIL n,k at the paper's 75 nm (and at the 108 nm of the Si sample).
Emitters: uniform in the perovskite volume (a sublayer weighs thickness x perovskite fraction).
Each scenario is averaged over a lateral spread of the perovskite thickness (different places emit
incoherently).  The internal emission spectrum is not taken from the paper (theirs was back-
calculated with their optics): it is solved so that the model's normal-incidence output is the
measured EL spectrum.  Light scattered in the lower layer counts as lost to the coherent pattern;
eta is the share of Lambertian light that would best complete a scenario (0 = none needed).
What the scattered light itself would look like, if it leaves at all, is the row
'fully scattered': isotropic radiance in the glass, refracted into air."""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gp_model as P, gp_device as Dv, gp_perov as PV

W = _os.path.join(ELLIPS_OUT, '')
WL = np.arange(505.0, 586.0, 1.0)
N_GLASS, N_FTO, N_ZADN, N_LIF = 1.5, 1.9 + 0.005j, 2.0, 1.4
D_FTO, D_ZADN, D_LIF, D_PER = 350.0, 45.0, 1.2, 270.0
K_PAPER = ([505, 510, 515, 520, 525, 530, 535, 540, 545, 550, 555, 560, 565, 570], [0.110, 0.113, 0.119, 0.125, 0.132, 0.125, 0.082, 0.035, 0.014, 0.006, 0.003, 0.0015, 0.0004, 0.0])
ONE = np.ones(len(WL), complex)


def measured():
    j = json.load(open(W + 'ed7g_digitized.json')); e = json.load(open(W + 'ed7c_digitized.json'))
    return np.array(j['theta'], float), np.array(j['I']), np.interp(WL, e['wl'], e['el'])


def diffuse(th_deg, n_glass=N_GLASS):
    """Light that was scattered on its way out: isotropic radiance in the glass, then glass -> air.
    I(theta) = cos(theta) T(theta) / T(0) - a Lambertian source times the Fresnel transmission."""
    ta = np.deg2rad(np.asarray(th_deg, float)); tg = np.arcsin(np.sin(ta) / n_glass); cg, ca = np.cos(tg), np.cos(ta)
    T = 1 - 0.5 * (((n_glass * cg - ca) / (n_glass * cg + ca))**2 + ((cg - n_glass * ca) / (cg + n_glass * ca))**2)
    return ca * T / T[0]


def fit_of(name='base'):
    return [r for r in json.load(open(W + 'perov_fits.json'))['variants'] if r['variant'] == name][0]


def perov(kind, fit, mod):
    """perovskite sublayers from the HIL side to the ETL side: [(N, d, emission weight per nm)], total D_PER"""
    if kind == 'paper':
        return [(2.3 + 1j * np.interp(WL, *K_PAPER), D_PER, 1.0)]
    Np = PV.nk_of(fit, WL) * mod.get('n_scale', 1.0)
    if kind == 'dense':
        return [(Np, D_PER, 1.0)]
    v = dict(fit['struct']); v['k_sc'] = v['k_sc'] * mod.get('ksc_scale', 1.0)
    if 'f_b' in mod: v['f_b'] = mod['f_b']
    if 'd_t' in mod: v['d_b'] = v['d_t'] + v['d_b'] - mod['d_t']; v['d_t'] = mod['d_t']
    pl = PV.perov_layers(WL, v, Np); tot = sum(x[1] for x in pl)
    zone = mod.get('zone', 'uniform')
    wt = lambda f, top: f if zone == 'uniform' else (1.0 if (zone == 'dense only') == top else 0.0)
    return [(N, d * D_PER / tot, wt(f, i == 0)) for i, (N, d, f) in enumerate(pl)][::-1]


def pattern(kind, fit=None, d_hil=75.0, sig=15.0, n_zadn=N_ZADN, spec=None, **mod):
    """-> I(theta) (not normalised), outcoupled spectrum at 0 deg per unit internal spectrum"""
    th, _, _ = measured(); spec = np.ones(len(WL)) if spec is None else spec
    Ng = 1.5 * ONE if kind == 'paper' else PV.grahil_N(WL)[0]
    al = Dv.n_al(WL); x, w = P._GH[5 if sig > 0 else 1]
    I = 0.0; ext0 = 0.0
    for xi, wi in zip(x, w):
        pl = perov(kind, fit, mod); pl = [(N, d + (sig * xi if i == 0 else 0.0), e) for i, (N, d, e) in enumerate(pl)]   # the spread sits on the HIL side
        st = [N_GLASS * ONE, N_FTO * ONE, Ng] + [p[0] for p in pl] + [n_zadn * ONE, N_LIF * ONE, al]
        ds = [D_FTO, d_hil] + [p[1] for p in pl] + [D_ZADN, D_LIF]
        for i, (N, d, e) in enumerate(pl):
            if e <= 0 or d <= 0: continue
            Ia, ext = Dv.angular(WL, spec, st, ds, 3 + i, th, nz=max(6, int(d / 6)))
            I = I + wi * e * d * Ia; ext0 = ext0 + wi * e * d * ext[0]
    return I, ext0


def run(label, kind, fit=None, **kw):
    th, Im, el = measured(); lam = np.cos(np.deg2rad(th)); use = th <= 70
    _, out0 = pattern(kind, fit, **kw)                     # outcoupling at 0 deg for a flat internal spectrum
    spec = el / out0; spec = spec / spec.max()             # internal spectrum that gives the measured EL at 0 deg
    I, _ = pattern(kind, fit, spec=spec, **kw); I = I / I[0]
    a = I[use] - lam[use]; eta = float(np.clip(np.sum((I[use] - Im[use]) * a) / np.sum(a * a), 0.0, 1.0))
    Imix = (1 - eta) * I + eta * lam; pk = float(WL[int(np.argmax(spec))])
    return dict(label=label, I=[float(v) for v in I], rms=float(np.sqrt(np.mean((I - Im)[use]**2))), eta=eta,
                rms_mix=float(np.sqrt(np.mean((Imix - Im)[use]**2))), internal_peak=pk)


def main():
    th, Im, el = measured(); lam = np.cos(np.deg2rad(th)); use = th <= 70
    base = fit_of('base'); rows = []
    for sig in (0.0, 15.0, 25.0):
        rows.append(run('paper optics                 spread %2.0f' % sig, 'paper', sig=sig))
        rows.append(run('dense n,k, one 270 nm layer  spread %2.0f' % sig, 'dense', base, sig=sig))
        rows.append(run('SE structure (base fit)      spread %2.0f' % sig, 'SE', base, sig=sig))
    rows.append(None)
    for name in ('GraHIL free', 'surface spread', 'organic matrix', 'graded lower', 'soft boundary', 'graded+spread', 'all'):
        rows.append(run('SE structure, fit variant "%s"' % name, 'SE', fit_of(name)))
    rows.append(None)
    for lab, kw in (('BufHIL 108 nm (as on the Si sample)', dict(d_hil=108.0)), ('ZADN n 1.85', dict(n_zadn=1.85)),
                    ('emission in the dense part only', dict(zone='dense only')), ('emission in the lower part only', dict(zone='lower only')),
                    ('scattering extinction off', dict(ksc_scale=0.0)), ('scattering extinction x 0.5', dict(ksc_scale=0.5)), ('scattering extinction x 1.5', dict(ksc_scale=1.5)),
                    ('dense n x 0.95', dict(n_scale=0.95)), ('dense n x 1.05', dict(n_scale=1.05)),
                    ('lower-layer perovskite fraction 0.35', dict(f_b=0.35)), ('lower-layer perovskite fraction 0.70', dict(f_b=0.70)),
                    ('lower-layer perovskite fraction 1.0 (no low-index part)', dict(f_b=1.0)),
                    ('dense part 90 nm of 284', dict(d_t=90.0)), ('dense part 150 nm of 284', dict(d_t=150.0))):
        rows.append(run('SE base, ' + lab, 'SE', base, **kw))
    dif = diffuse(th); rms = lambda I: float(np.sqrt(np.mean((I - Im)[use]**2)))
    json.dump(dict(theta=th.tolist(), measured=Im.tolist(), lambertian=dict(I=lam.tolist(), rms=rms(lam)), diffuse=dict(I=dif.tolist(), rms=rms(dif)),
                   rows=[r for r in rows if r]), open(W + 'device_check.json', 'w'), indent=1)
    print('%-58s ' % 'theta (deg)' + ' '.join('%5d' % t for t in th) + '  rms(0-70)   eta -> rms   internal peak')
    print('%-58s ' % 'measured (ED Fig. 7g)' + ' '.join('%5.3f' % v for v in Im)); print('%-58s ' % 'Lambertian' + ' '.join('%5.3f' % v for v in lam) + '    %.3f' % rms(lam))
    print('%-58s ' % 'fully scattered (isotropic in the glass)' + ' '.join('%5.3f' % v for v in dif) + '    %.3f' % rms(dif))
    for r in rows:
        if r is None: print(); continue
        print('%-58s ' % r['label'] + ' '.join('%5.3f' % v for v in r['I']) + '    %.3f     %.2f -> %.3f     %.0f nm' % (r['rms'], r['eta'], r['rms_mix'], r['internal_peak']))


if __name__ == '__main__':
    main()

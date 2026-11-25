"""Is the small surface roughness of the layer2 fit an artefact?  The fit of gp_perov repeated with the
roughness imposed instead of fitted.

Two kinds of roughness, which the optics treat differently:
  d_r     effective-medium layer (dense perovskite / void 50:50): roughness whose lateral scale is far
          below the wavelength.  For that kind, thickness ~ 1.5 x rms + 0.4 nm.
  sig_s   Gaussian spread of the surface height with a lateral scale too large for an effective medium
          (film undulation, grains of 100 nm and more): amplitudes are averaged coherently.  To first
          order it leaves the ratio r_p / r_s alone, which is why ellipsometry is nearly blind to it and
          only the fringe contrast and the intensity see it.
The device cross-section (SI Fig. 5 of Nature 611, 688) gives a perovskite top-surface rms of ~13.5 nm,
of which 2-3 nm at lateral periods below 300 nm and ~11 nm at 0.3-1 um; the top-view SEM of the
core/shell film gives close-packed ~70 nm cobbles with shallow grooves (rms 2-6 nm).

Each row refits everything else (thicknesses, lower layer, scattering law, the dense n,k basis),
starting from the base fit.  -> rough_scan.json"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gp_perov as PV

W = _os.path.join(ELLIPS_OUT, '')
SCAN = {'d_r 3': dict(d_r=3.0), 'd_r 6': dict(d_r=6.0), 'd_r 10': dict(d_r=10.0), 'd_r 15': dict(d_r=15.0), 'd_r 20': dict(d_r=20.0), 'd_r 30': dict(d_r=30.0),
        'sig_s 5': dict(sig_s=5.0), 'sig_s 10': dict(sig_s=10.0), 'sig_s 13.5': dict(sig_s=13.5), 'sig_s 20': dict(sig_s=20.0),
        'd_r 6, sig_s 13.5': dict(d_r=6.0, sig_s=13.5), 'd_r 10, sig_s 13.5': dict(d_r=10.0, sig_s=13.5)}
for k, f in SCAN.items(): PV.VARIANTS[k] = ([], dict(fixed=f))


def main():
    from multiprocessing import Pool
    old = json.load(open(W + 'perov_fits.json')); base = [r for r in old['variants'] if r['variant'] == 'base'][0]
    va = dict(old['stage_a'], **base['struct']); xb = np.array(base['x'][len(base['names']):])
    only = sys.argv[1:]; todo = [k for k in SCAN if not only or k in only]
    with Pool(4) as pool: res = pool.map(PV.run_variant, [(k, va, xb) for k in todo], chunksize=1)
    prev = json.load(open(W + 'rough_scan.json'))['scan'] if (only and _os.path.exists(W + 'rough_scan.json')) else []
    res = [r for r in prev if r['variant'] not in todo] + res
    json.dump(dict(base=base, scan=res), open(W + 'rough_scan.json', 'w'), indent=1)
    q = np.array([300.0, 400.0, 500.0, 540.0, 633.0])
    print('%-20s %6s %6s | %-27s | %5s %5s %5s %5s %5s %5s | n,k at 300 / 400 / 500 / 540 ; n(633)' % ('roughness', 'MSE', 'rms R', '205-300 300-450 450-570 570-1000', 'd_r', 'sig', 'd_t', 'd_b', 'f_b', 'k_sc'))
    for r in [dict(base, variant='fitted (base)')] + sorted(res, key=lambda r: list(SCAN).index(r['variant'])):
        N = PV.nk_of(r, q); v = r['struct']
        print('%-20s %6.2f %6.1f | %s | %5.1f %5.1f %5.0f %5.0f %5.2f %5.3f | ' % (r['variant'], r['mse_ncs'], r['rms_R'], '  '.join('%5.1f' % b for b in r['bands'].values()),
              v['d_r'], v.get('sig_s', 0.0), v['d_t'], v['d_b'], v['f_b'], v['k_sc']) + '  '.join('%.3f %.3f' % (n.real, n.imag) for n in N[:4]) + '  ; %.3f' % N[4].real)


if __name__ == '__main__':
    main()

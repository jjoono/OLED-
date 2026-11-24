"""Can each description of the perovskite reproduce the angular emission with a geometry that the
device cross-section (SI Fig. 5 of Nature 611, 688) allows?

The thicknesses in the Methods are nominal; the micrograph gives FTO 367 +- 15 (342-396), BufHIL
50 +- 11 (34-67), perovskite 254 +- 14 (231-273), ZADN 42 +- 5 (35-51) nm, plus a ~47 nm coat under
the FTO.  Within those ranges (the BufHIL up to its nominal 75 nm) a random search looks for the best
angular pattern of three families:

  one layer         the dense n,k as a single layer
  two, no scatter   dense on a low-index lower part, no scattering
  two, followed     the same with the scattering of the lower part followed (gp_device_model)

free: the four thicknesses; for the two-layer families the dense share (0.15-0.53) and x_d, the share
of the emission in the dense part.  The search runs at reduced resolution; the best few of each family
are re-evaluated at full resolution.  With eight angular points and five or six free numbers this is
a test of reachability, not a determination of the geometry.  -> device_fit.json"""
import os as _os
ELLIPS_DATA = _os.environ.get('ELLIPS_DATA', '.')
ELLIPS_OUT  = _os.environ.get('ELLIPS_OUT', '.')
import sys, json, numpy as np
sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import gp_device_model as DM, gp_device_check as DC

W = _os.path.join(ELLIPS_OUT, '')
BOUNDS = dict(d_fto=(342.0, 396.0), d_hil=(34.0, 75.0), d_per=(231.0, 273.0), d_zadn=(35.0, 51.0), dense_share=(0.15, 0.53))
FAM = {'one layer': ('dense', {}), 'two, no scatter': ('SE', dict(ksc_scale=0.0)), 'two, followed': ('SE', {})}
N_SEARCH, N_TOP = 140, 4
XD = [None, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
FULL = dict(WL=DM.WL.copy(), OFFS=DM.OFFS.copy(), NZ=DM.NZ)
FAST = dict(WL=np.arange(505.0, 585.1, 5.0), OFFS=np.array([-30.0, -15.0, 0.0, 15.0, 30.0]), NZ=4)


def run(a):
    fam, geo, fast = a
    for k, v in (FAST if fast else FULL).items(): setattr(DM, k, v)
    kind, kw = FAM[fam]; g = {k: v for k, v in geo.items() if not (k == 'dense_share' and kind == 'dense')}
    th, Im, el = DC.measured(); el = np.interp(DM.WL, DC.WL, el)
    C = DM.components(kind, DC.fit_of('base'), undercoat=True, **kw, **g)
    rows = [DM.evaluate(C, xd, 15.0, 'lost', 'follow', (th, Im, el)) for xd in (XD if C['has_L'] else [None])]
    best = min(rows, key=lambda r: r['rms']); vol = rows[0]
    kept = DM.evaluate(C, best['x_d'], 15.0, 'kept', 'follow', (th, Im, el))
    return dict(family=fam, geo=geo, best=best, volume=vol, eqe1_kept=kept['eqe_pl1'])


def main():
    from multiprocessing import Pool
    rng = np.random.default_rng(1); jobs = []
    for fam in FAM:
        for i in range(N_SEARCH):
            jobs.append((fam, {k: float(lo + rng.random() * (hi - lo)) for k, (lo, hi) in BOUNDS.items()}, True))
    with Pool(8) as pool:
        res = pool.map(run, jobs, chunksize=2)
        top = []
        for fam in FAM:
            top += sorted((r for r in res if r['family'] == fam), key=lambda r: r['best']['rms'])[:N_TOP]
        fin = pool.map(run, [(r['family'], r['geo'], False) for r in top], chunksize=1)
    out = {}
    print('measured            ' + ' '.join('%.3f' % v for v in DC.measured()[1]))
    print('%-16s %5s %5s %5s %5s %6s %5s | %6s %6s | %5s %5s %5s | %s' % ('family', 'FTO', 'HIL', 'per', 'ZADN', 'dense', 'x_d', 'rms70', 'rms60', 'EQE1', 'kept', 'air', 'I(10..80)'))
    for fam in FAM:
        rs = sorted((r for r in fin if r['family'] == fam), key=lambda r: r['best']['rms']); out[fam] = rs
        allr = sorted(r['best']['rms'] for r in res if r['family'] == fam)
        for r in rs[:2]:
            g, b = r['geo'], r['best']
            print('%-16s %5.0f %5.0f %5.0f %5.0f %6.2f %5.2f | %6.3f %6.3f | %5.3f %5.3f %5.3f | %s' % (fam, g['d_fto'], g['d_hil'], g['d_per'], g['d_zadn'], g['dense_share'], b['x_d'], b['rms'], b['rms60'],
                  b['eqe_pl1'], r['eqe1_kept'], b['fates']['air'], ' '.join('%.3f' % v for v in b['I'][1:])))
        print('%-16s search (reduced resolution): best %.3f, 10th percentile %.3f, median %.3f, worst %.3f' % ('', allr[0], allr[len(allr) // 10], allr[len(allr) // 2], allr[-1]))
    json.dump(dict(bounds=BOUNDS, n_search=N_SEARCH, families=out, search={fam: [dict(geo=r['geo'], rms=r['best']['rms'], x_d=r['best']['x_d']) for r in res if r['family'] == fam] for fam in FAM}),
              open(W + 'device_fit.json', 'w'), indent=1)


if __name__ == '__main__':
    main()

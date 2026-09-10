"""Geometry scan (d, roughness, angle offset) scored by the model-free chain
residual over 400-1300 nm - the range where a homogeneous layer is valid."""
import numpy as np, ce_fit as cf, ce_chain as cc, sys, json

LO, HI = 400.0, 1300.0

def run(sheet, ds, rgs, dths, stride=10):
    wl, Pm, Dm = cf.load(sheet)
    m = (wl >= LO) & (wl <= HI); wl, Pm, Dm = wl[m], Pm[m], Dm[m]
    ox, si = cf._mat('NTVE_JAW', wl), cf._mat('SI_JAW', wl)
    best = None; tab = []
    for d in ds:
        for rg in rgs:
            for dth in dths:
                c = cc.chain(wl, Pm, Dm, ox, si, d, rg, dth, stride)
                s = np.sqrt(np.mean(c[:, 3]**2))
                tab.append((s, d, rg, dth))
                if best is None or s < best[0]: best = (s, d, rg, dth)
    tab.sort()
    return best, tab

if __name__ == '__main__':
    sheet = sys.argv[1] if len(sys.argv) > 1 else '#1'
    b, tab = run(sheet, [47., 49., 51., 53., 55.], [0.0, 1.5, 3.0], [0.0, 0.4, 0.8])
    print('%s  top 8 of %d:' % (sheet, len(tab)))
    for s, d, rg, dth in tab[:8]:
        print('   rms=%.5f   d=%.1f  rough=%.1f  dth=%+.2f' % (s, d, rg, dth))
    json.dump(dict(d=b[1], rough=b[2], dth=b[3], rms=b[0]),
              open('ce_geo_%s.json' % sheet.strip('#'), 'w'))

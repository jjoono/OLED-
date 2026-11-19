"""Band summary (700/750/800/850 nm): reported record EQE vs. the EQE reachable with the low-loss multi-pass structure.
For each emitter the effective emission yield q_eff (PLQY x charge balance x spin factor, lumped) is calibrated so that the
optimised flat-glass ITO device reproduces the reported planar EQE: EQE(q) = q N / (1 - q + q F), N = F eta_sub eta_ext.
The same q_eff is then applied to the best low-loss structure (its F, eta_sub, eta_ext).  q_eff from an optimised flat
stack is a lower bound for the real material (the reported devices are not necessarily optimal), so the projected EQE is
conservative in that respect; it is optimistic in assuming a laterally infinite substrate and no electrical penalty."""
import json, glob, os
REPORTED = {   # emitter: (band, peak nm, reported max EQE %, source)
    'DTN': (700, 674, 20.1, 'Adv. Funct. Mater. 2025, e202529822 (Lumtec LT-N7047)'),
    'DCNSPTPA': (700, 696, 28.2, 'Adv. Funct. Mater. 2023, 2304398 (Lumtec LT-N7039)'),
    'BF2': (700, 721, 10.0, 'Nat. Photon. 12, 98 (2018) (Lumtec LT-N7003)'),
    'Pt-fprpz2': (750, 740, 24.0, 'Nat. Photon. 11, 63 (2017)'),
    'PyCNTPA-12': (750, 747, 16.3, 'TADF, CBP 12 wt%'),
    'PtTPBP': (750, 765, 8.5, 'Borek et al. 2007 (Lumtec LT-N734)'),
    'DCNSPTPA-neat': (800, 800, 2.61, 'Adv. Funct. Mater. 2023, 2304398, nondoped (Lumtec LT-N7039)'),
    'Pt794': (800, 794, 10.6, 'Adv. Funct. Mater. 2020, 2002173'),
    'PtazaTBP': (850, 848, 2.8, 'Appl. Phys. Lett. 109, 233302 (2016)'),
}
def load(em):
    rs = [json.load(open(f)) for f in glob.glob('results_emitter/%s_*.json' % em)]
    flat = [r for r in rs if r['extractor'] == 'flat']; best = [r for r in rs if r['extractor'] != 'flat']
    return (flat[0] if flat else None), best, rs
def eqe(r, q):
    N = r['F'] * r['eta_sub'] * r['eta_ext']; return q * N / (1 - q + q * r['F'])
rows = []
for em, (band, pk, rep, src) in REPORTED.items():
    flat, best, rs = load(em)
    if flat is None or not best: continue
    N = flat['F'] * flat['eta_sub'] * flat['eta_ext']; e = rep / 100
    q = e / (N - e * (flat['F'] - 1)); capped = q > 1; q = min(q, 1.0)
    best = max(best, key=lambda r: eqe(r, q))      # best structure at the calibrated q (Purcell matters at low q)
    rows.append(dict(band=band, emitter=em, peak=pk, reported=rep, flat_q1=100 * eqe(flat, 1.0), q_eff=q, capped=capped,
                     best=100 * eqe(best, q), best_q1=100 * eqe(best, 1.0), gain=eqe(best, q) / e,
                     struct='%s n%.1f DBR%d' % (best['TE'], best['n_sub'], best['ndbr']), F=best['F'],
                     eta_sub=best['eta_sub'], eta_ext=best['eta_ext'], src=src))
rows.sort(key=lambda r: (r['band'], r['peak']))
hdr = '| band | emitter | peak (nm) | reported EQE (%) | q_eff | flat model, q=1 (%) | low-loss structure | F | eta_sub | eta_ext | projected EQE (%) | gain | EQE at q=1 (%) |\n|' + '---|' * 13
lines = [hdr] + ['| %d | %s | %d | %.1f | %.2f%s | %.1f | %s | %.2f | %.3f | %.3f | **%.1f** | %.1fx | %.1f |' % (
    r['band'], r['emitter'], r['peak'], r['reported'], r['q_eff'], ' (cap)' if r['capped'] else '', r['flat_q1'], r['struct'],
    r['F'], r['eta_sub'], r['eta_ext'], r['best'], r['gain'], r['best_q1']) for r in rows]
open('band_summary.md', 'w').write(__doc__ + '\n\n' + '\n'.join(lines) + '\n\nSources:\n' + '\n'.join('- %s: %s' % (r['emitter'], r['src']) for r in rows) + '\n')
print('\n'.join(lines))

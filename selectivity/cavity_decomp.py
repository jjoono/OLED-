"""Lens versus cavity in the freeform-vs-hemisphere benchmark, in the full device model (birefringent CPS + TMM,
dispersive constants, spectrally weighted with the Ir(dmppy-ph)2tmd spectrum, 400-700 nm every 10 nm).
For each band/total objective j:  EQE_model = eta_rad * sum_lam s(lam) eta_sub(lam; cavity) eta_sta(lam; lens, cavity),
band EQE = same with the air power in the band.  The freeform/hemisphere ratio is split into a lens part (both lenses
at the hemisphere's cavity) and a cavity part (hemisphere lens at the freeform's cavity).  Output: data/cavity_decomp.json"""
import json, numpy as np, sel, cps2
MB = np.load('../multipass_reset/bsdfs.npz')
J = np.arange(0, 301, 10); LAMs = 400.0 + J
s = sel.D['s_Ir(dmppy-ph)2tmd'][J] * LAMs; s /= s.sum()
BANDS = {'0-20': (0, 20), '20-40': (20, 40), '40-60': (40, 60), '60-80': (60, 80), 'total': (0, 90)}
CACHE = {}


def device(etl, htl):
    key = (round(etl, 3), round(htl, 3))
    if key in CACHE:
        return CACHE[key]
    org = lambda m, d, j: (sel.nk(m + '_o', j).real, sel.nk(m + '_e', j).real, d)
    met = lambda m, d, j: (sel.nk(m, j), sel.nk(m, j), d)
    es, P, R = [], [], []
    for j, lam in zip(J, LAMs):
        S = cps2.Stack(float(lam), (sel.nk('EML_o', j).real, sel.nk('EML_e', j).real), 25.0, 12.5,
                       above=[org('B3', etl, j), met('Al', 100.0, j)], below=[org('TCTA', 10.0, j), org('TAPC', htl, j), met('ITO', 150.0, j)],
                       n_sub=1.51, h=0.865)
        r = cps2.solve_pol(S, npts=4000); es.append(r['air'] + r['sub'])
        th_f = np.linspace(0, 90, 1801)[:-1] + 0.025
        pf = np.clip(np.real(cps2.sub_angular(S, th_f)), 0, None); p = pf.reshape(90, 20).sum(1); P.append(p / p.sum())
        R.append(np.clip(np.real(cps2.stack_reflectance(S, np.arange(90) + 0.5)), 0, 1))
    CACHE[key] = (np.array(es), np.array(P), np.array(R)); return CACHE[key]


def eqe(lens, cav):
    BT, BR = MB[lens + '_BT'], MB[lens + '_BR'].T; es, P, R = device(*cav)
    out = {k: 0.0 for k in BANDS}
    for i in range(len(J)):
        air = np.linalg.solve(np.eye(90) - R[i][:, None] * BR, P[i]) @ BT
        for k, (a, b) in BANDS.items():
            out[k] += 0.98 * s[i] * es[i] * air[a:b].sum()
    return out


if __name__ == '__main__':
    cav = {k: (float(MB[k + '_x'][10]), float(MB[k + '_x'][11])) for k in ('hemi_20-40', 'hemi_40-60', 'hemi_60-80', 'hemi_total',
                                                                          'free_20-40', 'free_40-60', 'free_60-80', 'free_total')}
    res = {}
    for j in ('20-40', '40-60', '60-80', 'total'):
        h, f = f'hemi_{j}', f'free_{j}'
        own_h = eqe(h, cav[h])[j]; own_f = eqe(f, cav[f])[j]; lens = eqe(f, cav[h])[j]; cavity = eqe(h, cav[f])[j]
        res[j] = dict(cav_hemi=cav[h], cav_free=cav[f], hemi_own=own_h, free_own=own_f, free_at_hemi_cavity=lens, hemi_at_free_cavity=cavity,
                      ratio_own=own_f / own_h, ratio_lens=lens / own_h, ratio_cavity=cavity / own_h)
        print(j, 'cavities H %s F %s | ratio own %.3f  lens-only %.3f  cavity-only %.3f' % (np.round(cav[h], 1), np.round(cav[f], 1), own_f / own_h, lens / own_h, cavity / own_h), flush=True)
    json.dump(res, open('data/cavity_decomp.json', 'w'), indent=1)

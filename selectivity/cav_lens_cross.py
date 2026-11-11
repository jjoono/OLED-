"""Cross-check: 40-60 deg band EQE of the hemispherical and freeform 40-60 lenses at two cavities, at 550 nm and
spectrally weighted.  Output: data/cav_lens_cross.json"""
import json, numpy as np, cavity_decomp as C
out = {}
for cav in ((61.29770702992403, 97.51159318071265), (74.5974985592011, 81.66434821823951), (68.52461063471266, 114.37743028459539)):
    es, P, R = C.device(*cav)
    for lens in ('hemi_40-60', 'free_40-60', 'hemi_total'):
        BT, BR = C.MB[lens + '_BT'], C.MB[lens + '_BR'].T
        b = np.zeros(len(C.J)); t = np.zeros(len(C.J))
        for i in range(len(C.J)):
            air = np.linalg.solve(np.eye(90) - R[i][:, None] * BR, P[i]) @ BT
            b[i] = 0.98 * es[i] * air[40:60].sum(); t[i] = 0.98 * es[i] * air.sum()
        i550 = list(C.LAMs).index(550.0)
        out[f'{cav[0]:.1f}/{cav[1]:.1f}|{lens}'] = dict(b4060_550=float(b[i550]), b4060_spec=float(b @ C.s), total_spec=float(t @ C.s), eta_sub_550=float(es[i550]))
for cav in ('61.3/97.5', '74.6/81.7', '68.5/114.4'):
    h, f = out[f'{cav}|hemi_40-60'], out[f'{cav}|free_40-60']
    print(cav, 'free/hemi 40-60: 550 nm %.3f  spectral %.3f' % (f['b4060_550'] / h['b4060_550'], f['b4060_spec'] / h['b4060_spec']), ' hemi spec %.4f free spec %.4f' % (h['b4060_spec'], f['b4060_spec']))
json.dump(out, open('data/cav_lens_cross.json', 'w'), indent=1)
